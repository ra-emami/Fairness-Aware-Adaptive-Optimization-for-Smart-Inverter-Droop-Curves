# =============================================================================
# Scenario 6 — Variable (adaptive) Volt-Var breakpoints + priority-based fairness
#
# Multi-period enhanced IVACOPF for optimal smart-inverter Q-V droop settings
# on the 33-bus distribution feeder: 24-hour horizon at 15-minute resolution
# (96 steps). The four interior breakpoint voltages of each inverter's Q-V
# droop curve are decision variables, re-optimized once per hour and held
# constant across the four 15-minute intervals of that hour (Eqs. 25-27).
# PV curtailment follows the priority-based rule: units with the highest
# available capacity are curtailed first, enforced through the big-M
# reformulation of the conditional curtailment rules (Eqs. 37, 41-45).
#
# Reference:
#   R. Emami Mirak and A. Inaolaji, "Adaptive and fair optimization of smart
#   inverter droop curves in distribution grids," Electric Power Systems
#   Research, vol. 262, 2027, Art. no. 113613.
#   https://doi.org/10.1016/j.epsr.2026.113613
#   (Presented at the 24th Power Systems Computation Conference, PSCC 2026.)
# =============================================================================

using JuMP
using Gurobi
using DataFrames, CSV, Plots, Logging

Logging.disable_logging(Logging.Warn)

# ------------------------------------------------------------------ plotting
function Pdg_plot()
    hour_labels = [lpad(h, 2, "0") for h in 0:23]
    hour_indices = 1:4:96

    plt = plot(
        label="DG Active Power - IVACOPF",
        xlabel="Hour",
        ylabel="Active Power (p.u.)",
        title="DG Active Power - IVACOPF",
        legend=:topleft,
        legend_background_color = :transparent,
        legend_foreground_color = RGBA(0, 0, 0, 0.3),
        xticks=(hour_indices, hour_labels),
        size=(800, 600)
    )

    plP = reshape(permutedims(Pdg2, [1, 3, 2]), length(DG_SET), length(QUARTER_SET)*length(HOUR_SET))
    for i in 1:length(DG_SET)
        plP_max = reshape(Pdg_max_vary[DG_SET[i]]', length(QUARTER_SET)*length(HOUR_SET), 1)
        d = DG_SET[i]
        plot!(1:96, plP[i,:], label="DG $d", linewidth=2)
        plot!(1:96, plP_max, linestyle=:dash, linewidth=0.8, label=false)
    end

    for i in 1:length(DG_SET)
        plP_max = reshape(Pdg_max_vary[DG_SET[i]]', length(QUARTER_SET)*length(HOUR_SET), 1)
        d = DG_SET[i]
        plot!(1:96, plP_max-plP[i,:], label="DG $d Curtailment", linewidth=1)
    end

    plot!(1:96, reshape(permutedims(Pgen2, [1, 3, 2]), length(QUARTER_SET)*length(HOUR_SET)),
        label="Subs. Power",
        color=:gray40,
        linewidth=1
        )

    plot!(1:96, sum(reshape(permutedims(Pload, [1, 3, 2]), length(BUS_SET),length(QUARTER_SET)*length(HOUR_SET)), dims=1)',
        label="Load",
        color=:green,
        linewidth=1
        )

    plot!(1:96, reshape(sum(permutedims(Ploss2, [1, 3, 2]), dims=1), length(QUARTER_SET)*length(HOUR_SET)),
        label="Ploss",
        color=:red,
        linewidth=1
        )
    plot!(1:96, reshape(sum(permutedims(Ploss2, [1, 3, 2]), dims=1), length(QUARTER_SET)*length(HOUR_SET))
        + sum(reshape(permutedims(Pload, [1, 3, 2]), length(BUS_SET),length(QUARTER_SET)*length(HOUR_SET)), dims=1)',
        label="Total Gen",
        color=:purple,
        linestyle=:dash,
        linewidth=1
        )

    display(plt)
end

function vp_plot()
    plt = plot(
        xlabel="Bus",
        ylabel="Voltage (p.u.)",
        title="Voltage Profile - IVACOPF",
        xticks=(BUS_SET, BUS_SET),
        ylim=(0.89, 1.05),
        yticks=(0.9:0.02:1.05),
        legend=:bottomright
    )
    plot!(BUS_SET, V_max, label="Max Voltage", linewidth=2, color="darkorange2")
    plot!(BUS_SET, V_max, fillrange=V_min, fillalpha=0.2, label=false, linecolor=:transparent, color="dodgerblue")
    plot!(BUS_SET, V_min, label="Min Voltage", linewidth=2, color="dodgerblue4")
    annotate!(1, 1.05,
        text(voltage_range,
        :left,
        :top,
        6,
        color=:red))

    hline!([1], linestyle=:dash, color=:gray40, label=false, linewidth=0.5)
    hline!([V_limit[1]], linestyle=:dash, color=:gray40, label=false, linewidth=0.5)
    hline!([V_limit[2]], linestyle=:dash, color=:gray40, label=false, linewidth=0.5)

    display(plt)
end

# ---------------------------------------------------------------------- sets
HOUR_SET = 1:24
QUARTER_SET = 1:4
HQ_SET = vcat([collect((h*4-4).+QUARTER_SET) for h in HOUR_SET]...)

BUS_SET = 1:33
slack = 1
Vnom = 1

# ---------------------------------------------------------------- input data
DATA_DIR = joinpath(dirname(@__DIR__), "data")
System_Data = Matrix(DataFrame(CSV.File(joinpath(DATA_DIR, "33_bus_system_data.csv"), skipto=1)))
solar_profile = reshape(DataFrame(CSV.File(joinpath(DATA_DIR, "solar_profile_15min.csv")))[:,2],length(QUARTER_SET), length(HOUR_SET))'/100 # p.u.

# Base values
Sbase = parse(Float64, System_Data[2,9])*1e3 # in VA
Vbase = parse(Float64, System_Data[2,8])*1e3 # in V
Zbase = Vbase^2/Sbase # in Ω

Pload_std = reshape(parse.(Float64, String.(collect(skipmissing(System_Data[BUS_SET.+1, 4])))), :, 1) * 1e3 # in W
Qload_std = reshape(parse.(Float64, String.(collect(skipmissing(System_Data[BUS_SET.+1, 5])))), :, 1) * 1e3 # in VAr
Zdata = reshape(parse.(Float64, String.(collect(skipmissing(hcat(System_Data[BUS_SET[2:end].+1, 2:3], System_Data[BUS_SET[2:end].+1, 6:7]))))), :, 4) # in Ω

# 15-minute load variation profiles (% of peak) and bus -> load-class mapping:
# class 1 = industrial (buses 26-33), 2 = commercial (19-25), 3 = residential (2-18)
load_variation = DataFrame(CSV.File(joinpath(DATA_DIR, "load_variation_15min.csv"), skipto=1))
load_category = Dict(1 => "Industrial", 2 => "Commercial", 3 => "Residential")
class = Dict(parse.(Int, load_variation[2:34,10]) .=> parse.(Int, load_variation[2:34,11]))
load_percent_matrix = Matrix{Float64}(tryparse.(Float64, load_variation[HQ_SET.+1, 2:7])) *1e-2

load_percent_p = load_percent_matrix[:,1:3]
load_percent_q = load_percent_matrix[:,4:6]

Pload = zeros(length(BUS_SET),24,4)
Qload = zeros(length(BUS_SET),24,4)

for bus in BUS_SET[2:end]
    for h in 1:length(HOUR_SET)
        for q in 1:length(QUARTER_SET)
            Pload[bus, HOUR_SET[h], QUARTER_SET[q]] = Pload_std[bus] * load_percent_p[q+(h-1)*length(QUARTER_SET), class[bus]]/Sbase
            Qload[bus, HOUR_SET[h], QUARTER_SET[q]] = Qload_std[bus] * load_percent_q[q+(h-1)*length(QUARTER_SET), class[bus]]/Sbase
        end
    end
end

R = Dict()
X = Dict()

for k in 1:size(Zdata,1)
    R[Int(Zdata[k,1]),Int(Zdata[k,2])] = Zdata[k,3]/Zbase
    X[Int(Zdata[k,1]),Int(Zdata[k,2])] = Zdata[k,4]/Zbase
end

BRANCH_SET = [(Int(Zdata[k,1]), Int(Zdata[k,2])) for k in 1:size(Zdata,1)]
Bi_BRANCH_SET = vcat(BRANCH_SET, [(b,a) for (a,b) in BRANCH_SET])

# ------------------------------------------------------- PV / smart inverters
# Three PV units at buses 7, 18, 33 rated 2.2, 1.0, 1.5 MW; inverter
# apparent-power ratings oversized by 10% (S = 1.1 P).
DG_SET = [7 18 33]
DG_Psize = [2.2e6 1e6 1.5e6] # Max Pdg in W
DG_Ssize = 1.1 .* DG_Psize   # Max Sdg in VA
DG_SIZE =  [DG_Psize;
            DG_Ssize]

dg_number = length(DG_SET)
Pdg_max = Dict(DG_SET[i] => DG_SIZE[1,i]/Sbase for i in 1:dg_number) # p.u.
Sdg_max = Dict(DG_SET[i] => DG_SIZE[2,i]/Sbase for i in 1:dg_number) # p.u.
Qdg_max = Sdg_max # p.u.
NON_DG_SET = setdiff(BUS_SET,DG_SET,slack)

Pdg_max_vary = Dict(DG => Pdg_max[DG] * solar_profile for DG in DG_SET)

# Droop ordinates (Eq. 20): (+Qmax, +Qmax, 0, 0, -Qmax, -Qmax) with Qmax = Smax,
# since reactive capability remains fully available independently of active
# generation; the joint (P, Q) feasibility is enforced by the polygonal
# capability constraints (Eq. 17).
Qcap = Dict(DG => [Sdg_max[DG] for h in HOUR_SET] for DG in DG_SET)

qG_values = Dict(dg =>
    vcat(Qcap[dg]', Qcap[dg]', zeros(1, length(HOUR_SET)), zeros(1, length(HOUR_SET)), -Qcap[dg]', -Qcap[dg]' ) for dg in vec(DG_SET))

V_limit = [0.95 1.05] # bus voltage bounds in p.u. (Eq. 5)
global Vbp_limit = [0.80 1.20] # fixed extreme breakpoints V1, V6 in p.u. (Eq. 27)

# Linearization points (previous iterate); flat start at nominal voltage
global v_r_pr = ones(length(BUS_SET), 24, 4)
global v_im_pr = zeros(length(BUS_SET), 24, 4)
global Ibs_r_pr = zeros(length(BUS_SET), 24, 4)
global Ibs_im_pr = zeros(length(BUS_SET), 24, 4)
global Ibr_r_pr = Dict((branch, hour, quarter) => 0.0 for branch in BRANCH_SET, hour in 1:24, quarter in 1:4)
global Ibr_im_pr = Dict((branch, hour, quarter) => 0.0 for branch in BRANCH_SET, hour in 1:24, quarter in 1:4)

# ------------------------------------------------- iterative enhanced IVACOPF
@time begin
max_iter = 10       # iteration cap; convergence typically occurs in 3-4 iterations
conv_it = max_iter
converged = false
dif_h = zeros(max_iter)
tol = 1e-6          # convergence tolerance on the MLE (p.u.)
t_opt = zeros(max_iter)
for it in 1:max_iter
    global t_opt
    global v_r_pr
    global v_im_pr
    global Ibs_r_pr
    global Ibs_im_pr
    global Ibr_r_pr
    global Ibr_im_pr
    global Vbp_limit

    println("=============================== Iteration $it ===============================")

    global OPF = Model(Gurobi.Optimizer)
    Time_Limit = 3600 # in seconds; safety cap, never binding in the reported runs
    set_optimizer_attribute(OPF, "TimeLimit", Time_Limit)
    set_optimizer_attribute(OPF, "MIPGap", 0.001)  # 0.1% relative MIP optimality gap
    set_optimizer_attribute(OPF, "Threads", 0)     # all cores (Gurobi default)

    # Decision variables
    @variable(OPF, 0 <= Pgen[slack, HOUR_SET, QUARTER_SET])  # Active power generation at slack bus
    @variable(OPF, Qgen[slack, HOUR_SET, QUARTER_SET])  # Reactive power generation at slack bus
    @variable(OPF, V_limit[1] <= v[BUS_SET, HOUR_SET, QUARTER_SET] <= V_limit[2], start = Vnom) # voltage magnitude at bus i
    @variable(OPF, Δv[BUS_SET, HOUR_SET, QUARTER_SET])  # voltage deviation |v - Vnom|
    @variable(OPF, v_r[BUS_SET, HOUR_SET, QUARTER_SET], start = Vnom) # real part of the complex voltage
    @variable(OPF, v_im[BUS_SET, HOUR_SET, QUARTER_SET], start = 0) # imaginary part of the complex voltage
    @variable(OPF, Ibr_r[BRANCH_SET, HOUR_SET, QUARTER_SET])   # branch current, real part
    @variable(OPF, Ibr_im[BRANCH_SET, HOUR_SET, QUARTER_SET])  # branch current, imaginary part
    @variable(OPF, Psnd[Bi_BRANCH_SET, HOUR_SET, QUARTER_SET]) # sending-end active power (Eq. 9)
    @variable(OPF, Qsnd[Bi_BRANCH_SET, HOUR_SET, QUARTER_SET]) # sending-end reactive power (Eq. 10)
    @variable(OPF, 0 <= Ploss[BRANCH_SET, HOUR_SET, QUARTER_SET]) # branch active-power loss
    @variable(OPF, 0 <= Qloss[BRANCH_SET, HOUR_SET, QUARTER_SET]) # branch reactive-power loss
    @variable(OPF, Ibs_r[BUS_SET, HOUR_SET, QUARTER_SET])   # bus injection current, real part
    @variable(OPF, Ibs_im[BUS_SET, HOUR_SET, QUARTER_SET])  # bus injection current, imaginary part
    @variable(OPF, 0 <= Pdg[DG_SET, HOUR_SET, QUARTER_SET])  # PV active power generation
    @variable(OPF, Qdg[DG_SET, HOUR_SET, QUARTER_SET])  # PV reactive power generation
    @variable(OPF, 0 <= fairness_index[HOUR_SET, QUARTER_SET])  # ψ_t (Eq. 37)

    # ---------------- Volt-Var droop: Lambda method with variable breakpoints
    @variable(OPF, λ[1:6, DG_SET, HOUR_SET, QUARTER_SET] >= 0)  # interpolation weights (Eq. 21)
    @variable(OPF, δ[1:5, DG_SET, HOUR_SET, QUARTER_SET], Bin)  # SOS2 adjacency indicators (Eq. 24)
    @variable(OPF, V_break[1:6, DG_SET, HOUR_SET])  # breakpoint voltages, re-optimized hourly (Eq. 25)

    # Extreme breakpoints fixed (Eq. 27)
    @constraint(OPF, [d in DG_SET, h in HOUR_SET], V_break[1, d, h] == Vbp_limit[1])
    @constraint(OPF, [d in DG_SET, h in HOUR_SET], V_break[6, d, h] == Vbp_limit[2])
    # Interior breakpoint bounds (Eq. 26)
    Vb_lim = [
        Vbp_limit[1] Vbp_limit[1]; # V_breakpoint 1 limit
        0.82 0.98; # V_breakpoint 2 limit | V2 <= {V3 − 0.02 VN}
        0.97 1; # V_breakpoint 3 limit
        1 1.03; # V_breakpoint 4 limit
        1.02 1.18; # V_breakpoint 5 limit | {V4 + 0.02 VN} <= V5
        Vbp_limit[2] Vbp_limit[2] # V_breakpoint 6 limit
    ]

    @constraint(OPF, [i = 2, d in DG_SET, h in HOUR_SET], V_break[i, d, h] <= V_break[i+1, d, h] - 0.02)
    @constraint(OPF, [i = 5, d in DG_SET, h in HOUR_SET], V_break[i-1, d, h] + 0.02 <= V_break[i, d, h])
    margin = 0
    @constraint(OPF, [i = 1:5, d in DG_SET, h in HOUR_SET], V_break[i, d, h] + margin <= V_break[i+1, d, h])
    @constraint(OPF, [i = 2:5, d in DG_SET, h in HOUR_SET], Vb_lim[i,1] <= V_break[i, d, h])
    @constraint(OPF, [i = 2:5, d in DG_SET, h in HOUR_SET], V_break[i, d, h] <= Vb_lim[i,2])

    # PCC voltage interpolation (Eq. 22); bilinear in λ and V_break
    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], v[d,h,m] == sum(λ[i, d, h, m] * V_break[i, d, h] for i in 1:6))

    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], sum(λ[i, d, h, m] for i in 1:6) == 1)

    # SOS2 structure via binary adjacency indicators (Eq. 24)
    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], λ[1, d, h, m] <= δ[1, d, h, m])
    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], λ[2, d, h, m] <= δ[1, d, h, m] + δ[2, d, h, m])
    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], λ[3, d, h, m] <= δ[2, d, h, m] + δ[3, d, h, m])
    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], λ[4, d, h, m] <= δ[3, d, h, m] + δ[4, d, h, m])
    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], λ[5, d, h, m] <= δ[4, d, h, m] + δ[5, d, h, m])
    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], λ[6, d, h, m] <= δ[5, d, h, m])
    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], sum(δ[:, d, h, m]) == 1)  # only one pair of adjacent segments can be active

    # Reactive output from the droop curve (Eq. 23)
    @constraint(OPF, Reactive_DG[d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            Qdg[d, h, m] == sum(λ[i, d, h, m] * qG_values[d][i, h] for i in 1:6))

    # -------------- SI capability region: 32-sided polygonal semicircle (Eq. 17)
    k = 16  # number of constraint pairs
    γ = π / k  # angle step

    for l in 1:k
        θ = l * γ
        @constraint(OPF, [d in 1:dg_number, h in HOUR_SET, m in QUARTER_SET], cos(θ) * Pdg[DG_SET[d], h, m] + sin(θ) * Qdg[DG_SET[d], h, m] <= Sdg_max[DG_SET[d]])
        @constraint(OPF, [d in 1:dg_number, h in HOUR_SET, m in QUARTER_SET], cos(θ) * Pdg[DG_SET[d], h, m] + sin(θ) * Qdg[DG_SET[d], h, m] >= -Sdg_max[DG_SET[d]])
    end

    # Active/reactive output limits (Eq. 18)
    @constraint(OPF, DG_max1[DG in DG_SET, h in HOUR_SET, m in QUARTER_SET], Pdg[DG, h, m] <= Pdg_max_vary[DG][h,m])
    @constraint(OPF, [i in DG_SET, h in HOUR_SET, m in QUARTER_SET], Pdg[i,h,m] <= Pdg_max[i])
    @constraint(OPF, [i in DG_SET, h in HOUR_SET, m in QUARTER_SET], Qdg[i,h,m] <= Qdg_max[i])
    @constraint(OPF, [i in DG_SET, h in HOUR_SET, m in QUARTER_SET], -Qdg_max[i] <= Qdg[i,h,m])

    # ------------------------------------------- IVACOPF network constraints
    # Slack bus reference
    @constraint(OPF, [i in slack, h in HOUR_SET, m in QUARTER_SET], v_r[i,h,m] == Vnom)
    @constraint(OPF, [i in slack, h in HOUR_SET, m in QUARTER_SET], v_im[i,h,m] == 0)

    # Branch voltage-current relations (Eqs. 1-2)
    @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], v_r[i,h,m]-v_r[j,h,m] == R[i, j]*Ibr_r[(i, j),h,m] - X[i, j]*Ibr_im[(i, j),h,m] )
    @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], v_im[i,h,m]-v_im[j,h,m] == R[i, j]*Ibr_im[(i, j),h,m] + X[i, j]*Ibr_r[(i, j),h,m] )

    # Bus current injections (Eqs. 3-4)
    @constraint(OPF, [bus in BUS_SET, h in HOUR_SET, m in QUARTER_SET], Ibs_r[bus,h,m] ==  sum(Ibr_r[(bus, j),h,m] for (i, j) in BRANCH_SET if i == bus)
                                                    - sum(Ibr_r[(i, bus),h,m] for (i, j) in BRANCH_SET if j == bus))

    @constraint(OPF, [bus in BUS_SET, h in HOUR_SET, m in QUARTER_SET], Ibs_im[bus,h,m] ==  sum(Ibr_im[(bus, j),h,m] for (i, j) in BRANCH_SET if i == bus)
                                                    - sum(Ibr_im[(i, bus),h,m] for (i, j) in BRANCH_SET if j == bus))

    # Linearized power-balance constraints (Eqs. 7-8), first-order Taylor
    # around the previous iterate (v_r_pr, v_im_pr, Ibs_r_pr, Ibs_im_pr)
    @constraint(OPF, [i in slack, h in HOUR_SET, m in QUARTER_SET], Pgen[i,h,m] == v_r_pr[i,h,m]*Ibs_r[i,h,m] + v_im_pr[i,h,m]*Ibs_im[i,h,m]
                         + Ibs_r_pr[i,h,m]*v_r[i,h,m] + Ibs_im_pr[i,h,m]*v_im[i,h,m] - v_r_pr[i,h,m]*Ibs_r_pr[i,h,m] - v_im_pr[i,h,m]*Ibs_im_pr[i,h,m])

    @constraint(OPF, [i in slack, h in HOUR_SET, m in QUARTER_SET], Qgen[i,h,m] == v_im_pr[i,h,m]*Ibs_r[i,h,m] - v_r_pr[i,h,m]*Ibs_im[i,h,m]
                        + Ibs_r_pr[i,h,m]*v_im[i,h,m] - Ibs_im_pr[i,h,m]*v_r[i,h,m] - v_im_pr[i,h,m]*Ibs_r_pr[i,h,m] + v_r_pr[i,h,m]*Ibs_im_pr[i,h,m])

    @constraint(OPF, [i in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET], -Pload[i,h,m] == v_r_pr[i,h,m]*Ibs_r[i,h,m] + v_im_pr[i,h,m]*Ibs_im[i,h,m]
                        + Ibs_r_pr[i,h,m]*v_r[i,h,m] + Ibs_im_pr[i,h,m]*v_im[i,h,m] - v_r_pr[i,h,m]*Ibs_r_pr[i,h,m] - v_im_pr[i,h,m]*Ibs_im_pr[i,h,m])

    @constraint(OPF, [i in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET], -Qload[i,h,m] == v_im_pr[i,h,m]*Ibs_r[i,h,m] - v_r_pr[i,h,m]*Ibs_im[i,h,m]
                        + Ibs_r_pr[i,h,m]*v_im[i,h,m] - Ibs_im_pr[i,h,m]*v_r[i,h,m] - v_im_pr[i,h,m]*Ibs_r_pr[i,h,m] + v_r_pr[i,h,m]*Ibs_im_pr[i,h,m])

    @constraint(OPF, [i in DG_SET, h in HOUR_SET, m in QUARTER_SET], -Pload[i,h,m] + Pdg[i,h,m] == v_r_pr[i,h,m]*Ibs_r[i,h,m] + v_im_pr[i,h,m]*Ibs_im[i,h,m]
                        + Ibs_r_pr[i,h,m]*v_r[i,h,m] + Ibs_im_pr[i,h,m]*v_im[i,h,m] - v_r_pr[i,h,m]*Ibs_r_pr[i,h,m] - v_im_pr[i,h,m]*Ibs_im_pr[i,h,m])

    @constraint(OPF, [i in DG_SET, h in HOUR_SET, m in QUARTER_SET], -Qload[i,h,m] + Qdg[i,h,m] == v_im_pr[i,h,m]*Ibs_r[i,h,m] - v_r_pr[i,h,m]*Ibs_im[i,h,m]
                        + Ibs_r_pr[i,h,m]*v_im[i,h,m] - Ibs_im_pr[i,h,m]*v_r[i,h,m] - v_im_pr[i,h,m]*Ibs_r_pr[i,h,m] + v_r_pr[i,h,m]*Ibs_im_pr[i,h,m])

    # Linearized branch losses (Eqs. 15-16)
    @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Ploss[(i, j), h, m] == R[i, j] * (
                                                                        2 * Ibr_r_pr[(i, j), h, m]  * Ibr_r[(i, j), h, m] - Ibr_r_pr[(i, j), h, m]^2
                                                                        + 2 * Ibr_im_pr[(i, j), h, m] * Ibr_im[(i, j), h, m] - Ibr_im_pr[(i, j), h, m]^2))

    @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Qloss[(i, j), h, m] == X[i, j] * (
                                                                2 * Ibr_r_pr[(i, j), h, m]  * Ibr_r[(i, j), h, m] - Ibr_r_pr[(i, j), h, m]^2
                                                                + 2 * Ibr_im_pr[(i, j), h, m] * Ibr_im[(i, j), h, m] - Ibr_im_pr[(i, j), h, m]^2))

    # Sending-end power balance and consistency constraints (Eqs. 9-12)
    @constraint(OPF, [bus in slack, h in HOUR_SET, m in QUARTER_SET], Pgen[bus,h,m] == sum(Psnd[(i, j),h,m] for (i, j) in Bi_BRANCH_SET if i == bus) + Pload[bus,h,m] )
    @constraint(OPF, [bus in DG_SET, h in HOUR_SET, m in QUARTER_SET], Pdg[bus,h,m] == sum(Psnd[(i, j),h,m] for (i, j) in Bi_BRANCH_SET if i == bus) + Pload[bus,h,m] )
    @constraint(OPF, [bus in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET], 0 == sum(Psnd[(i, j),h,m] for (i, j) in Bi_BRANCH_SET if i == bus) + Pload[bus,h,m] )

    @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Psnd[(i, j),h,m] == Ploss[(i, j),h,m] - Psnd[(j, i),h,m])

    @constraint(OPF, [bus in slack, h in HOUR_SET, m in QUARTER_SET], Qgen[bus,h,m] == sum(Qsnd[(i, j),h,m] for (i, j) in Bi_BRANCH_SET if i == bus) + Qload[bus,h,m] )
    @constraint(OPF, [bus in DG_SET, h in HOUR_SET, m in QUARTER_SET], Qdg[bus,h,m] == sum(Qsnd[(i, j),h,m] for (i, j) in Bi_BRANCH_SET if i == bus) + Qload[bus,h,m] )
    @constraint(OPF, [bus in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET], 0 == sum(Qsnd[(i, j),h,m] for (i, j) in Bi_BRANCH_SET if i == bus) + Qload[bus,h,m] )

    @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Qsnd[(i, j),h,m] == Qloss[(i, j),h,m] - Qsnd[(j, i),h,m])

    # Linearized voltage magnitude (Eq. 6)
    @constraint(OPF, [i in BUS_SET, h in HOUR_SET, m in QUARTER_SET], v[i,h,m] == ( (v_r_pr[i,h,m] / sqrt(v_r_pr[i,h,m]^2 + v_im_pr[i,h,m]^2)) * v_r[i,h,m] )
                        + ( (v_im_pr[i,h,m] / sqrt(v_r_pr[i,h,m]^2 + v_im_pr[i,h,m]^2)) * v_im[i,h,m] ))

    # Voltage deviation |v - Vnom| for the AVD objective term (Eq. 34)
    @constraint(OPF, [i in BUS_SET, h in HOUR_SET, m in QUARTER_SET], v[i,h,m]-Vnom <= Δv[i,h,m])
    @constraint(OPF, [i in BUS_SET, h in HOUR_SET, m in QUARTER_SET], -v[i,h,m]+Vnom <= Δv[i,h,m])

    # ---------------- fairness: priority-based curtailment (Eqs. 37, 41-45)
    # Units with the highest available capacity are curtailed first. The
    # conditional rules of Eqs. 39-40 are reformulated with the binaries
    # fb1, fb2, fb3 and the disjunctive constant M (Eqs. 41-45):
    #   fb1 = 1 <=> Pdg[DG] <= Pdg_max_vary[d]
    #   fb2 = 1 <=> Pdg_max_vary[d] <= Pdg_max_vary[DG]
    #   fb3 = 1 <=> Pdg[DG] <= Pdg[d]
    fM = maximum(values(Pdg_max))   # M = max Pmax, the tightest valid choice
    @variable(OPF, fb1[DG in DG_SET, d in DG_SET, h in HOUR_SET, m in QUARTER_SET], Bin)
    @variable(OPF, fb2[DG in DG_SET, d in DG_SET, h in HOUR_SET, m in QUARTER_SET], Bin)
    @variable(OPF, fb3[DG in DG_SET, d in DG_SET, h in HOUR_SET, m in QUARTER_SET], Bin)
    for DG in DG_SET
        @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
                        Pdg[DG, h, m] <= Pdg[d, h, m] + (2-fb1[DG, d, h, m]-fb2[DG, d, h, m]) * fM)
        @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
                        Pdg[DG, h, m] >= Pdg[d, h, m] - (2-fb1[DG, d, h, m]-fb2[DG, d, h, m]) * fM)

        @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
                        Pdg[DG, h, m] <= Pdg_max_vary[DG][h,m] + (1+fb2[DG, d, h, m]-fb3[DG, d, h, m]) * fM)
        @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
                        Pdg[DG, h, m] >= Pdg_max_vary[DG][h,m] - (1+fb2[DG, d, h, m]-fb3[DG, d, h, m]) * fM)

        @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
                        Pdg[DG, h, m] <= Pdg_max_vary[d][h,m] + (1-fb1[DG, d, h, m]) * fM)
        @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
                        - fb1[DG, d, h, m] * fM + Pdg_max_vary[d][h,m] <= Pdg[DG, h, m])

        @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
                        Pdg_max_vary[d][h,m] <= Pdg_max_vary[DG][h,m] + (1-fb2[DG, d, h, m]) * fM)
        @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
                        - fb2[DG, d, h, m] * fM + Pdg_max_vary[DG][h,m] <= Pdg_max_vary[d][h,m])

        @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
                                        Pdg[DG, h, m] <= Pdg[d, h, m] + (1-fb3[DG, d, h, m]) * fM)
        @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
                                        - fb3[DG, d, h, m] * fM + Pdg[d, h, m] <= Pdg[DG, h, m])
    end

    # Egalitarian curtailment-margin constraint (Eq. 37)
    @constraint(OPF, Fairness[DG in DG_SET, h in HOUR_SET, m in QUARTER_SET], fairness_index[h,m] >=
                                                                    Pdg_max_vary[DG][h,m]-Pdg[DG, h, m])

    # ------------------------------------------------------- objective (Eq. 33)
    # Weighted sum of the normalized fairness term Ψ and the normalized
    # voltage-deviation term. In this implementation α weights Ψ and β = 1 - α
    # weights the voltage-deviation term, matching the (α, β) labeling of the
    # Pareto front in Fig. 9 of the paper. f1 and f2 are the normalization
    # constants η obtained from the respective single-objective optima; the
    # 1/(|B|·|T|) factor of the AVD (Eq. 34) is absorbed into f2.
    global scenario = 6
    α = 0.5
    β = 1-α

    f1 = 1.0323678    # η for the priority-based fairness term (variable breakpoints)
    f2 = 26.8814161   # η for the voltage-deviation term (variable breakpoints)

    @objective(OPF, Min, α*sum(fairness_index[h,m] for h in HOUR_SET, m in QUARTER_SET)/f1 +
                                β*sum(Δv[i,h,m] for i in BUS_SET, h in HOUR_SET, m in QUARTER_SET)/f2)

    println("----------------------------------------------------------------------------")
    println("Number of Variables: ", num_variables(OPF))
    println("Total Constraints: ", sum(length(all_constraints(OPF, F, S)) for (F, S) in list_of_constraint_types(OPF)))
    println("Number of Binary Variables: ", count(is_binary, all_variables(OPF)))
    println("Number of Integer Variables: ", count(is_integer, all_variables(OPF)))
    println("----------------------------------------------------------------------------\n")

    t_opt[it] = @elapsed status = optimize!(OPF)
    println("Iteration $(it) - optimization time: ", t_opt[it], " s")

    # Maximum linearization error (MLE): largest absolute difference between
    # the bilinear branch loss and its linearized counterpart (Eqs. 13 vs 15)
    max_diff = maximum([abs(value((v_r[i,h,m]-v_r[j,h,m])*Ibr_r[(i, j),h,m] + (v_im[i,h,m]-v_im[j,h,m])*Ibr_im[(i, j),h,m]
                - Ploss[(i, j),h,m])) for (i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET])
    dif_h[it] = max_diff

    # Update the linearization points with the current solution
    [v_r_pr[b, h, q] = value(v_r[b, h, q]) for b in BUS_SET, h in HOUR_SET, q in QUARTER_SET]
    [v_im_pr[b, h, q] = value(v_im[b, h, q]) for b in BUS_SET, h in HOUR_SET, q in QUARTER_SET]
    [Ibs_r_pr[b, h, q] = value(Ibs_r[b, h, q]) for b in BUS_SET, h in HOUR_SET, q in QUARTER_SET]
    [Ibs_im_pr[b, h, q] = value(Ibs_im[b, h, q]) for b in BUS_SET, h in HOUR_SET, q in QUARTER_SET]
    [Ibr_r_pr[(br, h, q)] = value(Ibr_r[br, h, q]) for br in BRANCH_SET, h in HOUR_SET, q in QUARTER_SET]
    [Ibr_im_pr[(br, h, q)] = value(Ibr_im[br, h, q]) for br in BRANCH_SET, h in HOUR_SET, q in QUARTER_SET]

    # ---------------------------------------------------- iteration reporting
    println("Iteration $it - Substation P (kWh): ", round.(sum(Array(value.(Pgen)), dims=(1,3))[:,:,1]*Sbase*1e-3, digits=3))
    println("Iteration $it - Substation Q (kVArh): ", round.(sum(Array(value.(Qgen)), dims=(1,3))[:,:,1]*Sbase*1e-3, digits=3))
    println("Iteration $it - PV P (kWh): ", round.(sum(Array(value.(Pdg)), dims=(1,3))[:,:,1]*Sbase*1e-3, digits=3))
    println("Iteration $it - PV Q (kVArh): ", round.(sum(Array(value.(Qdg)), dims=(1,3))[:,:,1]*Sbase*1e-3, digits=3))
    println("Iteration $it - Total P loss (kW): ", round.(sum(value.(Ploss[(i,j),h,m] for (i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET), dims=(1,3))[:,:,1]*Sbase*1e-3, digits=6))
    println("Iteration $it - Total Q loss (kVAr): ", round.(sum(value.(Qloss[(i,j),h,m] for (i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET), dims=(1,3))[:,:,1]*Sbase*1e-3, digits=6))
    println("Iteration $it - Total P load (kWh): ", round.(sum(Pload, dims=(1,3))[:,:,1]*Sbase*1e-3, digits=3))
    println("Iteration $it - Total Q load (kVArh): ", round.(sum(Qload, dims=(1,3))[:,:,1]*Sbase*1e-3, digits=3))

    # PV curtailment metrics (Eq. 35)
    global PVC = reshape([Pdg_max_vary[DG][h,q] - value(Pdg[DG, h, q]) for DG in DG_SET, h in HOUR_SET, q in QUARTER_SET], length(DG_SET), length(HOUR_SET), length(QUARTER_SET))
    global PVC_sum = sum(PVC[d,h,m] for d in 1:dg_number, h in HOUR_SET, m in QUARTER_SET)
    global PVC_sep = Dict(DG => [sum(PVC[findfirst(isequal(DG), vec(DG_SET)),h,m] for h in HOUR_SET, m in QUARTER_SET)] for DG in DG_SET)

    PVC_pct = 100*PVC_sum/sum(Pdg_max_vary[DG][h,m] for DG in DG_SET, h in HOUR_SET, m in QUARTER_SET)
    println("PV curtailment: ", round(PVC_sum/4, digits=7), " p.u.h = ",
            round(Sbase*PVC_sum/4e6, digits=4), " MWh  |  PVC = ", round(PVC_pct, digits=4), " %")
    println("Curtailment by PV unit (MWh): ", [[DG, round.(Sbase*PVC_sep[DG]/4e6, digits=4)] for DG in DG_SET ])

    println("alpha: ", α, "  beta: ", β)
    fairness_prnt = "Fairness index: " *string(round(value(sum(fairness_index[h,m] for h in HOUR_SET, m in QUARTER_SET)), digits=7))
    println(fairness_prnt)

    # Voltage metrics
    V = [value(v[i,h,m]) for i in sort(BUS_SET), h in HOUR_SET, m in QUARTER_SET]
    global MaxΔV = maximum(abs.(V .- Vnom))
    println("Max |V - Vnom|: ", round(MaxΔV, digits=7), " p.u.")
    AVD = sum(abs.(V .- Vnom))/(length(BUS_SET)*length(HQ_SET))
    println("AVD (Eq. 34): ", round(AVD*1e3, digits=4), " x 1e-3 p.u.")
    println("Total |V - Vnom|: ", round(sum(abs.(V .- Vnom)), digits=7), " p.u.")
    Ploss_prnt = "Total Ploss: " *string(round(value(sum(Ploss[(i, j), h, m] for (i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET)), digits=5))
    println(Ploss_prnt)
    global V_max = maximum(V, dims=[2,3])[:, 1, 1]  # max over hours and quarters
    global V_min = minimum(V, dims=[2,3])[:, 1, 1]  # min over hours and quarters

    println("MLE: ", max_diff)
    if max_diff < tol
        global conv_it = it
        global converged = true
        global Vbp = value.(V_break[:,:,:])
        global V_ac_l2 = [value(v[i, h, m]) for i in BUS_SET, h in HOUR_SET, m in QUARTER_SET]
        global Pgen2 = reshape([ value(Pgen[i,h,m]) for i in (slack), h in HOUR_SET, m in QUARTER_SET],
                                                    length(slack), length(HOUR_SET), length(QUARTER_SET))
        global Qgen2 = reshape([ value(Qgen[i,h,m]) for i in (slack), h in HOUR_SET, m in QUARTER_SET],
                                                    length(slack), length(HOUR_SET), length(QUARTER_SET))
        global Pdg2 = reshape([ value(Pdg[i,h,m]) for i in (DG_SET), h in HOUR_SET, m in QUARTER_SET],
                                                    length(DG_SET), length(HOUR_SET), length(QUARTER_SET))
        global Ploss2 = reshape([ value(Ploss[(i,j),h,m]) for (i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
                                                    length(BRANCH_SET), length(HOUR_SET), length(QUARTER_SET))
        global Qloss2 = reshape([ value(Qloss[(i,j),h,m]) for (i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
                                                    length(BRANCH_SET), length(HOUR_SET), length(QUARTER_SET))
        global Qdg2 = reshape([ value(Qdg[i,h,m]) for i in (DG_SET), h in HOUR_SET, m in QUARTER_SET],
                                            length(DG_SET), length(HOUR_SET), length(QUARTER_SET))
        println("Converged at iteration $it (MLE = $max_diff)")
        println("Scenario $(scenario) - total optimization time: $(sum(t_opt)) s")
        break
    end

    global V_ac_l2 = [value(v[i, h, m]) for i in BUS_SET, h in HOUR_SET, m in QUARTER_SET]
    global Pgen2 = reshape([ value(Pgen[i,h,m]) for i in (slack), h in HOUR_SET, m in QUARTER_SET],
                                                length(slack), length(HOUR_SET), length(QUARTER_SET))
    global Pdg2 = reshape([ value(Pdg[i,h,m]) for i in (DG_SET), h in HOUR_SET, m in QUARTER_SET],
                                                length(DG_SET), length(HOUR_SET), length(QUARTER_SET))
    global Ploss2 = reshape([ value(Ploss[(i,j),h,m]) for (i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
                                                length(BRANCH_SET), length(HOUR_SET), length(QUARTER_SET))
end
end # @time

if !converged
    println("WARNING: no convergence within $max_iter iterations (last MLE = $(dif_h[max_iter])).")
end

t_opt_c = zeros(max_iter)
for i in 1:conv_it
    t_opt_c[i] = sum(t_opt[1:i])
end
println("Optimization times (s): ", t_opt)
println("Cumulative optimization time (s): ", t_opt_c)
println("MLE per iteration (p.u.): ", dif_h)
println("MLE per iteration (x 1e-3 p.u.): ", round.(dif_h*1000, digits=7))

voltage_range = "Voltage Range: " * string(round(minimum(V_min), digits=4)) * " : " * string(round(maximum(V_max), digits=4)) * " (p.u.)"

Logging.disable_logging(Logging.Info)

Pdg_plot()
vp_plot()
