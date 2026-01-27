# Linear IVACOPF 1.0 - 24 Hours with 15 Minutes resolution - One Phase - 33 Bus
# Added breakpoints voltages
# Non Linear λ Method

using JuMP
using Gurobi
# using GLPK
# using Ipopt
# using AmplNLWriter, Bonmin_jll
using DataFrames, CSV, Plots, Logging

Logging.disable_logging(Logging.Warn)
# print("\033c")
#%%

function Pdg_plot()
    hour_labels = [lpad(h, 2, "0") for h in 0:23]
    hour_indices = 1:4:96  # [1, 5, 9, ..., 93]

    plt = plot(
        label="DG Active Power - IVACOPF", 
        xlabel="Hour", 
        ylabel="Active Power (p.u.)", 
        title="DG Active Power - IVACOPF",
        legend=:topleft,
        legend_background_color = :transparent,
        legend_foreground_color = RGBA(0, 0, 0, 0.3),
        xticks=(hour_indices, hour_labels),
        # xrotation=45,  # Rotate labels for readability
        # ylim=(0, 0.1)
        size=(800, 600)  # (width=800px, height=600px)
    )
    
    plP = reshape(permutedims(Pdg2, [1, 3, 2]), length(DG_SET), length(QUARTER_SET)*length(HOUR_SET))
    for i in 1:length(DG_SET)
        plP_max = reshape(Pdg_max_vary[DG_SET[i]]', length(QUARTER_SET)*length(HOUR_SET), 1)
        d = DG_SET[i]
        plot!(1:96, plP[i,:], label="DG $d", linewidth=2)
        # hline!([Pdg_max[d]], linestyle=:dash, label=false, linewidth=0.5)
        # plot!(1:96, plP_max, linestyle=:dash, linewidth=0.8, label="DG $d MAX")
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
        # linestyle=:dash,
        linewidth=1
        )

    plot!(1:96, sum(reshape(permutedims(Pload, [1, 3, 2]), length(BUS_SET),length(QUARTER_SET)*length(HOUR_SET)), dims=1)',
        label="Load",
        color=:green,
        # linestyle=:dash,
        linewidth=1
        )
    
    plot!(1:96, reshape(sum(permutedims(Ploss2, [1, 3, 2]), dims=1), length(QUARTER_SET)*length(HOUR_SET)),
        label="Ploss",
        color=:red,
        # linestyle=:dash,
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
    # savefig(plt, "D:\\Education\\Apply\\AmB\\Codes_LinDistFlow\\Pieces\\Saved_Results\\LD_P3.4.3_α$(round(α, digits=2))β$(round(β, digits=2)).svg")
end

function vp_plot()
    plt = plot(
        #label="Max Voltage", 
        xlabel="Bus", 
        ylabel="Voltage (p.u.)", 
        title="Voltage Profile - IVACOPF", 
        xticks=(BUS_SET, BUS_SET),
        ylim=(0.89, 1.05),
        yticks=(0.9:0.02:1.05),
        legend=:bottomright
        # ,size=(800, 600)  # (width=800px, height=600px)
    )
    plot!(BUS_SET, V_max, label="Max Voltage", linewidth=2, color="darkorange2")
    plot!(BUS_SET, V_max, fillrange=V_min, fillalpha=0.2, label=false, linecolor=:transparent, color="dodgerblue")
    plot!(BUS_SET, V_min, label="Min Voltage", linewidth=2, color="dodgerblue4")
    annotate!(1, 1.05, 
        # text(fairness_prnt * "\n" *DG_Miss_prnt * "\n" * voltage_limit * "\n" * voltage_range * "\n" * ΣΔVoltage * "\n" * ΣΔVoltage2,
        text(voltage_range,
        :left,
        :top,
        6,
        color=:red))

    hline!([1], linestyle=:dash, color=:gray40, label=false, linewidth=0.5)
    hline!([V_limit[1]], linestyle=:dash, color=:gray40, label=false, linewidth=0.5)
    hline!([V_limit[2]], linestyle=:dash, color=:gray40, label=false, linewidth=0.5)
    # for i in 1:6
    #     hline!([value.(V_break[i, DG_SET[1], 6])], linestyle=:dash, color=:gray40, label=false, linewidth=0.5)
    # end

    display(plt)
end

HOUR_SET = 1:24
QUARTER_SET = 1:4
# HQ_SET = 1:length(QUARTER_SET)*length(HOUR_SET)
HQ_SET = vcat([collect((h*4-4).+QUARTER_SET) for h in HOUR_SET]...)

BUS_SET = 1:33
slack = 1
Vnom = 1

Path = @__DIR__
System_Data = Matrix(DataFrame(CSV.File(Path * "/33_bus_syst_data.csv", skipto=1)))
solar_profile = reshape(DataFrame(CSV.File(Path * "/_solar_profile.csv"))[:,2],length(QUARTER_SET), length(HOUR_SET))'/100 # p.u.

# --- Base values
Sbase = parse(Float64, System_Data[2,9])*1e3 # in VA
Vbase = parse(Float64, System_Data[2,8])*1e3 # in V
Zbase = Vbase^2/Sbase # in Ω

Pload_std = reshape(parse.(Float64, String.(collect(skipmissing(System_Data[BUS_SET.+1, 4])))), :, 1) * 1e3 # in w
Qload_std = reshape(parse.(Float64, String.(collect(skipmissing(System_Data[BUS_SET.+1, 5])))), :, 1) * 1e3 # in VAr
Zdata = reshape(parse.(Float64, String.(collect(skipmissing(hcat(System_Data[BUS_SET[2:end].+1, 2:3], System_Data[BUS_SET[2:end].+1, 6:7]))))), :, 4) # in Ω


_15min_variation = DataFrame(CSV.File(Path * "/_15min_variation.csv", skipto=1))
load_category = Dict(1 => "Industrial", 2 => "Commercial", 3 => "Residential")
class = Dict(parse.(Int, _15min_variation[2:34,10]) .=> parse.(Int, _15min_variation[2:34,11]))
_15min_percent_matrix = Matrix{Float64}(tryparse.(Float64, _15min_variation[HQ_SET.+1, 2:7])) *1e-2
# _15min_percent_matrix = Matrix{Float64}(tryparse.(Float64, _15min_variation[2:97, 2:7])) *1e-2
# _15min_percent_matrix = Matrix{Float64}([tryparse(Float64, x) for x in values(_15min_variation[HQ_SET.+1, 2:7])]') * 1e-2

_15min_percent_p = _15min_percent_matrix[:,1:3]
_15min_percent_q = _15min_percent_matrix[:,4:6]

Pload = zeros(length(BUS_SET),24,4)
Qload = zeros(length(BUS_SET),24,4)

for bus in BUS_SET[2:end]
    # 33×24×4
    for h in 1:length(HOUR_SET)
        for q in 1:length(QUARTER_SET)
            Pload[bus, HOUR_SET[h], QUARTER_SET[q]] = Pload_std[bus] * _15min_percent_p[q+(h-1)*length(QUARTER_SET), class[bus]]/Sbase
            Qload[bus, HOUR_SET[h], QUARTER_SET[q]] = Qload_std[bus] * _15min_percent_q[q+(h-1)*length(QUARTER_SET), class[bus]]/Sbase
        end
    end

end


R = Dict()
X = Dict()

# add line RX
for k in 1:size(Zdata,1)
    R[Int(Zdata[k,1]),Int(Zdata[k,2])] = Zdata[k,3]/Zbase
    X[Int(Zdata[k,1]),Int(Zdata[k,2])] = Zdata[k,4]/Zbase
end

BRANCH_SET = [(Int(Zdata[k,1]), Int(Zdata[k,2])) for k in 1:size(Zdata,1)]
Bi_BRANCH_SET = vcat(BRANCH_SET, [(b,a) for (a,b) in BRANCH_SET])

# DG_SET = [3]
# DG_SIZE =  [0; # in w
#             0] # in VA

# DG_SET = [6 30]
# DG_SIZE =  [2e6 1.5e6; # in w
#             2.2e6 1.8e6] # in VA

# DG_SET = [6]
# DG_Psize = [2e6] # Max Pdg in w
# DG_Ssize = 1.1 .* DG_Psize # Max Sdg in VA
# DG_SIZE =  [DG_Psize;  # Max Pdg in w
#             DG_Ssize] # Max Sdg in VA
# DG_SET = [6 30]
# DG_Psize = [2e6 1.7e6] # Max Pdg in w
# DG_Ssize = 1.1 .* DG_Psize # Max Sdg in VA
# DG_SIZE =  [DG_Psize;  # Max Pdg in w
#             DG_Ssize] # Max Sdg in VA
# DG_SET = [4 7 17 21 31]
# DG_SIZE =  [0.5e6 0.4e6 1e6 0.3e6 1e6;  # Max Pdg in w
#             0.8e6 0.6e6 1.2e6 0.5e6 1.4e6] # Max Sdg in VA

# DG_SET = [7 17 31]
DG_SET = [7 18 33]
# DG_Psize = [1.8e6 1.3e6 0.9e6] # Max Pdg in w
DG_Psize = [2.2e6 1e6 1.5e6] # Max Pdg in w
DG_Ssize = 1.1 .* DG_Psize # Max Sdg in VA
DG_SIZE =  [DG_Psize;  # Max Pdg in w
            DG_Ssize] # Max Sdg in VA
# DG_SET = [4 7 17 21 31]
# DG_Psize = [0.7e6 0.6e6 0.8e6 0.5e6 1e6] # Max Pdg in w
# DG_Ssize = 1.1 .* DG_Psize # Max Sdg in VA
# DG_SIZE =  [DG_Psize;  # Max Pdg in w
#             DG_Ssize] # Max Sdg in VA

dg_number = length(DG_SET)
Pdg_max = Dict(DG_SET[i] => DG_SIZE[1,i]/Sbase for i in 1:dg_number) # p.u.
Sdg_max = Dict(DG_SET[i] => DG_SIZE[2,i]/Sbase for i in 1:dg_number) # p.u.
# Qdg_max = Dict(DG => sqrt(Sdg_max[DG]^2 - Pdg_max[DG]^2) for DG in DG_SET) # (p.u.)
Qdg_max = Sdg_max # (p.u.)
NON_DG_SET=setdiff(BUS_SET,DG_SET,slack)

Pdg_max_vary = Dict(DG => Pdg_max[DG] * solar_profile for DG in DG_SET)

Pcap = Dict(DG => [maximum(Pdg_max_vary[DG][h, :]) for h in HOUR_SET] for DG in DG_SET)
Qcap = Dict(DG => sqrt.(Sdg_max[DG]^2 .- (Pcap[DG].^2)) for DG in DG_SET)

qG_values = Dict(dg => 
    vcat(Qcap[dg]', Qcap[dg]', zeros(1, length(HOUR_SET)), zeros(1, length(HOUR_SET)), -Qcap[dg]', -Qcap[dg]' ) for dg in vec(DG_SET))

V_limit = [0.95 1.05] # in p.u.
global Vbp_limit = [0.85 1.15] # in p.u.

# global v_r_pr = ones(length(BUS_SET), length(HOUR_SET), length(QUARTER_SET))
global v_r_pr = ones(length(BUS_SET), 24, 4)
# global v_r_pr = v_r_prX
# global v_im_pr = zeros(length(BUS_SET), length(HOUR_SET), length(QUARTER_SET))
global v_im_pr = zeros(length(BUS_SET), 24, 4)
# global v_im_pr = v_im_prX
# global v_im_pr = VTM
global Ibs_r_pr = zeros(length(BUS_SET), 24, 4)
global Ibs_im_pr =zeros(length(BUS_SET), 24, 4)
# global Ibs_r_pr = zeros(length(BUS_SET), length(HOUR_SET), length(QUARTER_SET))
# global Ibs_im_pr =zeros(length(BUS_SET), length(HOUR_SET), length(QUARTER_SET))
global Ibr_r_pr = Dict((branch, hour, quarter) => 0.0 for branch in BRANCH_SET, hour in 1:24, quarter in 1:4)
global Ibr_im_pr = Dict((branch, hour, quarter) => 0.0 for branch in BRANCH_SET, hour in 1:24, quarter in 1:4)
# global Ibr_r_pr = Dict((branch, hour, quarter) => 0.0 for branch in BRANCH_SET, hour in HOUR_SET, quarter in QUARTER_SET)
# global Ibr_im_pr = Dict((branch, hour, quarter) => 0.0 for branch in BRANCH_SET, hour in HOUR_SET, quarter in QUARTER_SET)

@time begin
#%% ▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️▫️
max_iter = 15
conv_it = max_iter
dif_h = zeros(max_iter)
tol = 1e-3
# tol = 0.003
# tol = exp10(-3.7)
for it in 1:max_iter
    global v_r_pr
    global v_im_pr
    global Ibs_r_pr
    global Ibs_im_pr
    global Ibr_r_pr
    global Ibr_im_pr
    global Vbp_limit

    println("🔻🔺🔻🔺🔻🔺🔻🔺🔻🔺🔻🔺🔻 Iteration $it 🔻🔺🔻🔺🔻🔺🔻🔺🔻🔺🔻🔺🔻")

    global OPF = Model(Gurobi.Optimizer)
    # set_optimizer_attribute(OPF, "OutputFlag", 0)
    # global OPF = Model(GLPK.Optimizer)
    # global OPF = Model(Ipopt.Optimizer)
    Time_Limit = 3600 # in seconds
    set_optimizer_attribute(OPF, "TimeLimit", Time_Limit) # Set a time limit
    # General knobs
    set_optimizer_attribute(OPF, "MIPGap", 0.001)  # Balance speed and accuracy (?% gap)
    # set_optimizer_attribute(OPF, "Cuts", 0)  # Use aggressive cuts to reduce problem size
    # set_optimizer_attribute(OPF, "Method", 2)  # Use barrier solver for faster performance
    # set_optimizer_attribute(OPF, "Crossover", 0)  # Disable crossover for barrier method
    # set_optimizer_attribute(OPF, "OutputFlag", 1)  # Print solver progress
    # set_optimizer_attribute(OPF, "DisplayInterval", 5)  # Show logs every 5 seconds
    set_optimizer_attribute(OPF, "Threads", 0)  # all cores
    # set_optimizer_attribute(OPF, "MIPFocus", 1)  # Focus on getting feasible solutions first
    # set_optimizer_attribute(OPF, "Heuristics", 0.2)  # More heuristics for faster feasibility
    # set_optimizer_attribute(OPF, "Presolve", 2)  # Strong presolve to simplify the mode

    #Initializing parameters
    @variable(OPF, 0 <= Pgen[slack, HOUR_SET, QUARTER_SET])  # Active power generation at slack bus
    @variable(OPF, Qgen[slack, HOUR_SET, QUARTER_SET])  # Reactive power generation at slack bus
    @variable(OPF, V_limit[1] <= v[BUS_SET, HOUR_SET, QUARTER_SET] <= V_limit[2], start = Vnom) # voltage at bus i
    @variable(OPF, Δv[BUS_SET, HOUR_SET, QUARTER_SET])
    @variable(OPF, v_r[BUS_SET, HOUR_SET, QUARTER_SET], start = Vnom) # Real part of the complex voltage at bus i
    @variable(OPF, v_im[BUS_SET, HOUR_SET, QUARTER_SET], start = 0) # Imaginary part of the complex voltage at bus i
    # @variable(OPF, 0 <= v_variation[BUS_SET, HOUR_SET, QUARTER_SET])
    @variable(OPF, Ibr_r[BRANCH_SET, HOUR_SET, QUARTER_SET])
    @variable(OPF, Ibr_im[BRANCH_SET, HOUR_SET, QUARTER_SET])
    @variable(OPF, Psnd[Bi_BRANCH_SET, HOUR_SET, QUARTER_SET])
    @variable(OPF, Qsnd[Bi_BRANCH_SET, HOUR_SET, QUARTER_SET])
    @variable(OPF, 0 <= Ploss[BRANCH_SET, HOUR_SET, QUARTER_SET])
    @variable(OPF, 0 <= Qloss[BRANCH_SET, HOUR_SET, QUARTER_SET])
    @variable(OPF, Ibs_r[BUS_SET, HOUR_SET, QUARTER_SET]) 
    @variable(OPF, Ibs_im[BUS_SET, HOUR_SET, QUARTER_SET])
    @variable(OPF, 0 <= Pdg[DG_SET, HOUR_SET, QUARTER_SET])  # DG Active power generation
    @variable(OPF, Qdg[DG_SET, HOUR_SET, QUARTER_SET])  # DG Reactive power generation
    @variable(OPF, 0 <= fairness_index[HOUR_SET, QUARTER_SET])

    # λ method for v
    #----------------------------------------------------------------------------------------------------------
    @variable(OPF, λ[1:6, DG_SET, HOUR_SET, QUARTER_SET] >= 0)
    @variable(OPF, δ[1:5, DG_SET, HOUR_SET, QUARTER_SET], Bin)
    @variable(OPF, V_break[1:6, DG_SET, HOUR_SET])

    @constraint(OPF, [d in DG_SET, h in HOUR_SET], V_break[1, d, h] == Vbp_limit[1])
    @constraint(OPF, [d in DG_SET, h in HOUR_SET], V_break[6, d, h] == Vbp_limit[2])
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

    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], v[d,h,m] == sum(λ[i, d, h, m] * V_break[i, d, h] for i in 1:6))

    # @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in MONTH_SET], λ[:, d, h, m] in MOI.SOS2([1, 2, 3, 4, 5, 6]))
    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], sum(λ[i, d, h, m] for i in 1:6) == 1)

    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], λ[1, d, h, m] <= δ[1, d, h, m])
    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], λ[2, d, h, m] <= δ[1, d, h, m] + δ[2, d, h, m])
    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], λ[3, d, h, m] <= δ[2, d, h, m] + δ[3, d, h, m])
    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], λ[4, d, h, m] <= δ[3, d, h, m] + δ[4, d, h, m])
    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], λ[5, d, h, m] <= δ[4, d, h, m] + δ[5, d, h, m])
    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], λ[6, d, h, m] <= δ[5, d, h, m])
    @constraint(OPF, [d in DG_SET, h in HOUR_SET, m in QUARTER_SET], sum(δ[:, d, h, m]) == 1)  # Only one pair of adjacent segments can be active

    @constraint(OPF, Reactive_DG[d in DG_SET, h in HOUR_SET, m in QUARTER_SET],
            Qdg[d, h, m] == sum(λ[i, d, h, m] * qG_values[d][i, h] for i in 1:6))


    #----------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------
    @constraint(OPF, DG_max1[DG in DG_SET, h in HOUR_SET, m in QUARTER_SET], Pdg[DG, h, m] <= Pdg_max_vary[DG][h,m])
    @constraint(OPF, [i in DG_SET, h in HOUR_SET, m in QUARTER_SET], Pdg[i,h,m] <= Pdg_max[i])
    @constraint(OPF, [i in DG_SET, h in HOUR_SET, m in QUARTER_SET], Qdg[i,h,m] <= Qdg_max[i])

    #Defining constraints
    @constraint(OPF, [i in slack, h in HOUR_SET, m in QUARTER_SET], v_r[i,h,m] == Vnom)
    @constraint(OPF, [i in slack, h in HOUR_SET, m in QUARTER_SET], v_im[i,h,m] == 0)

    # Line Current Flow Constraints
    @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], v_r[i,h,m]-v_r[j,h,m] == R[i, j]*Ibr_r[(i, j),h,m] - X[i, j]*Ibr_im[(i, j),h,m] )
    @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], v_im[i,h,m]-v_im[j,h,m] == R[i, j]*Ibr_im[(i, j),h,m] + X[i, j]*Ibr_r[(i, j),h,m] )

    # Bus Current Injection Constraints
    @constraint(OPF, [bus in BUS_SET, h in HOUR_SET, m in QUARTER_SET], Ibs_r[bus,h,m] ==  sum(Ibr_r[(bus, j),h,m] for (i, j) in BRANCH_SET if i == bus) 
                                                    - sum(Ibr_r[(i, bus),h,m] for (i, j) in BRANCH_SET if j == bus))

    @constraint(OPF, [bus in BUS_SET, h in HOUR_SET, m in QUARTER_SET], Ibs_im[bus,h,m] ==  sum(Ibr_im[(bus, j),h,m] for (i, j) in BRANCH_SET if i == bus) 
                                                    - sum(Ibr_im[(i, bus),h,m] for (i, j) in BRANCH_SET if j == bus))

    # Power Balance Constraints
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

    # @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Ploss[(i, j),h,m] == (v_r[i,h,m]-v_r[j,h,m])*Ibr_r[(i, j),h,m] + (v_im[i,h,m]-v_im[j,h,m])*Ibr_im[(i, j),h,m])
    @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Ploss[(i, j), h, m] ==
                                                Ibr_r_pr[(i, j), h, m] * (v_r[i, h, m] - v_r[j, h, m]) +
                                                (v_r_pr[i, h, m] - v_r_pr[j, h, m]) * Ibr_r[(i, j), h, m] +
                                                Ibr_im_pr[(i, j), h, m] * (v_im[i, h, m] - v_im[j, h, m]) +
                                                (v_im_pr[i, h, m] - v_im_pr[j, h, m]) * Ibr_im[(i, j), h, m] -
                                                ((v_r_pr[i, h, m] - v_r_pr[j, h, m]) * Ibr_r_pr[(i, j), h, m] +
                                                (v_im_pr[i, h, m] - v_im_pr[j, h, m]) * Ibr_im_pr[(i, j), h, m]))

    # # @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Ploss[(i, j),h,m] == ((v_r[i,h,m]-v_r[j,h,m])^2+(v_im[i,h,m]-v_im[j,h,m])^2) * R[i, j]/(R[i, j]^2+X[i, j]^2) )
    # @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Ploss[(i, j), h, m] ==
    #                                                             (R[i,j] / (R[i,j]^2 + X[i,j]^2)) * (
    #                                                                 2 * ( (v_r_pr[i,h,m] - v_r_pr[j,h,m])  * (v_r[i,h,m] - v_r[j,h,m])
    #                                                                     + (v_im_pr[i,h,m] - v_im_pr[j,h,m]) * (v_im[i,h,m] - v_im[j,h,m]) )
    #                                                                 - (v_r_pr[i,h,m]  - v_r_pr[j,h,m])^2
    #                                                                 - (v_im_pr[i,h,m] - v_im_pr[j,h,m])^2))

    # # @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Ploss[(i, j),h,m] == (Ibr_r[(i, j),h,m]^2 + Ibr_im[(i, j),h,m]^2) * R[i, j])
    # @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Ploss[(i, j), h, m] == R[i, j] * (
    #                                                                     2 * Ibr_r_pr[(i, j), h, m]  * Ibr_r[(i, j), h, m] - Ibr_r_pr[(i, j), h, m]^2 +
    #                                                                     2 * Ibr_im_pr[(i, j), h, m] * Ibr_im[(i, j), h, m] - Ibr_im_pr[(i, j), h, m]^2))

    # @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Ploss[(i, j), h, m] ==
    #                                                                         R[i, j] * (2 * Ibr_r_pr[(i, j), h, m] * Ibr_r[(i, j), h, m] +
    #                                                                         2 * Ibr_im_pr[(i, j), h, m] * Ibr_im[(i, j), h, m] -
    #                                                                         (Ibr_r_pr[(i, j), h, m]^2 + Ibr_im_pr[(i, j), h, m]^2)))


    # @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Qloss[(i, j),h,m] == (v_im[i,h,m]-v_im[j,h,m])*Ibr_r[(i, j),h,m] - (v_r[i,h,m]-v_r[j,h,m])*Ibr_im[(i, j),h,m])
    @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Qloss[(i, j), h, m] ==
                                        # + (v_im[i] * Ibr_r) term
                                        v_im_pr[i, h, m] * Ibr_r[(i, j), h, m] + Ibr_r_pr[(i, j), h, m] * v_im[i, h, m] - v_im_pr[i, h, m] * Ibr_r_pr[(i, j), h, m]
                                        # - (v_im[j] * Ibr_r) term
                                        - v_im_pr[j, h, m] * Ibr_r[(i, j), h, m] - Ibr_r_pr[(i, j), h, m] * v_im[j, h, m] + v_im_pr[j, h, m] * Ibr_r_pr[(i, j), h, m]
                                        # - (v_r[i] * Ibr_im) term
                                        - v_r_pr[i, h, m] * Ibr_im[(i, j), h, m] - Ibr_im_pr[(i, j), h, m] * v_r[i, h, m] + v_r_pr[i, h, m] * Ibr_im_pr[(i, j), h, m]
                                        # + (v_r[j] * Ibr_im) term
                                        + v_r_pr[j, h, m] * Ibr_im[(i, j), h, m] + Ibr_im_pr[(i, j), h, m] * v_r[j, h, m] - v_r_pr[j, h, m] * Ibr_im_pr[(i, j), h, m])
    # @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Qloss[(i, j)] == (v_im_pr[i,h,m] - v_im_pr[j]) * Ibr_r_pr[(i, j),h,m]
    #                                                         - (v_r_pr[i,h,m] - v_r_pr[j]) * Ibr_im_pr[(i, j),h,m]
    #                                                         + Ibr_r_pr[(i, j),h,m] * (v_im[i] - v_im_pr[i,h,m])
    #                                                         - Ibr_r_pr[(i, j),h,m] * (v_im[j] - v_im_pr[j])
    #                                                         - Ibr_im_pr[(i, j),h,m] * (v_r[i] - v_r_pr[i,h,m])
    #                                                         + Ibr_im_pr[(i, j),h,m] * (v_r[j] - v_r_pr[j])
    #                                                         + (v_im_pr[i,h,m] - v_im_pr[j]) * (Ibr_r[(i, j)] - Ibr_r_pr[(i, j),h,m])
    #                                                         - (v_r_pr[i,h,m] - v_r_pr[j]) * (Ibr_im[(i, j)] - Ibr_im_pr[(i, j),h,m]))

    # @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Qloss[(i, j),h,m] == ((v_r[i,h,m]-v_r[j,h,m])^2+(v_im[i,h,m]-v_im[j,h,m])^2) * X[i, j]/(R[i, j]^2+X[i, j]^2))
    # @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Qloss[(i, j), h, m] == X[i, j] / (R[i, j]^2 + X[i, j]^2) * (
    #                                                             2 * ( (v_r_pr[i,h,m]  - v_r_pr[j,h,m])  * (v_r[i,h,m]  - v_r[j,h,m])
    #                                                             + (v_im_pr[i,h,m] - v_im_pr[j,h,m]) * (v_im[i,h,m] - v_im[j,h,m]) )
    #                                                             - (v_r_pr[i,h,m]  - v_r_pr[j,h,m])^2 - (v_im_pr[i,h,m] - v_im_pr[j,h,m])^2 ))

    
    # @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Qloss[(i, j),h,m] == (Ibr_r[(i, j),h,m]^2 + Ibr_im[(i, j),h,m]^2) * X[i, j])
    # @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Qloss[(i, j), h, m] == X[i, j] * (
    #                                                             2 * Ibr_r_pr[(i, j), h, m]  * Ibr_r[(i, j), h, m] - Ibr_r_pr[(i, j), h, m]^2
    #                                                             + 2 * Ibr_im_pr[(i, j), h, m] * Ibr_im[(i, j), h, m] - Ibr_im_pr[(i, j), h, m]^2))


    @constraint(OPF, [bus in slack, h in HOUR_SET, m in QUARTER_SET], Pgen[bus,h,m] == sum(Psnd[(i, j),h,m] for (i, j) in Bi_BRANCH_SET if i == bus) + Pload[bus,h,m] )
    @constraint(OPF, [bus in DG_SET, h in HOUR_SET, m in QUARTER_SET], Pdg[bus,h,m] == sum(Psnd[(i, j),h,m] for (i, j) in Bi_BRANCH_SET if i == bus) + Pload[bus,h,m] )
    @constraint(OPF, [bus in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET], 0 == sum(Psnd[(i, j),h,m] for (i, j) in Bi_BRANCH_SET if i == bus) + Pload[bus,h,m] )

    @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Psnd[(i, j),h,m] == Ploss[(i, j),h,m] - Psnd[(j, i),h,m])

    @constraint(OPF, [bus in slack, h in HOUR_SET, m in QUARTER_SET], Qgen[bus,h,m] == sum(Qsnd[(i, j),h,m] for (i, j) in Bi_BRANCH_SET if i == bus) + Qload[bus,h,m] )
    @constraint(OPF, [bus in DG_SET, h in HOUR_SET, m in QUARTER_SET], Qdg[bus,h,m] == sum(Qsnd[(i, j),h,m] for (i, j) in Bi_BRANCH_SET if i == bus) + Qload[bus,h,m] )
    @constraint(OPF, [bus in NON_DG_SET, h in HOUR_SET, m in QUARTER_SET], 0 == sum(Qsnd[(i, j),h,m] for (i, j) in Bi_BRANCH_SET if i == bus) + Qload[bus,h,m] )

    @constraint(OPF, [(i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], Qsnd[(i, j),h,m] == Qloss[(i, j),h,m] - Qsnd[(j, i),h,m])

    # Linearized Voltage Magnitude Constraints
    # @constraint(OPF, [i in BUS_SET], v[i] == sqrt(v_r[i]^2 + v_im[i]^2))
    @constraint(OPF, [i in BUS_SET, h in HOUR_SET, m in QUARTER_SET], v[i,h,m] == ( (v_r_pr[i,h,m] / sqrt(v_r_pr[i,h,m]^2 + v_im_pr[i,h,m]^2)) * v_r[i,h,m] )
                        + ( (v_im_pr[i,h,m] / sqrt(v_r_pr[i,h,m]^2 + v_im_pr[i,h,m]^2)) * v_im[i,h,m] ))

    # ρ = 0.01  # bound
    # @constraint(OPF, [(i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], -ρ <= (v_r[i,h,m]-v_r[j,h,m])  - (v_r_pr[i,h,m]-v_r_pr[j,h,m])  <= ρ)
    # @constraint(OPF, [(i,j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET], -ρ <= (v_im[i,h,m]-v_im[j,h,m]) - (v_im_pr[i,h,m]-v_im_pr[j,h,m]) <= ρ)

    @constraint(OPF, [i in BUS_SET, h in HOUR_SET, m in QUARTER_SET], v[i,h,m]-Vnom <= Δv[i,h,m]) 
    @constraint(OPF, [i in BUS_SET, h in HOUR_SET, m in QUARTER_SET], -v[i,h,m]+Vnom <= Δv[i,h,m]) 

    @constraint(OPF, Fairness[DG in DG_SET, h in HOUR_SET, m in QUARTER_SET], fairness_index[h,m] * Pdg_max_vary[DG][h,m] >=
                                                                                            (Pdg_max_vary[DG][h,m]-Pdg[DG, h, m]) )
    # @constraint(OPF, Fairness[DG in DG_SET, h in HOUR_SET, m in QUARTER_SET], fairness_index[h,m] >=
    #                                                                                         (Pdg_max_vary[DG][h,m]-Pdg[DG, h, m]) )
    
    # @objective(OPF, Min, sum(Ploss[(i, j), h, m] for (i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET))
    # @objective(OPF, Min, sum(Pdg_max[i]-Pdg[i] for i in DG_SET))
    # @objective(OPF, Min, sum(Δv[i,h,m] for i in BUS_SET, h in HOUR_SET, m in QUARTER_SET))
    # @objective(OPF, Min, sum(fairness_index[h,m] for h in HOUR_SET, m in QUARTER_SET))
    # @objective(OPF, Min, sum(fairness_index[h,m] for h in HOUR_SET, m in QUARTER_SET)
    #                     - 0.01*sum(v[i,h,m] for i in BUS_SET, h in HOUR_SET, m in QUARTER_SET))

    @variable(OPF, v_varx)
    @constraint(OPF, [i in BUS_SET, h in HOUR_SET, m in QUARTER_SET],  v_varx >= v[i,h,m] - Vnom)
    @constraint(OPF, [i in BUS_SET, h in HOUR_SET, m in QUARTER_SET],  v_varx >= -v[i,h,m] + Vnom)
    # @objective(OPF, Min, v_varx)
    # @objective(OPF, Max, sum(v[i,h,m] for i in BUS_SET, h in HOUR_SET, m in QUARTER_SET))

    α = 1
    # β = 1
    # Γ = 1-α-β
    β = 1-α
    # @constraint(OPF, sum(fairness_index[h,m] for h in HOUR_SET, m in QUARTER_SET) <= 5)

    # f1 = 0.5407519    # DG_SET = [6 30]# DG_Psize = [2e6 1.7e6] # Max Pdg in w
    # f2 = 47.4487479    # DG_SET = [6 30]# DG_Psize = [2e6 1.7e6] # Max Pdg in w
    # # f3 = 1
    # f1 = 0.4    # DG_SET = [4 7 17 21 31] # DG_Psize = [0.7e6 0.6e6 0.8e6 0.5e6 1e6] # Max Pdg in w
    # f2 = 34.7456195    # DG_SET = [4 7 17 21 31] # DG_Psize = [0.7e6 0.6e6 0.8e6 0.5e6 1e6] # Max Pdg in w
    # # f3 = 1
    # f1 = 0.0231253    # DG_SET = [7 17 31] # DG_Psize = [1.4e6 1.1e6 1e6] # Max Pdg in w
    # f2 = 29.7876619    # DG_SET = [7 17 31] # DG_Psize = [1.4e6 1.1e6 1e6] # Max Pdg in w
    # f3 = 1

    f1 = 0.4474876    # Final PVs - Variable V_break + Egal.
    f2 = 26.5214655    # Final PVs - Variable V_break + prop

    # @objective(OPF, Min, α*sum(fairness_index[h,m] for h in HOUR_SET, m in QUARTER_SET)/f1 + β*v_varx/f2)
    @objective(OPF, Min, α*sum(fairness_index[h,m] for h in HOUR_SET, m in QUARTER_SET)/f1 +
                                β*sum(Δv[i,h,m] for i in BUS_SET, h in HOUR_SET, m in QUARTER_SET)/f2)
    # @objective(OPF, Min, α*sum(fairness_index[h,m] for h in HOUR_SET, m in QUARTER_SET)/f1 + β*v_varx/f2 #)
    #         + Γ*sum(Ploss[(i, j), h, m] for (i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET)/f3)

    println("----------------------------------------------------------------------------")
    println("Number of Variables: ", num_variables(OPF))
    println("Total Constraints: ", sum(length(all_constraints(OPF, F, S)) for (F, S) in list_of_constraint_types(OPF)))
    println("Number of Binary Variables: ", count(is_binary, all_variables(OPF)))
    println("Number of Integer Variables: ", count(is_integer, all_variables(OPF)))
    # println("Objective Sense: ", objective_sense(OPF))
    println("----------------------------------------------------------------------------\n")

    # optimize!(OPF)
    @time status = optimize!(OPF)

    # max_diff = maximum([abs(value(Ibs_r[n]) - Ibs_r_pr[n]) for n in BUS_SET])
    # max_diff = maximum([abs(value(Ibr_r[(i, j)]) - Ibr_r_pr[(i, j),h,m]) for (i, j) in BRANCH_SET])
    # max_diff = maximum([abs(value(v_r[n]) - v_r_pr[n]) for n in BUS_SET])
    # max_diff = maximum([abs(value(v_im[n]) - v_im_pr[n]) for n in BUS_SET])
    max_diff = maximum([abs(value((v_r[i,h,m]-v_r[j,h,m])*Ibr_r[(i, j),h,m] + (v_im[i,h,m]-v_im[j,h,m])*Ibr_im[(i, j),h,m]
                - Ploss[(i, j),h,m])) for (i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET])
    # max_diff = maximum([abs(value((Ibr_r[(i, j),h,m]^2 + Ibr_im[(i, j),h,m]^2) * R[i, j] - Ploss[(i, j),h,m]))
    #                 for (i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET])

    # max_diff = maximum([abs(value(Qloss[(i, j)]) - Qloss_pr[(i, j)]) for (i, j) in BRANCH_SET])
    # Qloss_pr = Dict{Tuple{Int64, Int64}, Float64}(branch => value for (branch, value) in zip(BRANCH_SET, value.(Qloss)))
    dif_h[it] = max_diff
    # println("🔥 max_diff: ", max_diff )

    # v_r_pr = value.(v_r).data
    [v_r_pr[b, h, q] = value(v_r[b, h, q]) for b in BUS_SET, h in HOUR_SET, q in QUARTER_SET]
    # # v_r_pr = ones(length(BUS_SET), length(HOUR_SET), length(QUARTER_SET))
    # # v_im_pr = value.(v_im).data
    [v_im_pr[b, h, q] = value(v_im[b, h, q]) for b in BUS_SET, h in HOUR_SET, q in QUARTER_SET]
    # # v_im_pr = zeros(length(BUS_SET), length(HOUR_SET), length(QUARTER_SET))
    # # Ibs_r_pr = value.(Ibs_r).data
    [Ibs_r_pr[b, h, q] = value(Ibs_r[b, h, q]) for b in BUS_SET, h in HOUR_SET, q in QUARTER_SET]
    # # Ibs_im_pr = value.(Ibs_im).data
    [Ibs_im_pr[b, h, q] = value(Ibs_im[b, h, q]) for b in BUS_SET, h in HOUR_SET, q in QUARTER_SET]
    [Ibr_r_pr[(br, h, q)] = value(Ibr_r[br, h, q]) for br in BRANCH_SET, h in HOUR_SET, q in QUARTER_SET]
    [Ibr_im_pr[(br, h, q)] = value(Ibr_im[br, h, q]) for br in BRANCH_SET, h in HOUR_SET, q in QUARTER_SET]

    # println(v_r_pr)
    # println("Vr: ", round.(v_r_pr, digits=6), 
    #         "\nVim: ", round.(v_im_pr, digits=6),
    #         "\n|V|: ", round.([value(v[i]) for i in BUS_SET], digits=6),
    #         "\nIbs_r: ", round.(Ibs_r_pr, digits=6),
    #         "\nIbs_im: ", round.(Ibs_im_pr, digits=6),
    #         )

    # println("Ibr_r: ", round.([Ibr_r_pr[(i,j)] for (i, j) in BRANCH_SET], digits=6) )
    # println("Ibr_im: ", round.([Ibr_im_pr[(i,j)] for (i, j) in BRANCH_SET], digits=6) )

    # println("Ibs_r: ", round.(value.(Ibs_r[i] for i in BUS_SET), digits=6) )
    # println("Ibs_im: ", round.(value.(Ibs_im[i] for i in BUS_SET), digits=6) )

    # println("🔹🔸🔷🔶 Iteration $it 🔶🔷🔸🔹")
    println("🔹🔸🔷🔶 Iteration $it 🟨 Pslack: ", round.(sum(Array(value.(Pgen)), dims=(1,3))[:,:,1]*Sbase*1e-3, digits=3), " (kWh)")
    println("🔹🔸🔷🔶 Iteration $it 🟨 Qslack: ", round.(sum(Array(value.(Qgen)), dims=(1,3))[:,:,1]*Sbase*1e-3, digits=3), " (kVArh)")
    println("🔹🔸🔷🔶 Iteration $it 🟩 Pdg: ", round.(sum(Array(value.(Pdg)), dims=(1,3))[:,:,1]*Sbase*1e-3, digits=3), " (kWh)")
    # println("🔹🔸🔷🔶 Iteration $it 🟩 Pdg: ", round.(value.(Pdg[i,h,m] for i in DG_SET, h in HOUR_SET, m in QUARTER_SET)*Sbase*1e-3, digits=6), " (kW)")
    println("🔹🔸🔷🔶 Iteration $it 🟩 Qdg: ", round.(sum(Array(value.(Qdg)), dims=(1,3))[:,:,1]*Sbase*1e-3, digits=3), " (kVArh)")
    # println("🔹🔸🔷🔶 Iteration $it 🟩 Qdg: ", round.(value.(Qdg[i,h,m] for i in DG_SET, h in HOUR_SET, m in QUARTER_SET)*Sbase*1e-3, digits=6), " (kVAr)")
    println("🔹🔸🔷🔶 Iteration $it 🟥 ΣPloss: ", round.(sum(value.(Ploss[(i,j),h,m] for (i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET), dims=(1,3))[:,:,1]*Sbase*1e-3, digits=6), " (kW)")
    println("🔹🔸🔷🔶 Iteration $it 🟥 ΣQloss: ", round.(sum(value.(Qloss[(i,j),h,m] for (i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET), dims=(1,3))[:,:,1]*Sbase*1e-3, digits=6), " (kVAr)")
    println("🔹🔸🔷🔶 Iteration $it 🟦 ΣPload: ", round.(sum(Pload, dims=(1,3))[:,:,1]*Sbase*1e-3, digits=3), " (kWh)")
    println("🔹🔸🔷🔶 Iteration $it 🟦 ΣQload: ", round.(sum(Qload, dims=(1,3))[:,:,1]*Sbase*1e-3, digits=3), " (kVArh)")

    global MaxΔV = value(v_varx)
    fairness_prnt = "Fairness index: " *string(round(value(sum(fairness_index[h,m] for h in HOUR_SET, m in QUARTER_SET)), digits=7))
    println(fairness_prnt)

    DG_miss = value(sum(Pdg_max_vary[DG][h,m]-Pdg[DG, h, m] for DG in DG_SET, h in HOUR_SET, m in QUARTER_SET))
    DG_miss_prcnt_av = 100 * value(sum((Pdg_max_vary[DG][h,m]-Pdg[DG, h, m])/Pdg_max_vary[DG][h,m] 
                        for DG in DG_SET, h in HOUR_SET, m in QUARTER_SET if Pdg_max_vary[DG][h,m] != 0))/(dg_number*length(HQ_SET))
    DG_miss_prnt = "DG_Miss/|DG|.|T|: " * string(round(DG_miss/(dg_number*length(HQ_SET)), digits=7)) * " (p.u.) ➡️   " *
                string(round(DG_miss_prcnt_av, digits=7)) * " (%)"
    println(DG_miss_prnt)

    MaxΔV_prnt = "Max ΔV: " *string(round(MaxΔV, digits=7))
    println(MaxΔV_prnt)
    V = [value(v[i,h,m]) for i in sort(BUS_SET), h in HOUR_SET, m in QUARTER_SET]
    ΔVoltage2 = value(sum(Δv[i,h,m] for i in BUS_SET, h in HOUR_SET, m in QUARTER_SET))/(length(BUS_SET)*length(HQ_SET))
    ΣΔVoltage2 = "Σ|Vi - Vo|/|B|.|T|: " * string(round(ΔVoltage2, digits=7)*1e3) * " ×e-3 (p.u.)"
    println(ΣΔVoltage2)
    ΔVoltage1 = value(sum(Δv[i,h,m] for i in BUS_SET, h in HOUR_SET, m in QUARTER_SET))
    ΣΔVoltage1 = "Σ|Vi - Vo|*: " * string(round(ΔVoltage1, digits=7)) * " (p.u.)"
    println(ΣΔVoltage1)
    ΔVoltage = sum(abs.(V .- Vnom))
    ΣΔVoltage = "Σ|Vi - Vo|: " * string(round(ΔVoltage, digits=7)) * " (p.u.)"
    println(ΣΔVoltage)
    Ploss_prnt = "Total Ploss: " *string(round(value(sum(Ploss[(i, j), h, m] for (i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET)), digits=5))
    println(Ploss_prnt)
    global V_max = maximum(V, dims=[2,3])[:, 1, 1]  # Max over hours and months
    global V_min = minimum(V, dims=[2,3])[:, 1, 1]  # Min over hours and months

    println("🔥 max_diff: ", max_diff )
    if max_diff < tol
        global conv_it = it
        # global Vbp = value.(V_break[:,:,1])
        global V_ac_l2 = [value(v[i, h, m]) for i in BUS_SET, h in HOUR_SET, m in QUARTER_SET]
        global Pgen2 = reshape([ value(Pgen[i,h,m]) for i in (slack), h in HOUR_SET, m in QUARTER_SET],
                                                    length(slack), length(HOUR_SET), length(QUARTER_SET))
        global Pdg2 = reshape([ value(Pdg[i,h,m]) for i in (DG_SET), h in HOUR_SET, m in QUARTER_SET],
                                                    length(DG_SET), length(HOUR_SET), length(QUARTER_SET))
        global Ploss2 = reshape([ value(Ploss[(i,j),h,m]) for (i, j) in BRANCH_SET, h in HOUR_SET, m in QUARTER_SET],
                                                    length(BRANCH_SET), length(HOUR_SET), length(QUARTER_SET))
        # global Qdg = value.(Qdg)
        println("🚨🚨🚨🚨🚨🚨🚨🚨🚨 Converged on Iteration $it ➡  max_diff: $max_diff 🚨🚨🚨🚨🚨🚨🚨🚨🚨")
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
end #time ending
# V_min = minimum(V_ac_l2)
# V_max = maximum(V_ac_l2)

voltage_range = "Voltage Range: " * string(round(minimum(V_min), digits=4)) * " : " * string(round(maximum(V_max), digits=4)) * " (p.u.)"

Logging.disable_logging(Logging.Info)
# global_logger(ConsoleLogger(stderr, Logging.Info))

Pdg_plot()
vp_plot()

# plot(1:conv_it,dif_h[1:conv_it], title="Max. Diff.", xlabel="Iteration", ylabel="Max. Diff.", legend=false)
# plot(1:conv_it,log10.(dif_h[1:conv_it]), title="Max. Diff. Logarithm", xlabel="Iteration", ylabel="Max. Diff. logarithm", legend=false)
# plot(BUS_SET, (V_ac_l2), label="V IVACOPF Lin", xlabel="Bus", ylabel="Voltage (p.u.)", title="Voltage Profile - IVAC L 2.0", linewidth=1, color="darkorange1")
