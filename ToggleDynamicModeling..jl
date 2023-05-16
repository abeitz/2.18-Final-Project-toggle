__precompile__()

module ToggleDynamicModeling

## Load in packages
using DifferentialEquations
using DiffEqJump
using Plots

## Define structure and function exports
export generate_jump_problem
export Toggle_Parameters, DEFAULT_TOGGLE_PARAMS
export Cell_Parameters, DEFAULT_CELL_PARAMS
export Simulation_Parameters, DEFAULT_SIM_PARAMS

## Define parameter structures
mutable struct Toggle_Parameters
    α0::Float64 #Basal transcription rate constant [1/hr]
    α::Float64 #Activated transcription rate constant [1/hr]
    a::Float64 #Activator-DNA association rate constant [1/num*hr]
    d::Float64 #Activator-DNA dissociation rate constant [1/hr]
    a′::Float64 #miRNA-mRNA association rate constant [1/num*hr]
    d′::Float64 #miRNA-mRNA dissociation rate constant [1/hr]
    ks::Float64 #Splicing rate constant [1/hr]
    km::Float64 #miRNA-mediated degradation rate constant [1/hr]
    κ::Float64 #Translation rate constant [1/hr]
    δ::Float64 #RNA degradation rate constant [1/hr]
    γ::Float64 #Protein degradation rate constant [1/hr]
end

DEFAULT_TOGGLE_PARAMS = Toggle_Parameters(1,2,1,1,10,1,5,10,5,1,1)

mutable struct Cell_Parameters
    D1tot::Float64 #Node 1 DNA copy number
    D2tot::Float64 #Node 2 DNA copy number
    D3tot::Float64 #Downstream target DNA copy number
    D4tot::Float64 #Downstream target DNA copy number
    Ω::Float64 #System volume
end

DEFAULT_CELL_PARAMS = Cell_Parameters(10,10,0,0,1)

mutable struct Simulation_Parameters
    Toggle_Params::Toggle_Parameters
    Cell_Params::Cell_Parameters
end

DEFAULT_SIM_PARAMS = Simulation_Parameters(DEFAULT_TOGGLE_PARAMS,DEFAULT_CELL_PARAMS)

## Define stochastic model equations 
# Propensities
function rate_1(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.α0*u[1]
end

function rate_2(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.α0*u[2]
end

function rate_3(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.a*u[9]*u[1]/p.Cell_Params.Ω
end

function rate_4(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.d*u[11]
end

function rate_5(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.α*u[11]
end

function rate_6(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.a*u[10]*u[2]/p.Cell_Params.Ω
end

function rate_7(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.d*u[12]
end

function rate_8(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.α*u[12]
end

function rate_9(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.ks*u[3]
end

function rate_10(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.ks*u[4]
end

function rate_11(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.a′*u[5]*u[8]/p.Cell_Params.Ω
end

function rate_12(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.d′*u[13]
end

function rate_13(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.km*u[13]
end

function rate_14(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.a′*u[6]*u[7]/p.Cell_Params.Ω
end

function rate_15(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.d′*u[14]
end

function rate_16(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.km*u[14]
end

function rate_17(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.a′*u[3]*u[8]/p.Cell_Params.Ω
end

function rate_18(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.d′*u[15]
end

function rate_19(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.km*u[15]
end

function rate_20(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.a′*u[4]*u[7]/p.Cell_Params.Ω
end

function rate_21(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.d′*u[16]
end

function rate_22(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.km*u[16]
end

function rate_23(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.κ*u[5]
end

function rate_24(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.κ*u[6]
end

function rate_25(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.δ*u[3]
end

function rate_26(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.δ*u[4]
end

function rate_27(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.δ*u[5]
end

function rate_28(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.δ*u[6]
end

function rate_29(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.δ*u[7]
end

function rate_30(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.δ*u[8]
end

function rate_31(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.γ*u[9]
end

function rate_32(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.γ*u[10]
end

function rate_33(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.a*u[17]*u[9]/p.Cell_Params.Ω
end

function rate_34(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.d*u[18]
end

function rate_35(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.a*u[19]*u[10]/p.Cell_Params.Ω
end

function rate_36(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.d*u[20]
end

# Extents
function rate_1_affect!(integrator)
    integrator.u[3] += 1
end

function rate_2_affect!(integrator)
    integrator.u[4] += 1
end

function rate_3_affect!(integrator)
    integrator.u[1] -= 1
    integrator.u[9] -= 1
    integrator.u[11] += 1
end

function rate_4_affect!(integrator)
    integrator.u[1] += 1
    integrator.u[9] += 1
    integrator.u[11] -= 1
end

function rate_5_affect!(integrator)
    integrator.u[3] += 1
end

function rate_6_affect!(integrator)
    integrator.u[2] -= 1
    integrator.u[10] -= 1
    integrator.u[12] += 1
end

function rate_7_affect!(integrator)
    integrator.u[2] += 1
    integrator.u[10] += 1
    integrator.u[12] -= 1
end

function rate_8_affect!(integrator)
    integrator.u[4] += 1
end

function rate_9_affect!(integrator)
    integrator.u[3] -= 1
    integrator.u[5] += 1
    integrator.u[7] += 1
end

function rate_10_affect!(integrator)
    integrator.u[4] -= 1
    integrator.u[6] += 1
    integrator.u[8] += 1
end

function rate_11_affect!(integrator)
    integrator.u[5] -= 1
    integrator.u[8] -= 1
    integrator.u[13] += 1
end

function rate_12_affect!(integrator)
    integrator.u[5] += 1
    integrator.u[8] += 1
    integrator.u[13] -= 1
end

function rate_13_affect!(integrator)
    integrator.u[8] += 1
    integrator.u[13] -= 1
end

function rate_14_affect!(integrator)
    integrator.u[6] -= 1
    integrator.u[7] -= 1
    integrator.u[14] += 1
end

function rate_15_affect!(integrator)
    integrator.u[6] += 1
    integrator.u[7] += 1
    integrator.u[14] -= 1
end

function rate_16_affect!(integrator)
    integrator.u[7] += 1
    integrator.u[14] -= 1
end

function rate_17_affect!(integrator)
    integrator.u[3] -= 1
    integrator.u[8] -= 1
    integrator.u[15] += 1
end

function rate_18_affect!(integrator)
    integrator.u[3] += 1
    integrator.u[8] += 1
    integrator.u[15] -= 1
end

function rate_19_affect!(integrator)
    integrator.u[8] += 1
    integrator.u[15] -= 1
end

function rate_20_affect!(integrator)
    integrator.u[4] -= 1
    integrator.u[7] -= 1
    integrator.u[16] += 1
end

function rate_21_affect!(integrator)
    integrator.u[4] += 1
    integrator.u[7] += 1
    integrator.u[16] -= 1
end

function rate_22_affect!(integrator)
    integrator.u[7] += 1
    integrator.u[16] -= 1
end

function rate_23_affect!(integrator)
    integrator.u[9] += 1
end

function rate_24_affect!(integrator)
    integrator.u[10] += 1
end

function rate_25_affect!(integrator)
    integrator.u[3] -= 1
end

function rate_26_affect!(integrator)
    integrator.u[4] -= 1
end

function rate_27_affect!(integrator)
    integrator.u[5] -= 1
end

function rate_28_affect!(integrator)
    integrator.u[6] -= 1
end

function rate_29_affect!(integrator)
    integrator.u[7] -= 1
end

function rate_30_affect!(integrator)
    integrator.u[8] -= 1
end

function rate_31_affect!(integrator)
    integrator.u[9] -= 1
end

function rate_32_affect!(integrator)
    integrator.u[10] -= 1
end

function rate_33_affect!(integrator)
    integrator.u[9] -= 1
    integrator.u[17] -= 1
    integrator.u[18] += 1
end

function rate_34_affect!(integrator)
    integrator.u[9] += 1
    integrator.u[17] += 1
    integrator.u[18] -= 1
end

function rate_35_affect!(integrator)
    integrator.u[10] -= 1
    integrator.u[19] -= 1
    integrator.u[20] += 1
end

function rate_36_affect!(integrator)
    integrator.u[10] += 1
    integrator.u[19] += 1
    integrator.u[20] -= 1
end

## Generate stochastic simulation problem
function generate_jump_problem(ic::Vector{Float64}, tspan::Tuple{Float64,Float64}, p::Simulation_Parameters)
    rxn_1_jump = ConstantRateJump(rate_1, rate_1_affect!)
    rxn_2_jump = ConstantRateJump(rate_2, rate_2_affect!)
    rxn_3_jump = ConstantRateJump(rate_3, rate_3_affect!)
    rxn_4_jump = ConstantRateJump(rate_4, rate_4_affect!)
    rxn_5_jump = ConstantRateJump(rate_5, rate_5_affect!)
    rxn_6_jump = ConstantRateJump(rate_6, rate_6_affect!)
    rxn_7_jump = ConstantRateJump(rate_7, rate_7_affect!)
    rxn_8_jump = ConstantRateJump(rate_8, rate_8_affect!)
    rxn_9_jump = ConstantRateJump(rate_9, rate_9_affect!)
    rxn_10_jump = ConstantRateJump(rate_10, rate_10_affect!)
    rxn_11_jump = ConstantRateJump(rate_11, rate_11_affect!)
    rxn_12_jump = ConstantRateJump(rate_12, rate_12_affect!)
    rxn_13_jump = ConstantRateJump(rate_13, rate_13_affect!)
    rxn_14_jump = ConstantRateJump(rate_14, rate_14_affect!)
    rxn_15_jump = ConstantRateJump(rate_15, rate_15_affect!)
    rxn_16_jump = ConstantRateJump(rate_16, rate_16_affect!)
    rxn_17_jump = ConstantRateJump(rate_17, rate_17_affect!)
    rxn_18_jump = ConstantRateJump(rate_18, rate_18_affect!)
    rxn_19_jump = ConstantRateJump(rate_19, rate_19_affect!)
    rxn_20_jump = ConstantRateJump(rate_20, rate_20_affect!)
    rxn_21_jump = ConstantRateJump(rate_21, rate_21_affect!)
    rxn_22_jump = ConstantRateJump(rate_22, rate_22_affect!)
    rxn_23_jump = ConstantRateJump(rate_23, rate_23_affect!)
    rxn_24_jump = ConstantRateJump(rate_24, rate_24_affect!)
    rxn_25_jump = ConstantRateJump(rate_25, rate_25_affect!)
    rxn_26_jump = ConstantRateJump(rate_26, rate_26_affect!)
    rxn_27_jump = ConstantRateJump(rate_27, rate_27_affect!)
    rxn_28_jump = ConstantRateJump(rate_28, rate_28_affect!)
    rxn_29_jump = ConstantRateJump(rate_29, rate_29_affect!)
    rxn_30_jump = ConstantRateJump(rate_30, rate_30_affect!)
    rxn_31_jump = ConstantRateJump(rate_31, rate_31_affect!)
    rxn_32_jump = ConstantRateJump(rate_32, rate_32_affect!)
    rxn_33_jump = ConstantRateJump(rate_33, rate_33_affect!)
    rxn_34_jump = ConstantRateJump(rate_34, rate_34_affect!)
    rxn_35_jump = ConstantRateJump(rate_33, rate_33_affect!)
    rxn_36_jump = ConstantRateJump(rate_34, rate_34_affect!)

    return JumpProblem(
        DiscreteProblem(ic,tspan,p), 
        Direct(),
        rxn_1_jump,rxn_2_jump,rxn_3_jump,rxn_4_jump,rxn_5_jump,rxn_6_jump,rxn_7_jump,rxn_8_jump,rxn_9_jump,rxn_10_jump,
        rxn_11_jump,rxn_12_jump,rxn_13_jump,rxn_14_jump,rxn_15_jump,rxn_16_jump,rxn_17_jump,rxn_18_jump,rxn_19_jump,rxn_20_jump,
        rxn_21_jump,rxn_22_jump,rxn_23_jump,rxn_24_jump,rxn_25_jump,rxn_26_jump,rxn_27_jump,rxn_28_jump,rxn_29_jump,rxn_30_jump,
        rxn_31_jump,rxn_32_jump,rxn_33_jump,rxn_34_jump,rxn_35_jump,rxn_36_jump
        )
end

end