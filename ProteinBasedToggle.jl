module ProteinBasedToggle

using DifferentialEquations.JumpProcesses

export generate_jump_problem
export TF_Toggle_Parameters, DEFAULT_TOGGLE_PARAMS
export Cell_Parameters, DEFAULT_CELL_PARAMS
export Simulation_Parameters, DEFAULT_SIM_PARAMS

mutable struct TF_Toggle_Parameters
    α0::Float64 #Basal transcription rate constant [1/hr]
    α::Float64 #Activated transcription rate constant [1/hr]
    a1::Float64 #Activator-DNA association rate constant TALE1 [1/num*hr]
    d1::Float64 #Activator-DNA dissociation rate constant TALE1 [1/hr]
    a2::Float64 #Activator-DNA association rate constant TALE2 [1/num*hr]
    d2::Float64 #Activator-DNA dissociation rate constant TALE2 [1/hr]
    κ::Float64 #Translation rate constant [1/hr]
    δ::Float64 #RNA degradation rate constant [1/hr]
    γ::Float64 #Protein degradation rate constant [1/hr]
end

DEFAULT_TOGGLE_PARAMS = TF_Toggle_Parameters(0.25,2,1,1,1,1,5,1,1)

mutable struct Cell_Parameters
    D1tot::Float64 #Node 1 DNA copy number
    D2tot::Float64 #Node 2 DNA copy number
    D3tot::Float64 #Downstream target DNA copy number
    Ω::Float64 #System volume
end

DEFAULT_CELL_PARAMS = Cell_Parameters(10,10,0,1)

mutable struct Simulation_Parameters
    Toggle_Params::TF_Toggle_Parameters
    Cell_Params::Cell_Parameters
end

DEFAULT_SIM_PARAMS = Simulation_Parameters(DEFAULT_TOGGLE_PARAMS,DEFAULT_CELL_PARAMS)
# u = [ u1 = D1
#       u2 = D2
#       u3 = m1
#       u4 = m2
#       u5 = R1 
#       u6 = R2 
#       u7 = A1 
#       u8 = A2 
#       u9 = C1
#       u10 = C2 
#       u11 = C3 
#       u12 = C4 
#       u13 = D3
#       u14 = C5
# ]



function rate_1(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.α0*u[1]
end

function rate_2(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.α0*u[2]
end

function rate_3(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.a1*u[7]*u[1]/p.Cell_Params.Ω
end

function rate_4(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.d1*u[9]
end

function rate_5(u, p::Simulation_Parameters, t)
    return p.Toggle_Params.α * u[9]
end

function rate_6(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.a2*u[8]*u[2]/p.Cell_Params.Ω
end

function rate_7(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.d2*u[10]
end

function rate_8(u, p::Simulation_Parameters, t)
    return p.Toggle_Params.α * u[10]
end

function rate_9(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.a1*u[5]*u[1]/p.Cell_Params.Ω
end

function rate_10(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.d1*u[11]
end

function rate_11(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.a2*u[6]*u[2]/p.Cell_Params.Ω
end

function rate_12(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.d2*u[12]
end

function rate_13(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.κ*u[3]
end

function rate_14(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.κ*u[4]
end

function rate_15(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.δ*u[3]
end

function rate_16(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.δ*u[4]
end

function rate_17(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.γ*u[5]
end

function rate_18(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.γ*u[6]
end

function rate_19(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.γ*u[7]
end

function rate_20(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.γ*u[8]
end

function rate_21(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.a1*u[7]*u[13]/p.Cell_Params.Ω
end

function rate_22(u,p::Simulation_Parameters,t)
    return p.Toggle_Params.d1*u[14]
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
    integrator.u[7] -= 1
    integrator.u[9] += 1
end

function rate_4_affect!(integrator)
    integrator.u[1] += 1
    integrator.u[7] += 1
    integrator.u[9] -= 1
end

function rate_5_affect!(integrator)
    integrator.u[3] += 1
end

function rate_6_affect!(integrator)
    integrator.u[2] -= 1
    integrator.u[8] -= 1
    integrator.u[10] += 1
end

function rate_7_affect!(integrator)
    integrator.u[2] += 1
    integrator.u[8] += 1
    integrator.u[10] -= 1
end

function rate_8_affect!(integrator)
    integrator.u[4] += 1
end

function rate_9_affect!(integrator)
    integrator.u[1] -= 1
    integrator.u[5] -= 1
    integrator.u[11] += 1
end

function rate_10_affect!(integrator)
    integrator.u[1] += 1
    integrator.u[5] += 1
    integrator.u[11] -= 1
end

function rate_11_affect!(integrator)
    integrator.u[2] -= 1
    integrator.u[6] -= 1
    integrator.u[12] += 1
end

function rate_12_affect!(integrator)
    integrator.u[2] += 1
    integrator.u[6] += 1
    integrator.u[12] -= 1
end

function rate_13_affect!(integrator)
    integrator.u[7] += 1
    integrator.u[6] += 1
end

function rate_14_affect!(integrator)
    integrator.u[8] += 1
    integrator.u[5] += 1
end

function rate_15_affect!(integrator)
    integrator.u[3] -= 1
end

function rate_16_affect!(integrator)
    integrator.u[4] -= 1
end

function rate_17_affect!(integrator)
    integrator.u[5] -= 1
end

function rate_18_affect!(integrator)
    integrator.u[6] -= 1
end

function rate_19_affect!(integrator)
    integrator.u[7] -= 1
end

function rate_20_affect!(integrator)
    integrator.u[8] -= 1
end

function rate_21_affect!(integrator)
    integrator.u[13] -= 1
    integrator.u[7] -= 1
    integrator.u[14] += 1
end

function rate_22_affect!(integrator)
    integrator.u[13] += 1
    integrator.u[7] += 1
    integrator.u[14] -= 1
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

    return JumpProblem(
        DiscreteProblem(ic,tspan,p), 
        Direct(),
        rxn_1_jump,rxn_2_jump,rxn_3_jump,rxn_4_jump,rxn_5_jump,rxn_6_jump,rxn_7_jump,rxn_8_jump,rxn_9_jump,rxn_10_jump,
        rxn_11_jump,rxn_12_jump,rxn_13_jump,rxn_14_jump,rxn_15_jump,rxn_16_jump,rxn_17_jump,rxn_18_jump,rxn_19_jump,rxn_20_jump,
        rxn_21_jump,rxn_22_jump
        )
end
end