module ProteinToggleTest

using DifferentialEquations
using DifferentialEquations.JumpProcesses
using Plots
using DifferentialEquations.EnsembleAnalysis
include("ProteinBasedToggle.jl")
using .ProteinBasedToggle
using Peaks
using KernelDensity
using StatsBase
# Initial testing to look at single trajectories, no D3

p = DEFAULT_SIM_PARAMS
ic = [DEFAULT_CELL_PARAMS.D1tot, #D1
      DEFAULT_CELL_PARAMS.D2tot, #D2 
      0,                         #m1
      0,                         #m2
      0,                         #R1
      0,                         #R2
      10,                         #A1
      0,                         #A2
      0,                         #C1 
      0,                         #C2
      0,                         #C3 
      0,                         #C4
      DEFAULT_CELL_PARAMS.D3tot, #D3
      0]                         #C5
tspan = (0.0,50.0)
jump_prob  = generate_jump_problem(ic, tspan, p)
sol = solve(jump_prob, SSAStepper())
plot_1 = plot(sol, idxs = [(0, 7), (0, 8)], 
    label =["A1" "A2"], lw = 3, lc = [:black :red ])
xlabel!("time [hr]")
ylabel!("Protein counts")
savefig(plot_1, "figures/sample_traj.svg")


#-------------------------------------------------
#100 trajectories starting with A1 = 30 (A1 high state), find time when A1 = A2
#as the time to first "State switch"


p = DEFAULT_SIM_PARAMS
ic = [DEFAULT_CELL_PARAMS.D1tot, #D1
      DEFAULT_CELL_PARAMS.D2tot, #D2 
      0,                         #m1
      0,                         #m2
      0,                         #R1
      0,                         #R2
      70.0,                         #A1
      0,                         #A2
      0,                         #C1 
      0,                         #C2
      0,                         #C3 
      0,                         #C4
      DEFAULT_CELL_PARAMS.D3tot, #D3
      0]                         #C5
tspan = (0.0,100.0)

function survival_plot_output_func(sol, i)
    a1 = [v[7] for v in sol.u]
    a2 = [v[8] for v in sol.u]
    switch_index = findfirst((a1-a2) .==0)
    if switch_index == nothing
        switch_time = tspan[end]
    else
        switch_time = sol.t[switch_index]
    end
    return (switch_time, false)
end


trajectories_num = 100
α_vec = [0.2, 2.0, 4.0]
D3_vec = [0.0, 30.0]
α_mat = vcat(α_vec', α_vec')

for i = 1:lastindex(α_vec)
    p.Toggle_Params.α = α_vec[i]
    for k = 1:lastindex(D3_vec)
        ic[13] = D3_vec[k]
        jump_prob  = generate_jump_problem(ic, tspan, p)
        ens_prob = EnsembleProblem(prob = jump_prob, output_func = survival_plot_output_func)
        sim = solve(ens_prob, SSAStepper(), EnsembleThreads(), trajectories=trajectories_num)

        time_vec = sort(sim.u)
        num_of_switching_traj = count(sim.u .< tspan[end])
        if num_of_switching_traj == trajectories_num
            y = (trajectories_num .- collect(0:trajectories_num)) ./ trajectories_num 
            time_vec = append!(zeros(1), time_vec)
        else
            time_vec = time_vec[1:num_of_switching_traj+1]
            
            time_vec = append!(zeros(1), time_vec)
            counter = collect(0:num_of_switching_traj)
            counter = append!(counter, [num_of_switching_traj])
            y = (trajectories_num .- counter) ./ trajectories_num
        end
    
        if i == 1
            if k==1
                decay_plot = plot(time_vec, y, line_z = α_vec[i] * ones(length(time_vec)), label = "", zcolor = α_vec, grid = false, c= cgrad(:viridis), colorbar_title = "α", ls=:dash)
                ylims!((0, 1))
                xlabel!("time")
                ylabel!("Fraction of simulations in\noriginal A1-high state")
            else
                plot!(time_vec, y, line_z = α_vec[i] * ones(length(time_vec)), label = "", zcolor = α_vec, grid = false, c= cgrad(:viridis), colorbar_title = "α")
            end
        else
            if k==1
                plot!(time_vec, y, line_z = α_vec[i] * ones(length(time_vec)), label = "", zcolor = α_vec, grid = false, c= cgrad(:viridis), colorbar_title = "α", ls=:dash)
            else
                plot!(time_vec, y, line_z = α_vec[i] * ones(length(time_vec)), label = "", zcolor = α_vec, grid = false, c= cgrad(:viridis), colorbar_title = "α")
            end
        end
    end
end
current()
savefig(current(), "figures/D3vsalpha_switching_fixed.svg")

p = DEFAULT_SIM_PARAMS
ic = [DEFAULT_CELL_PARAMS.D1tot, #D1
      DEFAULT_CELL_PARAMS.D2tot, #D2 
      0,                         #m1
      0,                         #m2
      0,                         #R1
      0,                         #R2
      0,                         #A1
      0,                         #A2
      0,                         #C1 
      0,                         #C2
      0,                         #C3 
      0,                         #C4
      DEFAULT_CELL_PARAMS.D3tot, #D3
      0]                         #C5
tspan = (0.0,100.0)
trajectories_num = 500
#α_vec = [2.0]
α_vec = [0.2, 2.0, 4.0]
D3_vec = [0.0, 1.0, 10.0, 30.0]
p = DEFAULT_SIM_PARAMS
a1_low_states = zeros(lastindex(D3_vec), lastindex(α_vec))
a1_high_states = zeros(lastindex(D3_vec), lastindex(α_vec))
cv_A1_high = zeros(lastindex(D3_vec), lastindex(α_vec))
frac_A1_high = zeros(lastindex(D3_vec), lastindex(α_vec))


function ss_output_func(sol, i)
    A1_end = sol[end][7]
    return (A1_end, false)
end

quant_color_pal = cgrad(:viridis, length(α_vec), categorical = true)

for i = 1:lastindex(α_vec)
    A1_vec = zeros(trajectories_num)
    p.Toggle_Params.α = α_vec[i]

    for k = 1:lastindex(D3_vec)
        ic[13] = D3_vec[k]    
        jump_prob  = generate_jump_problem(ic, tspan, p)
        ens_prob = EnsembleProblem(prob = jump_prob, output_func = ss_output_func)
        sim = solve(ens_prob, SSAStepper(), EnsembleThreads(), trajectories=trajectories_num)
        A1_vec = sim.u
        a1_kde = kde(A1_vec)
        peak_idxs, _ = findmaxima(a1_kde.density)
        if lastindex(peak_idxs) == 2
            a1_high_states[k, i] = a1_kde.x[peak_idxs[end]] 
            a1_low_states[k, i] = a1_kde.x[peak_idxs[1]]
        else
            a1_low_states[k, i] = a1_kde.x[peak_idxs[1]]
            a1_high_states[k, i] = a1_kde.x[peak_idxs[1]] 
        end

        A1_vec_A1_high = A1_vec[A1_vec.>10]
        cv_A1_high[k, i] = std(A1_vec_A1_high) / mean(A1_vec_A1_high)
        frac_A1_high[k, i] = length(A1_vec_A1_high) / length(A1_vec)
        
        hist1 = histogram(A1_vec, bins =20,  normalize = :pdf, legend = :none, fillcolor = :gray, lc = :gray ,fillalpha= 0.5, grid=false)
        plot!(a1_kde.x, a1_kde.density, lc = :black, lw = 2)
        xlims!((0, extrema(A1_vec)[2]))
        xlabel!("A1 [protein counts]")
        ylabel!("Probability Density")
        savefig(hist1, string("figures/A1_histogram_base_params_alpha_", p.Toggle_Params.α, 
                                "D3tot_", D3_vec[k] , "_.svg"))
    end
end


α_mat = vcat(α_vec', α_vec', α_vec', α_vec')
p1 = plot(D3_vec, a1_high_states, line_z = α_mat, label = "", zcolor = α_vec, grid = false, c= cgrad(:viridis), colorbar_title = "α" )
xlabel!("D3 total")
ylabel!("A1 Protein counts in A1-high state")
savefig(p1, "figures/D3vsalpha_A1_counts.svg")

p2 = plot(D3_vec, frac_A1_high, line_z = α_mat, label = "", zcolor = α_vec, grid = false, c= cgrad(:viridis) , colorbar_title = "α")
xlabel!("D3 total")
ylabel!("Fraction A1-high cells")
savefig(p2, "figures/D3vsalpha_A1_frac.svg")

p3 = plot(D3_vec, cv_A1_high, line_z = α_mat, label = "", zcolor = α_vec, grid = false, c= cgrad(:viridis), colorbar_title = "α" )
xlabel!("D3 total")
ylabel!("Coefficient of Variation\nof A1-high cells")
savefig(p3, "figures/D3vsalpha_A1_cv.svg")
display(p2)


#-------------------------------------------------
# Moving on from alpha and D3
# look atcopy number of D1, D2

#100 trajectories starting with A1 = 30 (A1 high state), find time when A1 = A2
#as the time to first "State switch"


p = DEFAULT_SIM_PARAMS
ic = [DEFAULT_CELL_PARAMS.D1tot, #D1
      DEFAULT_CELL_PARAMS.D2tot, #D2 
      0,                         #m1
      0,                         #m2
      0,                         #R1
      70.0,                         #R2
      70.0,                         #A1
      0,                         #A2
      0,                         #C1 
      0,                         #C2
      0,                         #C3 
      0,                         #C4
      DEFAULT_CELL_PARAMS.D3tot, #D3
      0]                         #C5
tspan = (0.0,100.0)


trajectories_num = 100
D1_D2_vec = [1.0, 3.0, 10.0, 30.0]
D3_vec = [0.0, 30.0]

for i = 1:lastindex(D1_D2_vec)
    ic[1] = D1_D2_vec[i]
    ic[2] = D1_D2_vec[i]
    for k = 1:lastindex(D3_vec)
        ic[13] = D3_vec[k]
        num_of_switching_traj = 0
        time_vec = zeros(trajectories_num)
        jump_prob  = generate_jump_problem(ic, tspan, p)
        ens_prob = EnsembleProblem(prob = jump_prob)
        sim = solve(ens_prob, SSAStepper(), EnsembleThreads(), trajectories=trajectories_num)
        for j = 1:lastindex(sim.u)
            a1 = [v[7] for v in sim[j].u]
            a2 = [v[8] for v in sim[j].u]
            switch_index = findfirst((a1-a2) .==0)
            if switch_index == nothing
                time_vec[j] = tspan[end]
            else
                time_vec[j] = sim[j].t[switch_index]
                num_of_switching_traj = num_of_switching_traj +1
            end
        end
        time_vec = sort(time_vec)
        if num_of_switching_traj == trajectories_num
            y = (trajectories_num .- collect(0:trajectories_num)) ./ trajectories_num 
            time_vec = append!(zeros(1), time_vec)
        else
            time_vec = time_vec[1:num_of_switching_traj+1]
            
            time_vec = append!(zeros(1), time_vec)
            counter = collect(0:num_of_switching_traj)
            counter = append!(counter, [num_of_switching_traj])
            y = (trajectories_num .- counter) ./ trajectories_num
        end
    
        if i == 1
            if k==1
                decay_plot = plot(time_vec, y, line_z = D1_D2_vec[i] * ones(length(time_vec)), label = "", zcolor = D1_D2_vec, grid = false, c= cgrad(:viridis), colorbar_title = "D1, D2 copies", ls=:dash)
                ylims!((0, 1))
                xlabel!("time")
                ylabel!("Fraction of simulations in\noriginal A1-high state")
            else
                plot!(time_vec, y, line_z = D1_D2_vec[i] * ones(length(time_vec)), label = "", zcolor = D1_D2_vec, grid = false, c= cgrad(:viridis), colorbar_title = "D1, D2 copies")
            end
        else
            if k==1
                plot!(time_vec, y, line_z = D1_D2_vec[i] * ones(length(time_vec)), label = "", zcolor = D1_D2_vec, grid = false, c= cgrad(:viridis), colorbar_title = "D1, D2 copies", ls=:dash)
            else
                plot!(time_vec, y, line_z = D1_D2_vec[i] * ones(length(time_vec)), label = "", zcolor = D1_D2_vec, grid = false, c= cgrad(:viridis), colorbar_title = "D1, D2 copies")
            end
        end
    end
end
current()
savefig(current(), "figures/D3vsD1_switching_fixed.svg")

p = DEFAULT_SIM_PARAMS
ic = [DEFAULT_CELL_PARAMS.D1tot, #D1
      DEFAULT_CELL_PARAMS.D2tot, #D2 
      0,                         #m1
      0,                         #m2
      0,                         #R1
      0,                         #R2
      0,                         #A1
      0,                         #A2
      0,                         #C1 
      0,                         #C2
      0,                         #C3 
      0,                         #C4
      DEFAULT_CELL_PARAMS.D3tot, #D3
      0]                         #C5
tspan = (0.0,100.0)
trajectories_num = 500
D1_D2_vec = [1.0, 3.0, 10.0, 30.0]
D3_vec = [0.0, 1.0, 10.0, 30.0]

p = DEFAULT_SIM_PARAMS
a1_low_states = zeros(lastindex(D3_vec), lastindex(D1_D2_vec))
a1_high_states = zeros(lastindex(D3_vec), lastindex(D1_D2_vec))
cv_A1_high = zeros(lastindex(D3_vec), lastindex(D1_D2_vec))
frac_A1_high = zeros(lastindex(D3_vec), lastindex(D1_D2_vec))

for i = 1:lastindex(D1_D2_vec)
    ic[1] = D1_D2_vec[i]
    ic[2] = D1_D2_vec[i]
    A1_vec = zeros(trajectories_num)
    for k = 1:lastindex(D3_vec)
        ic[13] = D3_vec[k]    
        jump_prob  = generate_jump_problem(ic, tspan, p)
        ens_prob = EnsembleProblem(prob = jump_prob)
        sim = solve(ens_prob, SSAStepper(), EnsembleThreads(), trajectories=trajectories_num)
        for j = 1:lastindex(sim.u)
            A1_vec[j] = sim[j].u[end][7]
        end

        a1_kde = kde(A1_vec)
        peak_idxs, _ = findmaxima(a1_kde.density)
        if lastindex(peak_idxs) == 2
            a1_high_states[k, i] = a1_kde.x[peak_idxs[end]] 
            a1_low_states[k, i] = a1_kde.x[peak_idxs[1]]
        else
            a1_low_states[k, i] = a1_kde.x[peak_idxs[1]]
            a1_high_states[k, i] = a1_kde.x[peak_idxs[1]] 
        end

        A1_vec_A1_high = A1_vec[A1_vec.>10]
        cv_A1_high[k, i] = std(A1_vec_A1_high) / mean(A1_vec_A1_high)
        frac_A1_high[k, i] = length(A1_vec_A1_high) / length(A1_vec)
        
        hist1 = histogram(A1_vec, bins =20,  normalize = :pdf, legend = :none, fillcolor = :gray, lc = :gray ,fillalpha= 0.5, grid=false)
        plot!(a1_kde.x, a1_kde.density, lc = :black, lw = 2)
        xlims!((0, extrema(A1_vec)[2]))
        xlabel!("A1 [protein counts]")
        ylabel!("Probability Density")
        savefig(hist1, string("figures/A1_histogram_base_params_D1_", D1_D2_vec[i], 
                                "D3tot_", D3_vec[k] , "_.svg"))
    end
    # if i == 1
    #     plot(D3_vec, a1_high_states[:, i], label = [string("α=", α_vec[i])], palette= quant_color_pal )
    #     plot!(D3_vec, a1_low_states[:, i], legend=:none, palette= quant_color_pal)
    # else
    #     plot!(D3_vec, a1_high_states[:, i], label = [string("α=", α_vec[i])], palette= quant_color_pal )
    #     plot!(D3_vec, a1_low_states[:, i], legend=:none, palette= quant_color_pal)
    # end
end


D1_D2_mat = vcat(D1_D2_vec', D1_D2_vec', D1_D2_vec', D1_D2_vec')
p1 = plot(D3_vec, a1_high_states, line_z = D1_D2_mat, label = "", zcolor = D1_D2_vec, grid = false, c= cgrad(:viridis), colorbar_title = "D1, D2" )
xlabel!("D3 total")
ylabel!("A1 Protein counts in A1-high state")
savefig(p1, "figures/D3vsD1_A1_counts.svg")

p2 = plot(D3_vec, frac_A1_high, line_z = D1_D2_mat, label = "", zcolor = D1_D2_vec, grid = false, c= cgrad(:viridis) , colorbar_title = "D1, D2")
xlabel!("D3 total")
ylabel!("Fraction A1-high cells")
savefig(p2, "figures/D3vsD1_A1_frac.svg")

p3 =plot(D3_vec, cv_A1_high, line_z = D1_D2_mat, label = "", zcolor = D1_D2_vec, grid = false, c= cgrad(:viridis), colorbar_title = "D1, D2" )
xlabel!("D3 total")
ylabel!("Coefficient of Variation\nof A1-high cells")
savefig(p3, "figures/D3vsD1_A1_cv.svg")

#-------------------------------------------------
# Moving on from alpha and D3
# look atcopy number of D1, D2

#100 trajectories starting with A1 = 30 (A1 high state), find time when A1 = A2
#as the time to first "State switch"


p = DEFAULT_SIM_PARAMS
ic = [DEFAULT_CELL_PARAMS.D1tot, #D1
      DEFAULT_CELL_PARAMS.D2tot, #D2 
      0,                         #m1
      0,                         #m2
      0,                         #R1
      0,                         #R2
      70.0,                         #A1
      0,                         #A2
      0,                         #C1 
      0,                         #C2
      0,                         #C3 
      0,                         #C4
      DEFAULT_CELL_PARAMS.D3tot, #D3
      0]                         #C5
tspan = (0.0,100.0)


trajectories_num = 100
kappa_vec = [5.0, 2.0, 1.0]
gamma_vec = [1.0, 0.4, 0.2]
D3_vec = [0.0, 30.0]

for i = 1:lastindex(kappa_vec)
    p.Toggle_Params.κ = kappa_vec[i]
    p.Toggle_Params.γ = gamma_vec[i]
    for k = 1:lastindex(D3_vec)
        ic[13] = D3_vec[k]
        num_of_switching_traj = 0
        time_vec = zeros(trajectories_num)
        jump_prob  = generate_jump_problem(ic, tspan, p)
        ens_prob = EnsembleProblem(prob = jump_prob)
        sim = solve(ens_prob, SSAStepper(), EnsembleThreads(), trajectories=trajectories_num)
        for j = 1:lastindex(sim.u)
            a1 = [v[7] for v in sim[j].u]
            a2 = [v[8] for v in sim[j].u]
            switch_index = findfirst((a1-a2) .==0)
            if switch_index == nothing
                time_vec[j] = tspan[end]
            else
                time_vec[j] = sim[j].t[switch_index]
                num_of_switching_traj = num_of_switching_traj +1
            end
        end
        time_vec = sort(time_vec)
        if num_of_switching_traj == trajectories_num
            y = (trajectories_num .- collect(0:trajectories_num)) ./ trajectories_num 
            time_vec = append!(zeros(1), time_vec)
        else
            time_vec = time_vec[1:num_of_switching_traj+1]
            
            time_vec = append!(zeros(1), time_vec)
            counter = collect(0:num_of_switching_traj)
            counter = append!(counter, [num_of_switching_traj])
            y = (trajectories_num .- counter) ./ trajectories_num
        end
    
        if i == 1
            if k==1
                decay_plot = plot(time_vec, y, line_z = 2.0 / kappa_vec[i] * ones(length(time_vec)), label = "", zcolor = 2.0 ./ kappa_vec, grid = false, c= cgrad(:viridis), colorbar_title = "α/κ", ls=:dash)
                ylims!((0, 1))
                xlabel!("time")
                ylabel!("Fraction of simulations in\noriginal A1-high state")
            else
                plot!(time_vec, y, line_z = 2.0 / kappa_vec[i] * ones(length(time_vec)), label = "", zcolor = 2.0 ./ kappa_vec, grid = false, c= cgrad(:viridis), colorbar_title = "α/κ")
            end
        else
            if k==1
                plot!(time_vec, y, line_z = 2.0 / kappa_vec[i] * ones(length(time_vec)), label = "", zcolor = 2.0 ./ kappa_vec, grid = false, c= cgrad(:viridis), colorbar_title = "α/κ", ls=:dash)
            else
                plot!(time_vec, y, line_z = 2.0 / kappa_vec[i] * ones(length(time_vec)), label = "", zcolor = 2.0 ./ kappa_vec, grid = false, c= cgrad(:viridis), colorbar_title = "α/κ")
            end
        end
    end
end
current()
savefig(current(), "figures/D3vskappa_switching_fixed.svg")

p = DEFAULT_SIM_PARAMS
ic = [DEFAULT_CELL_PARAMS.D1tot, #D1
      DEFAULT_CELL_PARAMS.D2tot, #D2 
      0,                         #m1
      0,                         #m2
      0,                         #R1
      0,                         #R2
      0,                         #A1
      0,                         #A2
      0,                         #C1 
      0,                         #C2
      0,                         #C3 
      0,                         #C4
      DEFAULT_CELL_PARAMS.D3tot, #D3
      0]                         #C5
tspan = (0.0,100.0)
trajectories_num = 500
kappa_vec = [5.0, 2.0, 1.0]
gamma_vec = [1.0, 0.4, 0.2]
D3_vec = [0.0, 1.0, 10.0, 30.0]

p = DEFAULT_SIM_PARAMS
a1_low_states = zeros(lastindex(D3_vec), lastindex(kappa_vec))
a1_high_states = zeros(lastindex(D3_vec), lastindex(kappa_vec))
cv_A1_high = zeros(lastindex(D3_vec), lastindex(kappa_vec))
frac_A1_high = zeros(lastindex(D3_vec), lastindex(kappa_vec))

for i = 1:lastindex(kappa_vec)
    p.Toggle_Params.κ = kappa_vec[i]
    p.Toggle_Params.γ = gamma_vec[i]
    A1_vec = zeros(trajectories_num)
    for k = 1:lastindex(D3_vec)
        ic[13] = D3_vec[k]    
        jump_prob  = generate_jump_problem(ic, tspan, p)
        ens_prob = EnsembleProblem(prob = jump_prob)
        sim = solve(ens_prob, SSAStepper(), EnsembleThreads(), trajectories=trajectories_num)
        for j = 1:lastindex(sim.u)
            A1_vec[j] = sim[j].u[end][7]
        end

        a1_kde = kde(A1_vec)
        peak_idxs, _ = findmaxima(a1_kde.density)
        if lastindex(peak_idxs) == 2
            a1_high_states[k, i] = a1_kde.x[peak_idxs[end]] 
            a1_low_states[k, i] = a1_kde.x[peak_idxs[1]]
        else
            a1_low_states[k, i] = a1_kde.x[peak_idxs[1]]
            a1_high_states[k, i] = a1_kde.x[peak_idxs[1]] 
        end

        A1_vec_A1_high = A1_vec[A1_vec.>10]
        cv_A1_high[k, i] = std(A1_vec_A1_high) / mean(A1_vec_A1_high)
        frac_A1_high[k, i] = length(A1_vec_A1_high) / length(A1_vec)
        
        hist1 = histogram(A1_vec, bins =20,  normalize = :pdf, legend = :none, fillcolor = :gray, lc = :gray ,fillalpha= 0.5, grid=false)
        plot!(a1_kde.x, a1_kde.density, lc = :black, lw = 2)
        xlims!((0, extrema(A1_vec)[2]))
        xlabel!("A1 [protein counts]")
        ylabel!("Probability Density")
        savefig(hist1, string("figures/A1_histogram_base_params_kappa_", kappa_vec[i], 
                                "D3tot_", D3_vec[k] , "_.svg"))
    end
    # if i == 1
    #     plot(D3_vec, a1_high_states[:, i], label = [string("α=", α_vec[i])], palette= quant_color_pal )
    #     plot!(D3_vec, a1_low_states[:, i], legend=:none, palette= quant_color_pal)
    # else
    #     plot!(D3_vec, a1_high_states[:, i], label = [string("α=", α_vec[i])], palette= quant_color_pal )
    #     plot!(D3_vec, a1_low_states[:, i], legend=:none, palette= quant_color_pal)
    # end
end


kappa_mat = 2.0 ./ vcat(kappa_vec', kappa_vec', kappa_vec', kappa_vec')
p1 = plot(D3_vec, a1_high_states, line_z = kappa_mat, label = "", zcolor = 2.0 ./ kappa_vec, grid = false, c= cgrad(:viridis), colorbar_title = "α/κ" )
xlabel!("D3 total")
ylabel!("A1 Protein counts in A1-high state")
savefig(p1, "figures/D3vskappa_A1_counts.svg")

p2 = plot(D3_vec, frac_A1_high, line_z = kappa_mat, label = "", zcolor = 2.0 ./ kappa_vec, grid = false, c= cgrad(:viridis) , colorbar_title = "α/κ")
xlabel!("D3 total")
ylabel!("Fraction A1-high cells")
savefig(p2, "figures/D3vskappa_A1_frac.svg")

p3 = plot(D3_vec, cv_A1_high, line_z = kappa_mat, label = "", zcolor = 2.0 ./ kappa_vec, grid = false, c= cgrad(:viridis), colorbar_title = "α/κ" )
xlabel!("D3 total")
ylabel!("Coefficient of Variation\nof A1-high cells")
savefig(p3, "figures/D3vskappa_A1_cv.svg")


end