## Load in packages
using DifferentialEquations
using DiffEqJump
using .ToggleDynamicModeling
using Plots
using Peaks
using KernelDensity
using StatsBase

function survival_plot_output_func(sol, i)
    a1 = [v[9] for v in sol.u]
    a2 = [v[10] for v in sol.u]
    switch_index = findfirst((a1-a2) .==0)
    if switch_index == nothing
        switch_time = tspan[end]
    else
        switch_time = sol.t[switch_index]
    end
    return (switch_time, false)
end

# Initial testing to look at single trajectories, no D3
p = DEFAULT_SIM_PARAMS
ic = [DEFAULT_CELL_PARAMS.D1tot, DEFAULT_CELL_PARAMS.D2tot, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, DEFAULT_CELL_PARAMS.D3tot, 0, 0, 0]
tspan = (0.0,500.0)
jump_prob  = generate_jump_problem(ic, tspan, p)
sol = solve(jump_prob, SSAStepper())
plot_1 = plot(sol.t,sol[9,:],label="A1", lw = 3, lc = :black)
plot!(sol.t,sol[10,:],label="A2", lw = 3, lc = :red)
xlabel!("time [hr]")
ylabel!("Protein counts")


#-------------------------------------------------
#100 trajectories starting with A1 = 20 (A1 high state), find time when A1 = A2
#as the time to first "State switch"

p = DEFAULT_SIM_PARAMS
ic = [DEFAULT_CELL_PARAMS.D1tot, DEFAULT_CELL_PARAMS.D2tot, 0, 0, 0, 0,
    0, 0, 70, 0, 0, 0,
    0, 0, 0, 0, DEFAULT_CELL_PARAMS.D3tot, 0, 0, 0]
tspan = (0.0,80.0)

trajectories_num = 100
α_vec = [1.0, 2.0]
#α_vec = [0.2, 2.0, 5.0]
time_vec = zeros(trajectories_num)
a1_low_states = zeros(0)
a1_high_states = zeros(0)
cv_A1_high = zeros(0)
frac_A1_high = zeros(0)
A1_vec = zeros(trajectories_num)

α_vec = [0.2, 2.0, 4.0]
D3_vec = [0.0, 10.0]
α_mat = vcat(α_vec', α_vec')

for i = 1:lastindex(α_vec)
    p.Toggle_Params.α = α_vec[i]
    for k = 1:lastindex(D3_vec)
        ic[17] = D3_vec[k]
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

savefig(current(), "jl_output/2023.05.10_toggle/D3vsalpha_switching_ic70.svg")

## D3 parameter sweep
p = DEFAULT_SIM_PARAMS
ic = [DEFAULT_CELL_PARAMS.D1tot, DEFAULT_CELL_PARAMS.D2tot, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, DEFAULT_CELL_PARAMS.D3tot, 0, 0, 0]
tspan = (0.0,30.0)
trajectories_num = 200
#α_vec = [2.0]
α_vec = [0.2, 2.0, 4.0]
D3_vec = [0.0, 1.0, 10.0, 30.0]
a1_low_states = zeros(lastindex(D3_vec), lastindex(α_vec))
a1_high_states = zeros(lastindex(D3_vec), lastindex(α_vec))
cv_A1_high = zeros(lastindex(D3_vec), lastindex(α_vec))
frac_A1_high = zeros(lastindex(D3_vec), lastindex(α_vec))

quant_color_pal = cgrad(:viridis, length(α_vec), categorical = true)

for i = 1:lastindex(α_vec)
    A1_vec = zeros(trajectories_num)
    p.Toggle_Params.α = α_vec[i]

    for k = 1:lastindex(D3_vec)
        ic[17] = D3_vec[k]    
        jump_prob  = generate_jump_problem(ic, tspan, p)
        ens_prob = EnsembleProblem(prob = jump_prob)
        sim = solve(ens_prob, SSAStepper(), EnsembleThreads(), trajectories=trajectories_num)
        for j = 1:lastindex(sim.u)
            A1_vec[j] = sim[j].u[end][9]
        end

        a1_kde = kde(A1_vec)
        peak_idxs, _ = findmaxima(a1_kde.density)
        if lastindex(peak_idxs) > 1
            a1_high_states[k, i] = a1_kde.x[peak_idxs[end]] 
            a1_low_states[k, i] = a1_kde.x[peak_idxs[1]]
        else
            a1_low_states[k, i] = a1_kde.x[peak_idxs[1]]
            a1_high_states[k, i] = a1_kde.x[peak_idxs[1]] 
        end

        A1_vec_A1_high = A1_vec[A1_vec.>10]
        cv_A1_high[k, i] = std(A1_vec_A1_high) / mean(A1_vec_A1_high)
        frac_A1_high[k, i] = length(A1_vec_A1_high) / length(A1_vec)
        
        hist1 = histogram(A1_vec, bins =10,  normalize = :pdf, legend = :none, fillcolor = :gray, lc = :gray ,fillalpha= 0.5, grid=false)
        plot!(a1_kde.x, a1_kde.density, lc = :black, lw = 2)
        xlims!((0, extrema(A1_vec)[2]))
        xlabel!("A1 [protein counts]")
        ylabel!("Probability Density")
        savefig(hist1, string("jl_output/2023.05.10_toggle/A1_histogram_base_params_alpha_", p.Toggle_Params.α, 
                                "D3tot_", D3_vec[k] , "_short.svg"))
    end
    # if i == 1
    #     plot(D3_vec, a1_high_states[:, i], label = [string("α=", α_vec[i])], palette= quant_color_pal )
    #     plot!(D3_vec, a1_low_states[:, i], legend=:none, palette= quant_color_pal)
    # else
    #     plot!(D3_vec, a1_high_states[:, i], label = [string("α=", α_vec[i])], palette= quant_color_pal )
    #     plot!(D3_vec, a1_low_states[:, i], legend=:none, palette= quant_color_pal)
    # end
end

qcp = cgrad(:viridis, length(α_vec), categorical = true)
plot(D3_vec, a1_high_states, label = ["α=0.2" "α=2" "α=4"], lc = [qcp[1] qcp[2] qcp[3]] , grid = false)
quant_color_pal = cgrad(:viridis, length(α_vec), categorical = true)
plot!(D3_vec, a1_low_states, label = ["" "" "" ], lc = [qcp[1] qcp[2] qcp[3]], ls=:dash )

α_mat = vcat(α_vec', α_vec', α_vec', α_vec')
p1 = plot(D3_vec, a1_high_states, line_z = α_mat, label = "", zcolor = α_vec, grid = false, c= cgrad(:viridis), colorbar_title = "α" )
xlabel!("D3 total")
ylabel!("A1 Protein counts in A1-high state")
savefig(p1, "jl_output/2023.05.10_toggle/D3vsalpha_A1_counts_short.svg")

p2 = plot(D3_vec, frac_A1_high, line_z = α_mat, label = "", zcolor = α_vec, grid = false, c= cgrad(:viridis) , colorbar_title = "α")
xlabel!("D3 total")
ylabel!("Fraction A1-high cells")
savefig(p2, "jl_output/2023.05.10_toggle/D3vsalpha_A1_frac_short.svg")

p3 = plot(D3_vec, cv_A1_high, line_z = α_mat, label = "", zcolor = α_vec, grid = false, c= cgrad(:viridis), colorbar_title = "α" )
xlabel!("D3 total")
ylabel!("Coefficient of Variation\nof A1-high cells")
savefig(p3, "jl_output/2023.05.10_toggle/D3vsalpha_A1_cv_short.svg")

##

#-------------------------------------------------
# Moving on from alpha and D3
# look atcopy number of D1, D2

#100 trajectories starting with A1 = 30 (A1 high state), find time when A1 = A2
#as the time to first "State switch"

p = DEFAULT_SIM_PARAMS
ic = [DEFAULT_CELL_PARAMS.D1tot, DEFAULT_CELL_PARAMS.D2tot, 0, 0, 0, 0,
    0, 0, 70, 0, 0, 0,
    0, 0, 0, 0, DEFAULT_CELL_PARAMS.D3tot, 0, 0, 0]
tspan = (0.0,80.0)

trajectories_num = 100
D1_D2_vec = [1.0, 3.0, 10.0, 30.0]
D3_vec = [0.0, 30.0]

for i = 1:lastindex(D1_D2_vec)
    ic[1] = D1_D2_vec[i]
    ic[2] = D1_D2_vec[i]
    for k = 1:lastindex(D3_vec)
        ic[17] = D3_vec[k]
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

savefig(current(), "jl_output/2023.05.10_toggle/D3vsD1_switching.svg")


p = DEFAULT_SIM_PARAMS
ic = [DEFAULT_CELL_PARAMS.D1tot, DEFAULT_CELL_PARAMS.D2tot, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, DEFAULT_CELL_PARAMS.D3tot, 0, 0, 0]
tspan = (0.0,30.0)

trajectories_num = 200
D1_D2_vec = [1.0, 3.0, 10.0, 30.0]
D3_vec = [0.0, 30.0]

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
        ic[17] = D3_vec[k]    
        jump_prob  = generate_jump_problem(ic, tspan, p)
        ens_prob = EnsembleProblem(prob = jump_prob)
        sim = solve(ens_prob, SSAStepper(), EnsembleThreads(), trajectories=trajectories_num)
        for j = 1:lastindex(sim.u)
            A1_vec[j] = sim[j].u[end][9]
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
        
        hist1 = histogram(A1_vec, bins =10,  normalize = :pdf, legend = :none, fillcolor = :gray, lc = :gray ,fillalpha= 0.5, grid=false)
        plot!(a1_kde.x, a1_kde.density, lc = :black, lw = 2)
        xlims!((0, extrema(A1_vec)[2]))
        xlabel!("A1 [protein counts]")
        ylabel!("Probability Density")
        savefig(hist1, string("jl_output/2023.05.10_toggle/A1_histogram_base_params_D1_", D1_D2_vec[i], 
                                "D3tot_", D3_vec[k] , "_short.svg"))
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
savefig(p1, "jl_output/2023.05.10_toggle/D3vsD1_A1_counts_short.svg")

p2 = plot(D3_vec, frac_A1_high, line_z = D1_D2_mat, label = "", zcolor = D1_D2_vec, grid = false, c= cgrad(:viridis) , colorbar_title = "D1, D2")
xlabel!("D3 total")
ylabel!("Fraction A1-high cells")
savefig(p2, "jl_output/2023.05.10_toggle/D3vsD1_A1_frac_short.svg")

p3 =plot(D3_vec, cv_A1_high, line_z = D1_D2_mat, label = "", zcolor = D1_D2_vec, grid = false, c= cgrad(:viridis), colorbar_title = "D1, D2" )
xlabel!("D3 total")
ylabel!("Coefficient of Variation\nof A1-high cells")
savefig(p3, "jl_output/2023.05.10_toggle/D3vsD1_A1_cv_short.svg")


##
p = DEFAULT_SIM_PARAMS
ic = [DEFAULT_CELL_PARAMS.D1tot, DEFAULT_CELL_PARAMS.D2tot, 0, 0, 0, 0,
    0, 0, 70, 0, 0, 0,
    0, 0, 0, 0, DEFAULT_CELL_PARAMS.D3tot, 0, 0, 0]
tspan = (0.0,80.0)

trajectories_num = 100
kappa_vec = [5.0, 2.0, 1.0]
gamma_vec = [1.0, 0.4, 0.2]
D3_vec = [0.0, 30.0]

for i = 1:lastindex(kappa_vec)
    p.Toggle_Params.κ = kappa_vec[i]
    p.Toggle_Params.γ = gamma_vec[i]
    for k = 1:lastindex(D3_vec)
        ic[17] = D3_vec[k]
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

savefig(current(), "jl_output/2023.05.10_toggle/D3vskappa_switching.svg")

p = DEFAULT_SIM_PARAMS
ic = [DEFAULT_CELL_PARAMS.D1tot, DEFAULT_CELL_PARAMS.D2tot, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, DEFAULT_CELL_PARAMS.D3tot, 0, 0, 0]
tspan = (0.0,30.0)
trajectories_num = 200
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
        ic[17] = D3_vec[k]    
        jump_prob  = generate_jump_problem(ic, tspan, p)
        ens_prob = EnsembleProblem(prob = jump_prob)
        sim = solve(ens_prob, SSAStepper(), EnsembleThreads(), trajectories=trajectories_num)
        for j = 1:lastindex(sim.u)
            A1_vec[j] = sim[j].u[end][9]
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
        
        hist1 = histogram(A1_vec, bins =10,  normalize = :pdf, legend = :none, fillcolor = :gray, lc = :gray ,fillalpha= 0.5, grid=false)
        plot!(a1_kde.x, a1_kde.density, lc = :black, lw = 2)
        xlims!((0, extrema(A1_vec)[2]))
        xlabel!("A1 [protein counts]")
        ylabel!("Probability Density")
        savefig(hist1, string("jl_output/2023.05.10_toggle/A1_histogram_base_params_kappa_", kappa_vec[i], 
                                "D3tot_", D3_vec[k] , "_short.svg"))
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
savefig(p1, "jl_output/2023.05.10_toggle/D3vskappa_A1_counts_short.svg")

p2 = plot(D3_vec, frac_A1_high, line_z = kappa_mat, label = "", zcolor = 2.0 ./ kappa_vec, grid = false, c= cgrad(:viridis) , colorbar_title = "α/κ")
xlabel!("D3 total")
ylabel!("Fraction A1-high cells")
savefig(p2, "jl_output/2023.05.10_toggle/D3vskappa_A1_frac_short.svg")

p3 = plot(D3_vec, cv_A1_high, line_z = kappa_mat, label = "", zcolor = 2.0 ./ kappa_vec, grid = false, c= cgrad(:viridis), colorbar_title = "α/κ" )
xlabel!("D3 total")
ylabel!("Coefficient of Variation\nof A1-high cells")
savefig(p3, "jl_output/2023.05.10_toggle/D3vskappa_A1_cv_short.svg")


end