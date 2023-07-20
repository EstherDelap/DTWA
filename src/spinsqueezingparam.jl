function spin_sqeezing_param(n, N, time_interval, number_repeats, Γ_deph, Γ_decay, Ω)
	collective_spin_all_traj = Vector{Vector{Vector{Float64}}}(undef, number_repeats)
	initial_traj =  euler(N, time_interval, initial_spin_array(n,1,1), Γ_deph, Γ_decay, Ω)[2]
	collective_spin_all_traj[1] = initial_traj
	average_collective_spin = initial_traj
	
	for i in 2:number_repeats
		traj = euler(N, time_interval, initial_spin_array(n,1,1), Γ_deph, Γ_decay, Ω)[2]
		collective_spin_all_traj[i] = traj
		average_collective_spin += traj
	end
	
	ξ_squared_time =Vector{Float64}()
	for t in 1:N
		mean_spin_vector = normalize(average_collective_spin[t]/number_repeats)
		collective_spin_all_traj_at_t = getindex.(collective_spin_all_traj[:],t)
		θ = acos(dot([0,0,1],mean_spin_vector))
		k= normalize(cross([0,0,1],mean_spin_vector))
		x=normalize((cos(θ)*[1,0,0]) + (sin(θ)*cross(k,[1,0,0])) + ((1-cos(θ))*k*dot(k,[1,0,0])))
		y=normalize((cos(θ)*[0,1,0]) + (sin(θ)*cross(k,[0,1,0])) + ((1-cos(θ))*k*dot(k,[0,1,0])))
		res = optim_perp(x,y,collective_spin_all_traj_at_t)
		ξ_squared = (minimum(res)^2)*(n/dot(mean_spin_vector,mean_spin_vector))
		push!(ξ_squared_time, ξ_squared)
	end
	return ξ_squared_time
end

function optim_perp(x,y,S)
	f = ϕ -> s_perp(ϕ, x, y, S)
	return optimize(f, 0,2*π)
end

function s_perp(ϕ,x,y,S)

	n = x*cos(ϕ)+y*sin(ϕ)
	
	s_perp_ave = 0
	s_perp_squared = 0
	for trajectory in S
		traj = normalize(trajectory)
		s_perp_ave+=dot(traj,n)
		s_perp_squared+=(dot(traj,n))^2
	end
	return sqrt((s_perp_squared/length(S))-((s_perp_ave/length(S))^2))
end