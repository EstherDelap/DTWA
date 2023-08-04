function spin_sqeezing_param_3D(collective_spin, average_spin, number_spins)
	#Find the spin sqeezing parameter as defined in equation 31 in Peter Rabls paper PhsRevA.105.013716. Returns an array of the time 
	#evolution of the square of the spin squeezing parameter, if the spins evolve in time as collective_spin
	ξ_squared_time =Vector{Float64}()
	N = size(average_spin)[1] #total number of timesteps
	for t in 1:N
		mean_spin_vector = normalize(average_spin[t]) #average across all trajectories, so this is a vector [sx, sy, sz]
		collective_spin_all_traj_at_t = getindex.(collective_spin[:],t) #[[sx(t),sy(t),sz(t)]_traj1,[sx(t),sy(t),sz(t)]_traj2, ...]
		θ = acos(dot([0,0,1],mean_spin_vector))
		k= normalize(cross([0,0,1],mean_spin_vector))
		x=normalize((cos(θ)*[1,0,0]) + (sin(θ)*cross(k,[1,0,0])) + ((1-cos(θ))*k*dot(k,[1,0,0])))
		y=normalize((cos(θ)*[0,1,0]) + (sin(θ)*cross(k,[0,1,0])) + ((1-cos(θ))*k*dot(k,[0,1,0])))
		res = optim_perp(x,y,collective_spin_all_traj_at_t)
		ξ_squared = (minimum(res)^2)*(number_spins/dot(mean_spin_vector,mean_spin_vector))
		push!(ξ_squared_time, ξ_squared)
	end
	return ξ_squared_time
end

function optim_perp(x,y,S)
	#needed to find the minimum in the spin squeezing parameter
	f = ϕ -> s_perp(ϕ, x, y, S)
	return optimize(f, 0,2*π)
end

function s_perp(ϕ,x,y,S)
	#needed to find the minimum in the spin squeezing parameter
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