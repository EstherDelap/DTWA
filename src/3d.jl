struct Interactions
	data::Array{Float64,6}
end

function Interactions(dim::NTuple{3,Int})
	h = Interactions(Array{Float64}(undef,dim[1],dim[2],dim[3],dim[1],dim[2],dim[3]))
	return h
end

function Base.getindex(interaction_matrix::Interactions, r1::Vector{Int64}, r2::Vector{Int64})
	return interaction_matrix.data[r1[1], r1[2], r1[3], r2[1], r2[2], r2[3]]
end

function Base.setindex!(interaction_matrix::Interactions, d::Float64, r1, r2)
	interaction_matrix.data[r1[1], r1[2], r1[3], r2[1], r2[2], r2[3]] = d
end

LinearAlgebra.norm(i::CartesianIndex) = sqrt(i[1]^2 + i[2]^2 + i[3]^2)

function J_matrix(dim, α)
	h = Interactions(Array{Float64}(undef,dim[1],dim[2],dim[3],dim[1],dim[2],dim[3]))
	n = dim[1]*dim[2]*dim[3]
    for r1 in CartesianIndices(h.data)
		for r2 in CartesianIndices(h.data)
            if r1 == r2
                h[r1,r2] = 0.0
            else
			    h[r1,r2] = 1/(n*(norm(r1-r2)^α))
            end
		end
	end
	return h
end

function spin_array_3D(dim, axis, dir)
	#dim is a tuple of length 3, 'axis' is the axis (x y or z) along which the spin is aligned, and dir is the direction (1 or -1) of this alignement. Produced a 3 dimensional matrix of vectors, so each point is a vector [Sx, Sy, Sz]
	spins = Array{Vector{Float64}}(undef, dim)
	for i in 1:dim[1]
		for j in 1:dim[2]
			for k in 1:dim[3]
				spins[i,j,k] = intial_Si(axis,dir)
			end
		end
	end
	return spins
end

function euler_3D(N, time_interval, S_0, Γ_deph, Γ_decay, Ω, α, J)
	#J is a 6 dimensional matrix of interaction strengths, so J(index1,index2) gives the interaction strength between spin at index 1 and spin at index 2
	#N is the numer of time steps, time_interval is a tuple of the form (t_0, t_final), S_0 is the initial spin vector of the form S_0 = []
	#N is the numer of time steps, time_interval is a tuple of the form (t_0, t_final), S_0 is the initial spin vector of the form S_0 = [sx_1 sy_1 sz_1; sx_2 sy_2 sz_2; ...] 
	collective_spin = Vector{Vector{Float64}}(undef,N) #collective_spin[t]=[sx(t),sy(t),sz(t)]

	dt = time_interval[2]/N
	S = deepcopy(S_0)
	dim = size(S)
	number_atoms = dim[1]*dim[2]*dim[3]
	S_new = Array{Vector{Float64}}(undef, dim)
	
	sx = 0
	sy = 0
	sz = 0
	for spin in S_0
		sx+=spin[1]
		sy+=spin[2]
		sz+=spin[3]
	end
	collective_spin[1] = [sx/2, sy/2, sz/2] #collective spin at time t=0
	for t in 2:N #we have already done the first one
		sx = 0
		sy = 0
		sz = 0
		for i in 1:dim[1]
			for j in 1:dim[2]
				for k in 1:dim[3]
					dW = randn()*dt #different for each spin
					sum_x = 0
					sum_y = 0
					for l in 1:dim[1]
						for m in 1:dim[2]
							for n in 1:dim[3]
								J_ij = J[[i,j,k],[l,m,n]]
								sum_x+=2*J_ij*S[i,j,k][2]*S[l,m,n][3]
								sum_y+=2*J_ij*S[i,j,k][1]*S[l,m,n][3]
							end
						end
					end
					x= S[i,j,k][1]  - (sum_x + Γ_deph *S[i,j,k][1]+ Γ_decay/2 * S[i,j,k][1])*dt - (sqrt(2*Γ_deph)+sqrt(Γ_decay))*S[i,j,k][2]*dW
					y= S[i,j,k][2] + (sum_y -(Ω*S[i,j,k][3]) - (Γ_deph+Γ_decay/2)*S[i,j,k][2])*dt +(sqrt(2*Γ_deph)+sqrt(Γ_decay))*S[i,j,k][1]*dW
					z = S[i,j,k][3] + (Ω*S[i,j,k][2] - Γ_decay*(S[i,j,k][3]+1))*dt + (sqrt(Γ_decay)*(S[i,j,k][3]+1))*dW
					S_new[i,j,k] = [x,y,z]

					sx += x
					sy += y
					sz += z
				end
			end
		end
		collective_spin[t] = [sx/2, sy/2, sz/2]
		S = deepcopy(S_new)
	end
	return collective_spin
end

function spin_sqeezing_param_3D(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, α)
	J = J_matrix(dim, α)
    number_spins = dim[1]*dim[2]*dim[3]
	J_bar = sum(J.data)/number_spins
    time_interval = (0.0, 10.0/J_bar)
	collective_spin_all_traj = Vector{Vector{Vector{Float64}}}(undef, number_repeats)
	initial_traj = euler_3D(N, time_interval, spin_array_3D(dim, 1, 1), Γ_deph, Γ_decay, Ω, α,J)[2]
	collective_spin_all_traj[1] = initial_traj
	average_collective_spin = initial_traj
	for i in 2:number_repeats
		traj = euler_3D(N, time_interval, spin_array_3D(dim, 1, 1), Γ_deph, Γ_decay, Ω,α,J)[2]
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
		ξ_squared = (minimum(res)^2)*(number_spins/dot(mean_spin_vector,mean_spin_vector))
		push!(ξ_squared_time, ξ_squared)
	end
	return ξ_squared_time
end

function z_magnetisation_3D(dim, N, time_interval, number_repeats, Γ_deph, Γ_decay, Ω, α)
	J = J_matrix(dim, α)
	S_0 =  spin_array_3D(dim, 3, -1)
	Sz_cummulative =  euler_3D(N, time_interval, S_0, Γ_deph, Γ_decay, Ω, α, J)[1]
	for i in 2:number_repeats #we have already done the first one
		Sz_cummulative += euler_3D(N, time_interval, spin_array_3D(dim, 3, -1), Γ_deph, Γ_decay, Ω, α, J)[1]
	end
	return Sz_cummulative/number_repeats
end