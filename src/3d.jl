struct Interactions
	data::Array{Float64,6}
end

Interactions(dim::NTuple{3,Int}) = Interactions(dim...)

function Interactions(d1, d2, d3)
	return Interactions(Array{Float64}(undef,d1,d2,d3,d1,d2,d3))
end

function Base.getindex(interaction_matrix::Interactions, r1::Vector{Int64}, r2::Vector{Int64})
	return interaction_matrix.data[r1[1], r1[2], r1[3], r2[1], r2[2], r2[3]]
end

function Base.setindex!(interaction_matrix::Interactions, d::Float64, r1, r2)
	interaction_matrix.data[r1[1], r1[2], r1[3], r2[1], r2[2], r2[3]] = d
end

dims(ints::Interactions) = size(ints.data)[1:3]

Base.CartesianIndices(ints::Interactions) = CartesianIndices(dims(ints))

LinearAlgebra.norm(i::CartesianIndex) = sqrt(i[1]^2 + i[2]^2 + i[3]^2)

function Jz_Ising(dims,α)
	n = *(dims...)
	h = Interactions(dims[1],dims[2],dims[3]) # undef array
	for r1 in CartesianIndices(h.data)
		for r2 in CartesianIndices(h.data)
			if norm(r1-r2) == 0
				h[r1,r2] = 0.0
			else
				h[r1,r2] = 1/(n*(norm(r1-r2)^α))
			end
		end
	end
	return h
end

function Jx_Ising(dim)
	h = Interactions(zeros(Float64,dim[1],dim[2],dim[3],dim[1],dim[2],dim[3]))
	return h
end

function J_XYZ(dims,j)
	h = Interactions(dims[1], dims[2], dims[3])
	for r1 in CartesianIndices(h)
		for r2 in CartesianIndices(h)
			if norm(r1-r2) == 1.0
				h[r1,r2] = j
			else
				h[r1,r2] = 0.0
			end
		end
		#now the part for periodic boundary conditions...
		for (index,value) in enumerate(Tuple(r1))
			if value==1
				r2 = [r1[1],r1[2],r1[3]]
				r2[index] = dims[index]
				h[r1,CartesianIndex(Tuple(r2))] = j
			elseif value == dims[index] #ie if its at the edge
				r2 = [r1[1],r1[2],r1[3]]
				r2[index] = 1
				h[r1,CartesianIndex(Tuple(r2))] = j
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

function euler_3D(N, time_interval, S_0, Γ_deph, Γ_decay, Ω, Jx, Jy, Jz)
	
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
					sum_z = 0
					for l in 1:dim[1]
						for m in 1:dim[2]
							for n in 1:dim[3]
								
								Jy_ij = Jx[[i,j,k],[l,m,n]]
								Jz_ij = Jy[[i,j,k],[l,m,n]]
								Jx_ij = Jz[[i,j,k],[l,m,n]]

								sum_x+=2*Jy_ij*S[i,j,k][3]*S[l,m,n][2]-2*Jx_ij*S[i,j,k][2]*S[l,m,n][3]
								
								sum_y+=2*Jz_ij*S[i,j,k][1]*S[l,m,n][3]-2*Jx_ij*S[i,j,k][3]*S[l,m,n][1]
								
								sum_z+=2*Jx_ij*S[i,j,k][2]*S[l,m,n][1]-2*Jy_ij*S[i,j,k][1]*S[l,m,n][2]
							end
						end
					end
					x= S[i,j,k][1]  - (sum_x + Γ_deph *S[i,j,k][1]+ Γ_decay/2 * S[i,j,k][1])*dt - (sqrt(2*Γ_deph)+sqrt(Γ_decay))*S[i,j,k][2]*dW
					
					y= S[i,j,k][2] + (sum_y -(Ω*S[i,j,k][3]) - (Γ_deph+Γ_decay/2)*S[i,j,k][2])*dt +(sqrt(2*Γ_deph)+sqrt(Γ_decay))*S[i,j,k][1]*dW
					
					z = S[i,j,k][3] + (sum_z + Ω*S[i,j,k][2] - Γ_decay*(S[i,j,k][3]+1))*dt + (sqrt(Γ_decay)*(S[i,j,k][3]+1))*dW
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



function repeated_euler(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, α, method, axis, dir)
	#α is the α is the dissipative model and α = [jx,jy,jz] in the XYZ model
	if method == "Ising"
		Jx = Jx_Ising(dim)
		Jy = Jx_Ising(dim)
		Jz = Jz_Ising(dim,α[1])
	else
		Jx =J_XYZ(dim,α[1])
		Jy =J_XYZ(dim,α[2])
		Jz =J_XYZ(dim,α[3])
	end

    number_spins = dim[1]*dim[2]*dim[3]
	J_bar = sum(Jz.data)/number_spins
    time_interval = (0.0, 10.0/J_bar)
	collective_spin_all_traj = Vector{Vector{Vector{Float64}}}(undef, number_repeats)
	initial_traj = euler_3D(N, time_interval, spin_array_3D(dim, axis, dir), Γ_deph, Γ_decay, Ω, Jx, Jy, Jz)
	collective_spin_all_traj[1] = initial_traj
	average = initial_traj
	for i in 2:number_repeats
		traj = euler_3D(N, time_interval, spin_array_3D(dim, axis, dir), Γ_deph, Γ_decay, Ω, Jx, Jy, Jz)
		collective_spin_all_traj[i] = traj
		average += traj
	end
	return collective_spin_all_traj, average/number_repeats
end

function spin_sqeezing_param_3D(collective_spin, average_spin, number_spins)
	ξ_squared_time =Vector{Float64}()
	N = size(average_spin)[1]
	for t in 1:N
		mean_spin_vector = normalize(average_spin[t])
		collective_spin_all_traj_at_t = getindex.(collective_spin[:],t)
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

function magnetisation_3D(average, axis)
	#axis refers to 1 for x, 2 for y or 3 for z, depending along where you want to measure magnetisation
	return getindex.(average[:],axis)
end