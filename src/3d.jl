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


"""
   Jz_Ising(dim::Vector{Int64},α::Int64)

Produces the interaction matrix in the z direction for the transverse Ising model.

Returns an instance h of the struct Interactions, where h[r1,r2]=0.0 for all Cartesian Indeces that are the same, and 
for all other Cartesian Indeces. d(r1,r2) is the Euclidean distance between the cartesian indeces r1 and r2, and N is the total number of spins. 
"""
function Jz_Ising(dims,α)
	#the transverse Ising Hamiltonian. the diagonal is zeroes. the interactions between two particles depends on their distance
	#J_ij = 1/N * 1/d(spin i, spin j)^α
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

"""
   Jx_Ising(dim::Vector{Int64})

Returns an instance h of the struct Interactions where h[r1,r2]=0.0 for all Cartesian Indeces r1, r2. 
	
In the transverse Ising model, spins dont interact in the x or y directions, so all interactions are zero.
"""
function Jx_Ising(dim)
	#these are for the x and y directions for the ising model, so that the same equation can be used for the xyz model, dim is a vector
	h = Interactions(zeros(Float64,dim[1],dim[2],dim[3],dim[1],dim[2],dim[3]))
	return h
end

"""
    J_XYZ(dims::Vector{Int64},j::Float64)

Returns the interaction matrix for the XYZ model with periodic boundary conditions.

For the XYZ model, the interactions matrix for a specific axis is a nearest neighbour interaction of strength j with periodic boundary conditions.
Returns a six dimensional matrix h that is an instance of the struct Interactions, such that h[r1,r2] = j if r1 and r2 are cartesian indeces that
are neighbours, or if they are cartesian indeces at opposite ends of the matrix (to create all periodic boundary conditions), and h[r1,r2] = 0.0 
for all other cartesian indices.

"""
function J_XYZ(dims,j)
	#j needs to be a floating number. Diagonals are zero. Nearest neighbour interactions are equal to j, so those next to the diagonal.
	#we use periodic boundary conditions so the edges are also equal to j. 
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
				r2= [r1[1],r1[2],r1[3]]
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

"""
    intial_Si(i::Int64,dir::Int64) 

Returns an initial spin for a single position in a spin array, as a vector of length 3, [sx, sy, sz]. The spin is initially aligned along the axis
indicated by i (1->x axis, 2->y axis, 3->z axis) in the direction dir, which is either 1 or -1.  
"""
function intial_Si(i,dir)
	#initial spin at sight i. S[i] is set to -1. Sx and Sy are chosen randomly from [-1,1] Returns [Sx,Sy,Sz,S0]
	S = [rand([-1,1]), rand([-1,1]), rand([-1,1])]
	S[i] = 1.0*dir
	return Vector{Float64}(S)
end

"""
   spin_array_3D(dim::Vector{Int64}, axis::Int64, dir::Int64)

Returns a 3D array of spins of dimensions dim, where the spins are aligned aling axis in the direction dir. Each element is a vector of length 3 
indicating the x, y and z spins [sx, sy, sz]. 

# Arguments
- dim::Vector{Int64}: the dimensions of the spin array  
- axis::Int64: the axis along which the spins are aligned, either 1, 2 or 3 to indicate the x, y or z axis respectively
- dir::Int64: the direction in which the spins are aligned along axis, either 1 or -1
"""
function spin_array_3D(dim::Vector{Int64}, axis, dir)
	#dim is a tuple of length 3, 'axis' is the axis (x y or z) along which the spin is aligned, and dir is the direction (1 or -1) of this alignement.
	#Produces a 3 dimensional matrix of vectors, so each point is a vector [Sx, Sy, Sz]
	spins = Array{Vector{Float64}}(undef, dim[1], dim[2], dim[3])
	for i in 1:dim[1]
		for j in 1:dim[2]
			for k in 1:dim[3]
				spins[i,j,k] = intial_Si(axis,dir)
			end
		end
	end
	return spins
end
"""
    euler_3D(N::Int64, time_interval::Tuple{Int64}, S_0::Array{Vector{Float64}}, Γ_deph::Float64, Γ_decay::Float64, Ω::Float64, Jx::Interactions, Jy::Interactions, Jz::Interactions)

Computes the time evolution of the initial spin array S_0 over the specified time_interval with N timesteps.

Returns the time evolution of the collective spin, defined as. This is in the form of 
a vector of length N where each element is a vecotr of length 3.

# Arguments
- N::Int64: the number of time-steps
- time_interval::Tuple{Int64}: the time interval in the form (t_0,t_final)
- S_0::Array{Vector{Float64}}: A three dimensional initial spin array, each element being a vector of length 3 indicating the [sx, sy, sz] spin values for that position
- Γ_deph::Float64: the rate of local dephasing
- Γ_decay::Float64: the rate of local Γ_decay
- Ω::Float64: the external field strength
- Jx::Interactions: an Interactions matrix indeication the interactions between spins in the x direction. For the Ising model, this is zero.
- Jy::Interactions: an Interactions matrix indeication the interactions between spins in the y direction. For the Ising model, this is zero.
- Jz::Interactions: an Interactions matrix indeication the interactions between spins in the z direction
"""
function euler_3D(N, time_interval, S_0, Γ_deph, Γ_decay, Ω, Jx, Jy, Jz)
	#returns the time evolution of the spins with the initial spin array S_0. Works with both Ising model and XYZ model. 
	#N is the number of time steps, Jx, Jy and Jz are the interactioon matrices (Jx=Jy=0 for the Ising model). 
	#return collective_spin, a vector of length N such that  collective_spin[t]=[ave(sx[t]),ave(sy[t]),ave(sz[t])] 
	#is the average value of all three spins at time t, ie averaged over all the spins. 
	collective_spin = Vector{Vector{Float64}}(undef,N) #collective_spin[t]=[sx(t),sy(t),sz(t)]

	dt = time_interval[2]/N
	S = deepcopy(S_0)
	dim = size(S)
	number_atoms = dim[1]*dim[2]*dim[3]
	S_new = Array{Vector{Float64}}(undef, dim[1],dim[2],dim[3])
	
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
								
								Jx_ij = Jx[[i,j,k],[l,m,n]]
								Jy_ij = Jy[[i,j,k],[l,m,n]]
								Jz_ij = Jz[[i,j,k],[l,m,n]]

								sum_x+=2*Jy_ij*S[i,j,k][3]*S[l,m,n][2]-2*Jz_ij*S[i,j,k][2]*S[l,m,n][3]
								
								sum_y+=2*Jz_ij*S[i,j,k][1]*S[l,m,n][3]-2*Jx_ij*S[i,j,k][3]*S[l,m,n][1]
								
								sum_z+=2*Jx_ij*S[i,j,k][2]*S[l,m,n][1]-2*Jy_ij*S[i,j,k][1]*S[l,m,n][2]
							end
						end
					end
					x= S[i,j,k][1]  + (sum_x - (Γ_deph + Γ_decay/2)* S[i,j,k][1])*dt - (sqrt(2*Γ_deph)+sqrt(Γ_decay))*S[i,j,k][2]*dW
					
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
	return collective_spin #vector of length N, each item is a vector of length 3
end

"""
    repeated_euler(dim::Vector{Int64}, N::Int64,number_repeats::Int64, Γ_deph::Float64, Γ_decay::Float64,Ω::Float64, α::AbstractVector, method::String, axis::Int64, dir::Int64)

Repeatedly computes the time evolution of the spin system of dimensions dim for a certain time interval.

Returns two arrays: the first is the time evolution of the collective spin, defined as 
for each trajectory, in the form of a vecotr of length number_trajectories where each value is a vector of length N decribing the time
evolution of the collective spin. The second array returned is the time evolution of the collective spin averaged across all trajectories, in the form of a 
vector of length N where each element is the averaged [sx, sy, sz] spin vector. 
 
# Arguments
- dim::Vector{Int64}: the dimensions of the spin system. Assumes a 3d Cartesian spin system, but can be made two dimensional by using 1
- N::Int64: the number of time-steps
- number_repeats::Int64: the number of repeated trajectories
- Γ_deph::Float64: the rate of local dephasing
- Γ_decay::Float64: the rate of local Γ_decay
- Ω::Float64: the external field strength
- α::AbstractVector: A vector of length 3, for the Ising model interactions this is three of the same integers, [α,α,α] representing the interaction decay between negihbouring spins,
 while for the xyz model this is three possibly different floats, [jx,jy,jz] indicating the strength of interaction in the x, y and z directions
- method::String: either "Ising" or "XYZ" indicating the model of interactions
- axis::Int64: the axis along which the spins are intially aligned, either 1, 2 or 3 for the x, y or z axis respectively
- dir::Int64: the direction along which the spins are intially aligned
"""
function repeated_euler(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, α, method, axis, dir)
	#dim is a vector of integers, N is the number of timesteps, α is a vector of length 3; in the Ising model, it is three values,
	#in the XYZ model, α is the jx, jy, jz values. method id wither "Ising" of "XYZ". axis is the axis along with the spins are
	#initially aligned and dir is the direction along this axis. Performs euler_3D for number_repeats. Return the 
	#collective_spin_all_traj, a vector of length number_repeats uch that collective_spin_all_traj[n] is the nth trajectory 
	#of the time evolution of the spins, so each element in collective_spin_all_traj is a vector of length N. Thus to get all the
	#average spin values for a specific time t across all the all the trajectories we getindex.(collective_spin_all_traj[:],t)
	#also return the average, a vector of length N, so that average[t]=[<sx[t]>,<sy[t]>,<sz[t]>] (<sx[t]> is average over all traj)
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
    #time_interval = (0.0, 10.0/J_bar)"
	time_interval = (0.0,3000.0)
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