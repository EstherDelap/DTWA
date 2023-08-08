var documenterSearchIndex = {"docs":
[{"location":"#DTWA.jl","page":"DTWA.jl","title":"DTWA.jl","text":"","category":"section"},{"location":"","page":"DTWA.jl","title":"DTWA.jl","text":"Documentation for DTWA.jl","category":"page"},{"location":"","page":"DTWA.jl","title":"DTWA.jl","text":"DTWA.repeated_euler\nDTWA.euler_3D\nDTWA.spin_array_3D\nDTWA.intial_Si\nDTWA.J_XYZ\nDTWA.Jx_Ising","category":"page"},{"location":"#DTWA.repeated_euler","page":"DTWA.jl","title":"DTWA.repeated_euler","text":"repeated_euler(dim::Vector{Int64}, N::Int64,number_repeats::Int64, Γ_deph::Float64, Γ_decay::Float64,Ω::Float64, α::AbstractVector, method::String, axis::Int64, dir::Int64)\n\nRepeatedly computes the time evolution of the spin system of dimensions dim for a certain time interval.\n\nReturns two arrays: the first is the time evolution of the collective spin, defined as left sum_i^n s_i^x2 sum_i^n s_i^y2 sum_i^n s_i^z2right for each trajectory, in the form of a vector of length number_trajectories where each value is a vector of length N decribing the time evolution of the collective spin. The second array returned is the time evolution of the collective spin averaged across all trajectories,  in the form of a vector of length N where each element is the averaged [sx, sy, sz] spin vector. \n\nArguments\n\ndim::Vector{Int64}: the dimensions of the spin system. Assumes a 3d Cartesian spin system, but can be made two dimensional by using 1\n`N::Int64: the number of time-steps\nnumber_repeats::Int64: the number of repeated trajectories\nΓ_deph::Float64: the rate of local dephasing\nΓ_decay::Float64: the rate of local decay\nΩ::Float64: the external field strength\nα::AbstractVector: A vector of length 3, for the Ising model interactions this is three of the same integers, [α,α,α] representing the interaction decay between negihbouring spins,\n\nwhile for the XYZ model this is three possibly different floats, [jx,jy,jz] indicating the strength of interaction in the x, y and z directions\n\nmethod::String: either \"Ising\" or \"XYZ\" indicating the model of interactions\naxis::Int64: the axis along which the spins are intially aligned, either 1, 2 or 3 for the x, y or z axis respectively\ndir::Int64: the direction along which the spins are intially aligned\n\n\n\n\n\n","category":"function"},{"location":"#DTWA.euler_3D","page":"DTWA.jl","title":"DTWA.euler_3D","text":"euler_3D(N::Int64, time_interval::Tuple{Int64}, S_0::Array{Vector{Float64}}, Γ_deph::Float64, Γ_decay::Float64, Ω::Float64, Jx::Interactions, Jy::Interactions, Jz::Interactions)\n\nComputes the time evolution of the initial spin array S_0 over the specified time_interval with N time steps.\n\nReturns the time evolution of the collective spin, defined asleft sum_i^n s_i^x2 sum_i^n s_i^y2 sum_i^n s_i^z2right.  This is in the form of a vector of length N where each element is a vector of length 3.\n\nArguments\n\nN::Int64: the number of time-steps\ntime_interval::Tuple{Int64}: the time interval in the form (t_0,t_final)\nS_0::Array{Vector{Float64}}: A three dimensional initial spin array, each element being a vector of length 3 indicating the [sx, sy, sz] spin values for that position\nΓ_deph::Float64: the rate of local dephasing\nΓ_decay::Float64: the rate of local Γ_decay\nΩ::Float64: the external field strength\nJx::Interactions: an Interactions matrix indeication the interactions between spins in the x direction. For the Ising model, this is zero.\nJy::Interactions: an Interactions matrix indeication the interactions between spins in the y direction. For the Ising model, this is zero.\nJz::Interactions: an Interactions matrix indeication the interactions between spins in the z direction\n\n\n\n\n\n","category":"function"},{"location":"#DTWA.spin_array_3D","page":"DTWA.jl","title":"DTWA.spin_array_3D","text":"spinarray3D(dim::Vector{Int64}, axis::Int64, dir::Int64)\n\nReturns a 3D array of spins of dimensions dim, where the spins are aligned along axis in the direction dir. Each element is a vector of length 3 indicating the x, y and z spins [sx, sy, sz]. \n\nArguments\n\ndim::Vector{Int64}: the dimensions of the spin array  \naxis::Int64: the axis along which the spins are aligned, either 1, 2 or 3 to indicate the x, y or z axis respectively\ndir::Int64: the direction in which the spins are aligned along axis, either 1 or -1\n\n\n\n\n\n","category":"function"},{"location":"#DTWA.intial_Si","page":"DTWA.jl","title":"DTWA.intial_Si","text":"intial_Si(i::Int64,dir::Int64)\n\nReturns an initial spin for a single position in a spin array, as a vector of length 3, [sx, sy, sz]. The spin is initially  aligned along the axis indicated by i (equal 1, 2 or 3 to indeicate x, y or z respectively) in the direction dir, which is either 1 or -1.  \n\n\n\n\n\n","category":"function"},{"location":"#DTWA.J_XYZ","page":"DTWA.jl","title":"DTWA.J_XYZ","text":"J_XYZ(dims::Vector{Int64},j::Float64)\n\nReturns the interaction matrix for the XYZ model with periodic boundary conditions.\n\nFor the XYZ model, the interactions matrix for a specific axis is a nearest neighbour interaction of strength j with periodic boundary conditions. Returns a six dimensional matrix h that is an instance of the struct Interactions, such that h[r1,r2] = j if r1 and r2 are cartesian indeces that are neighbours, or if they are cartesian indeces at opposite ends of the matrix (to create all periodic boundary conditions), and h[r1,r2] = 0.0  for all other cartesian indices.\n\n\n\n\n\n","category":"function"},{"location":"#DTWA.Jx_Ising","page":"DTWA.jl","title":"DTWA.Jx_Ising","text":"Jx_Ising(dim::Vector{Int64})\n\nReturns an instance h of the struct Interactions where h[r1,r2]=0.0 for all Cartesian Indeces r1, r2. \n\nIn the transverse Ising model, spins dont interact in the x or y directions, so all interactions are zero.\n\n\n\n\n\n","category":"function"}]
}
