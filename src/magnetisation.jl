function euler(N, time_interval, S_0, Γ_deph, Γ_decay, Ω)
	#N is the numer of time steps, time_interval is a tuple of the form (t_0, t_final), S_0 is the initial spin vector of the form S_0 = []
	#N is the numer of time steps, time_interval is a tuple of the form (t_0, t_final), S_0 is the initial spin vector of the form S_0 = [sx_1 sy_1 sz_1; sx_2 sy_2 sz_2; ...] 
	collective_spin = Vector{Vector{Float64}}(undef,N)
	J=1
	dt = time_interval[2]/N
	S = deepcopy(S_0)
	n = size(S)[1]
	S_new = Array{Float64}(undef,n,3)
	
	sx = 0
	sy = 0
	sz = 0
	for spin in eachrow(S_0)
		sx+=spin[1]
		sy+=spin[2]
		sz+=spin[3]
	end
	collective_spin[1] = [sx/2, sy/2, sz/2]
	Sz = [sz]
	for t in 2:N #we have already done the first one
		sx = 0
		sy = 0
		sz = 0
		for (i,row) in enumerate(eachrow(S))
			dW = randn()*dt
			sum_x = 0
			sum_y = 0
			for (j,neighbours) in enumerate(eachrow(S))
				if i!=j
					sum_x+=2*J/n*row[2]*neighbours[3]
					sum_y+=2*J/n*row[1]*neighbours[3]
				end
			end
			S_new[i,1] = row[1] - (sum_x + Γ_deph *row[1]+ Γ_decay/2 * row[1])*dt -(sqrt(2*Γ_deph)+sqrt(Γ_decay))*row[2]*dW
			S_new[i,2] = row[2] + (sum_y -(Ω*row[3]) - (Γ_deph+Γ_decay/2)*row[2])*dt + (sqrt(2*Γ_deph)+sqrt(Γ_decay))*row[1]*dW
			S_new[i,3] = row[3] + (Ω*row[2] - Γ_decay*(row[3]+1))*dt + (sqrt(Γ_decay)*(row[3]+1))*dW
			
			sx+=S_new[i,1]
			sy+=S_new[i,2]
			sz+=S_new[i,3]
		end
		
		collective_spin[t] = [sx/2, sy/2, sz/2]
		copy!(S, S_new)
		push!(Sz, sz)
	end
	return Sz, collective_spin
end

function intial_Si(i,dir)
	#initial spin at sight i. S[i] is set to -1. Sx and Sy are chosen randomly from [-1,1] Returns [Sx,Sy,Sz,S0]
	Sy = rand([-1,1])
	Sz = rand([-1,1])
	S = [rand([-1,1]), rand([-1,1]), rand([-1,1])]
	S[i] = 1.0*dir
	return Vector{Float64}(S)
end

function initial_spin_array(n,i,dir)
	#n is the number of spins, produces a n*3 matrix of the form S = [sx_1 sy_1 sz_1; sx_2 sy_2 sz_2; sx_3 sy_3 sz_3; ...], i is the axis along wich the spin is aligned
	S = Vector{Float64}([])
	for j in 1:n
		S = vcat(S,intial_Si(i,dir))
	end
	S_new = reshape(S,3,:)
	return Transpose(S_new)
end

function z_magnetisation(n, N, time_interval, number_repeats, Γ_deph, Γ_decay, Ω)
	#n is the number of spins, N is the number of time steps, time interval is given as (t_0, t_final)
	S_0 = initial_spin_array(n,3,-1)
	Sz_cummulative =  euler(N, time_interval, S_0, Γ_deph, Γ_decay, Ω)[1]
	for i in 2:number_repeats #we have already done the first one
		Sz_cummulative += euler(N, time_interval,initial_spin_array(n,3,-1), Γ_deph, Γ_decay, Ω)[1]
	end
	return Sz_cummulative/number_repeats
end