
function update(S::Vector,γ, g, dt)
	Sx_0=deepcopy(S[1])
	Sy_0=deepcopy(S[2])
	Sz_0=deepcopy(S[3])
	Sx = Sx_0 - (γ/2)*Sx_0*dt
	Sy = Sy_0 + (-(g*Sz_0)-((γ/2)*Sy_0))*dt
	Sz = Sz_0 + ((g*Sy_0)-((γ/2)*Sz_0))*dt
	return Sx,Sy,Sz
end

function spin_at_T(T,N,S,γ,g)
	dt = T/N
	spins = []
	Sx = S[1]; Sy = S[2]; Sz = S[3]
	for i in 1:N
		Sx, Sy, Sz = update([Sx,Sy,Sz],γ, g, dt)
		push!(spins,[Sx,Sy,Sz])
	end
	return spins
end

function different_gamma(S::Vector,m,g,γ,T,N)
	#S_1 = [1.0,1.0,-1.0]
	#S_2 = [-1.0,1.0,-1.0]
	#S_3 = [1.0,-1.0,-1.0]
	#S_4 = [-1.0,-1.0,-1.0]
	spins_1 = spin_at_T(T,N,S[1],γ,g)
	spins_2 = spin_at_T(T,N,S[2],γ,g)
	spins_3 = spin_at_T(T,N,S[3],γ,g)
	spins_4 = spin_at_T(T,N,S[4],γ,g)
	Sz_1=[-1.0]
	Sz_2=[-1.0]
	Sz_3=[-1.0]
	Sz_4=[-1.0]
	for i in 1:length(spins_1)
		push!(Sz_1,spins_1[i][3])
		push!(Sz_2,spins_2[i][3])
		push!(Sz_3,spins_3[i][3])
		push!(Sz_4,spins_4[i][3])
	end	
	return Sz_1, Sz_2, Sz_3, Sz_4
end