using DTWA
using Plots
using JLD2, FileIO


function plot_different_gamma(S::Vector{Vector{Float64}};m=3, g = 1.0, γ = g/m, T = 10/γ,N = 100.0 )
    t= range(0,stop = T,step = T/N)
    S1,S2,S3,S4 = DTWA.different_gamma(S,m,g, γ, T,N)
	P1 = plot(t,[S1,S2,S3,S4], title=string("Sz(t) for g/γ=",m))
	xlabel!("t")
	ylabel!("Sz")
	P2 = plot(t,sum(+,(S1,S2,S3,S4))/4, title = string("<Sz(t)> for g/γ=",m))
	xlabel!("t")
	ylabel!("Sz")
	return P1, P2
end

function plot_magnetisation(n, N, time_interval, number_repeats,Γ_deph, Γ_decay, Ω)
	P =plot(DTWA.z_magnetisation(n, N, time_interval, number_repeats,Γ_deph, Γ_decay, Ω))
	return P
end

#function plot_spin_sqeezing_3D(dim,N, number_repeats,Γ_deph, Γ_decay,Ω)
	#P1 = DTWA.spin_sqeezing_param_3D(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, 0)
	#P2 = DTWA.spin_sqeezing_param_3D(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, 1)
	#P3 = DTWA.spin_sqeezing_param_3D(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, 2)
	#P4 = DTWA.spin_sqeezing_param_3D(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, 3)
	#P5 = DTWA.spin_sqeezing_param_3D(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, 4)
	#P6 = DTWA.spin_sqeezing_param_3D(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, 5)
	#P = plot([-10*log.(10,P1), -10*log.(10,P2), -10*log.(10,P3), -10*log.(10,P4), -10*log.(10,P5), -10*log.(10,P6)])
	#P = plot([-10*log.(10,P1),-10*log.(10,P3), -10*log.(10,P5)])
	#return P
#end

#function plot_magnetisation_3D(dim, N, time_interval, number_repeats, Γ_deph, Γ_decay, Ω, α)
	#P =plot(DTWA.z_magnetisation_3D(dim, N, time_interval, number_repeats, Γ_deph, Γ_decay, Ω, α))
	#return P
#end

function plot_spin_squeezing_2(file_dir, number_spins)
	data = load(file_dir,"collective_spin")
	averag = load(file_dir,"average")
	p = plot(DTWA.spin_sqeezing_param_3D(data, averag, number_spins))
end

function plot_magnetisation_2(file_dir, axis)
	#axis is 1,2 or 3 (ie x, y, or z) and file_dir is the file directory of the .jld2 file
	average = load(file_dir,"average")
	p = plot(DTWA.magnetisation_3D(average, axis))
end

