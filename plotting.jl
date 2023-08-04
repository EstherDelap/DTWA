using DTWA
using Plots
using JLD2, FileIO

function load_data(file_name, number, intro::String= "")
	#number is a range, filename is a string with a slash at the end
	#loads the data in file_name (from your local computer). file_name should be a file directory e.g. "/Documents/UCL_summer_internship/RESULTS/results_3/"
	#The content of this file need number of .jld2 files, all of the format "intro""number"".jld2"
	#For example, if number =4 and intro="take1_" then incside filename there should be: "take1_1.jld2", "take1_2.jld2", "take1_3.jld2", "take1_4.jld2"
	data=Vector{Dict}(undef,size(number)[1])
	for (i, x) in enumerate(number)
		data[i] = load(file_name*intro*string(x)*".jld2")
	end
	return data
end



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

