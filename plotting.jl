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

function plot_magnetisation_loading_data(file_name, number, intro::String= "")
	#plot(plot_magnetisation_not_loading(...)) to acc plot the data
	#plots the magnetisation if the data has been loaded into the file filename. number is a range, and file_name mist end in a slash
	data = load_data(file_name, number, intro) #data is an array of dictionaries
	
	plots = Vector{Vector{Float64}}(undef,size(data)[1])
	label = Array{String}(undef, 1,size(data)[1])
	for (i,dict) in enumerate(data)
		plots[i] = DTWA.magnetisation_3D(dict["average"],dict["axis"])
		label[1,i] = string("Ω = ",dict["Ω"])
	end
	return plots, label
end

function plot_magnetisation(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, α, method, axis, dir)
	#plot(plot_magnetisation_not_loading(...)) to acc plot the data
	#dim is a vector of length 3, N is the number of time intervals, α is a vector of length 3 that, for the Ising model, is just
	#one integer repeated three times, and for the XYZ model is [Jx,Jy,Jz] and must be floats. method is either "Ising" or "XYZ".
	#axis is the axis along which the spins are initially aligned ie 1, 2 or 3 for x, y or z respectively, and dir is the direction (-1 or +1)
	data = DTWA.repeated_euler(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, α, method, axis, dir)[2] #this is the average spin over time
	P =plot(DTWA.magnetisation_3D(data, axis))
	return P
end

function plot_spin_sqeezing_loading_data(file_name, number, intro::String= "")
	#to plot the data, do: plot(plot_spin_sqeezing_loading_data(...)[1], label = plot_spin_sqeezing_loading_data(...)[2])
	#plots the spin squeezing parameter if the data has been loaded into the file filename. number is a range, and file_name must end in a slash
	data = load_data(file_name, number, intro) #data is an array of dictionaries
	dim= data[1]["dim"]
	plots = Vector{Vector{Float64}}(undef,size(data)[1])
	labels = Array{String}(undef,1, size(data)[1])
	for (i,dict) in enumerate(data)
		plots[i]= -10*log.(10,DTWA.spin_sqeezing_param(dict["collective_spin"], dict["average"], dim[1]*dim[2]*dim[3]))
		labels[1,i] = string("α=", dict["α"][1])
	end
	return plots, labels
end

function plot_spin_sqeezing(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, α, method, axis, dir)
	#returns the plot for only a single simulation, so do plot(plot_spin_sqeezing_not_loading(...))
	#dim is a vector of length 3, N is the number of time intervals, α is a vector of length 3 that, for the Ising model, is just
	#one integer repeated three times, and for the XYZ model is [Jx,Jy,Jz] and must be floats. method is either "Ising" or "XYZ".
	#axis is the axis along which the spins are initially aligned ie 1, 2 or 3 for x, y or z respectively, and dir is the direction (-1 or +1)
	collectvie_spin,average = DTWA.repeated_euler(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, α, method, axis, dir)
	p= -10*log.(10,DTWA.spin_sqeezing_param(collectvie_spin, average, dim[1]*dim[2]*dim[3] ))
	return p
end
