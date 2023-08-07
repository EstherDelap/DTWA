# Realistic simulations of open quantum systems using DDTWA
[![Runtests](https://github.com/EstherDelap/DTWA/actions/workflows/Runtests.yml/badge.svg)][def]

[def]: https://github.com/EstherDelap/DTWA/actions/workflows/Runtests.yml

## Description
Julia scripts for simulating open, interacting quantum systems using the Disipative Discrete Truncated Wigner Approximation. This is based on Peter Rabl's paper "Realistic simulations of spin squeezing and cooperative coupling effects in large ensembles
of interacting two-level systems". Interactions of the system with the environment are modeled as local decay and dephasing. This is achieved by replacing the mean field equations used in the Discrete Truncated Wigner Approximation by stochastic trajectories. The two interaction hamiltonians that can currently be used in this script are the transverse Ising Hamiltonian and the dissipative XYZ Hamiltonian. The time evolution of the average magnetisation can be calculated and plotted, as can the time evolution of the spin squeezing parameter. Of particular interest is simulating the XYZ model with different Jx, Jy and Jz parameters.

## Instalation

Download [julia](https://julialang.org/) (at least version 1.9) and clone my repository. There is also a julia extension to [VS Code](https://code.visualstudio.com/docs/languages/julia) (for which you still need julia downloaded). Then in the terminal (once in the desired file directory) clone the directory by writing "git clone (http or ssh code)".

## Usage

### Time evolution of a spin system

To generate the time evolution of the collective spin for a time interval time_interval, a total number of timesteps N and an initial spin array S_0:
```julia
DTWA.euler_3D(N, time_interval, S_0, Γ_deph, Γ_decay, Ω, Jx, Jy, Jz) 
```
The collective spin is $$\left[\sum_i^n s_i^x/2, \sum_i^n s_i^y/2, \sum_i^n s_i^z/2\right]$$ so this function returns a vector of length N, each element being  vector of length 3. Γ_deph and Γ_decay are the dephasing and decay rates respectively, and Ω is the driving field term. S_0 needs to be a 3 dimensional array where each element is a vector of length 3, containing [s_i^x, s_i^y, s_i^z] for that position in the array. 

To produce a 3D array of spins, each aligned in direction 'dir' along 'axis':
```julia
DTWA.spin_array_3D(dim, axis, dir)
```
dim must be a vector of integers of length 3, axis is either 1, 2 r 3 for x, y or z respectively and dir is +1 or -1. Jx, Jy, Jz are the interaction matrices. To produce interaction matrices for the transverse Ising model, use:
```julia
DTWA.Jx_Ising(dim)
DTWA.Jz_Ising(dim,α)
```
The first function is used for both Jx and Jy and produces an array of zeros, since spins don't interact along those directions in the transverse Ising model. The second produces interactions between spins that depend on their distance, according to the equation $$J_{ij}=\frac{1}{N} \frac{1}{|r_i - r_j|^α}$$


To produce interaction matrices for the XYZ model, where spins interact along the x, y and z directions with strength Jx, Jy and Jz respectively:
```julia
DTWA.J_XYZ(dims,j)
```
where j is either jx, jy or jz.

### Multiple time evolutions of a spin system

To generate the time evolution of the collective spin for multiple trajectories:
```julia
DTWA.repeated_euler(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, α, method, axis, dir)
```
This calls DTWA.euler_3D a number_repeats number of times. Every time it calls euler_3D it produces a new initial spin array S_0 using DTWA.spin_array_3D where the spins are aligned along axis in direction dir. 'method' is either "Ising" or "XYZ" so that the function produces the correct interaction matrices to go as an argument in DTWA.euler_3D. α is a vector of length 3. If Ising is the model used, it is 3 of the same integers, and if XYZ is the model used, it is 3 floating point numbers, [jx, jy. jz]. Returns a vector of lengths numner_repeats, each elemnt of which contains the contecnts of DTWA.euler_3D for a single trajectory with a randomly generated initial spin array, aligned in a specific direction. So each element of the vector returned by repeat_euler is a vector itself of length N decribing the time evolution of the average spins [sx,sy,sz]. DTWA.repeated_euler also returns average, which is a vector of length N, and is the time evolution of the average [sx,sy,sz], averaged over all the trajectories. 

### Plotting Magnetisation and the spin squeezing parameter

To plot the time evolution of the magnetisation of a spin system, we do:
```julia
plot_magnetisation(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, α, method, axis, dir)
```
This calls DTWA.repeated_euler and uses the second value it returns, the time evolution of [sx,sy,sz] spins averaged over all spins and all trajectories. The time evolution of the magnetisation along 'axis' is then given by
```julia
DTWA.magnetisation(average, axis)
```
where axis is 1, 2 or 3 to indicate the x, y or z direction. This is called by plot_magnetisation, which then reuturns the plot of the result.

The spin squeezing parameter is defined as $$\xi ^2 = min_{\phi}\left( \Delta S_{\phi}^{perp} \right)^2 \times \frac{N}{|\langle \vec{S} \rangle |^N}$$, whose full definition and use can be read in Peter Rable's paper ["Realistic simulations of spin squeezing and cooperative coupling effects in large ensembles of interacting two-level systems"](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.105.013716). To plot the time evolution of the spin squeezing parameter, one can do
```julia
plot_spin_sqeezing(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, α, method, axis, dir)
```
which, similarly to plotting the magnetisation, calls DTWA.repeated_euler to get the data needed to calculate the spin sqeesing parameter using the function
```julia
DTWA.spin_sqeezing_param_3D(collective_spin, average_spin, number_spins)
```
Which uses the optim module to get the minimum value in the spin squeezing parameter. This returns a vector of length N of the time evolution of the spin squeezing parameter. plot_spin_sqeezing then returns a plot of this.

### Using a cluster

You can then use DTWA.repeated_euler on a cluster. This is particularly usefull for running a batch script, for example when you want to produce multiple plots on one graph where one parameter varies. An example of the bash script needed to run the file is:

```bash
#!/bin/bash -l

# Wall clock time
#$ -l h_rt=72:00:00

# Request a single core
#$ -pe smp 1

# Request RAM
#$ -l mem=150G

# Request TMPDIR space
#$ -l tmpfs=500G


# Set up the job array.  In this instance we have requested 7 tasks equal to the number of lines in the parameter file
#-t 1-7

#$ -o output

#$ -e errors

# Set working dir
#$ -wd your/working/direcory

# Job name
#$ -N your_job_name

# Parse parameter file to get variables.
number=$SGE_TASK_ID

BASE=$PWD
paramfile=$BASE/input/params.txt
SAVE=$Base/results

module load julia/1.9.1
filename="`sed -n ${number}p $paramfile | awk '{print $1}'`"
v1="`sed -n ${number}p $paramfile | awk '{print $2}'`" #dim
v2="`sed -n ${number}p $paramfile | awk '{print $3}'`" #dim
v3="`sed -n ${number}p $paramfile | awk '{print $4}'`" #dim

v4="`sed -n ${number}p $paramfile | awk '{print $5}'`" #N
v5="`sed -n ${number}p $paramfile | awk '{print $6}'`" #number_repeats
v6="`sed -n ${number}p $paramfile | awk '{print $7}'`" #Γ_deph 
v8="`sed -n ${number}p $paramfile | awk '{print $9}'`" #Ω
v7="`sed -n ${number}p $paramfile | awk '{print $8}'`" #Γ_decay

v9="`sed -n ${number}p $paramfile | awk '{print $10}'`" #α
v10="`sed -n ${number}p $paramfile | awk '{print $11}'`" #α 
v11="`sed -n ${number}p $paramfile | awk '{print $12}'`" #α

v12="`sed -n ${number}p $paramfile | awk '{print $13}'`" #method
v13="`sed -n ${number}p $paramfile | awk '{print $14}'`" #axis 
v14="`sed -n ${number}p $paramfile | awk '{print $15}'`" #dir

cd $TMPDIR
cp $BASE/main.jl $TMPDIR/
julia --project=$BASE main.jl $v1 $v2 $v3 $v4 $v5 $v6 $v7 $v8 $v9 $v10 $v11 $v12 $v13 $v14
mv results.jld2 $SAVE/$filename

```

The julia script for this file, called main.jl, is:

```julia
using DTWA
using JLD2

function main(dim, N, number_repeats, Γ_deph, Γ_decay, Ω, α1,α2, α3, method, axis, dir)
    if method === "Ising"
        a = parse(Int64,α1)
        b = parse(Int64,α2)
        c = parse(Int64,α3)
    else
        a = parse(Float64,α1)
        b = parse(Float64,α2)
        c = parse(Float64,α3)
    end
    outfile = "results.jld2"
    rv1, rv2= DTWA.repeated_euler(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, [a,b,c], method, axis, dir)
    jldsave(outfile;dim=dim, α=[a,b,c], Γ_deph=Γ_deph, Γ_decay=Γ_decay, Ω = Ω, axis = axis, dir = dir, collective_spin = rv1, average = rv2)
end

arg1 = parse(Int64,ARGS[9])
arg2 = parse(Int64,ARGS[10])
arg3 = parse(Int64,ARGS[11])

main([parse(Int64,ARGS[1]),parse(Int64,ARGS[2]),parse(Int64,ARGS[3])],parse(Int64,ARGS[4]),parse(Int64,ARGS[5]), parse(Float64,ARGS[6]), parse(Float64,ARGS[7]), parse(Float64,ARGS[8]), ARGS[9], ARGS[10], ARGS[11], ARGS[12], parse(Int64,ARGS[13]), parse(Int64,ARGS[14]))
```
The paramter file, called params.txt in the filder input, is:
```txt
take1_1.jld2 4 4 4 1000 1000 0.0 0.0025 0 0 0 0 Ising 1 1
take1_2.jld2 4 4 4 1000 1000 0.0 0.0025 0 1 1 1 Ising 1 1
take1_3.jld2 4 4 4 1000 1000 0.0 0.0025 0 2 2 2 Ising 1 1
take1_4.jld2 4 4 4 1000 1000 0.0 0.0025 0 3 3 3 Ising 1 1
take1_5.jld2 4 4 4 1000 1000 0.0 0.0025 0 4 4 4 Ising 1 1
take1_6.jld2 4 4 4 1000 1000 0.0 0.0025 0 5 5 5 Ising 1 1
take1_7.jld2 4 4 4 1000 1000 0.0 0.0025 0 6 6 6 Ising 1 1
```
Running submit.sh will run the file main.jl 7 times (because of the jobscript). Each time main.jl is run, the input parameters correspond to the line in the parameter file, achieved with the command 'sed -n ${number}p' which chooses only the line corresponding to 'number' which is defined as the number of the jobscript (ie starts at 1 and ends at 7). Each line is then split into its different arguments, separated by spaces, using the command awk.  main.jl then calls DTWA.repeated_euler, after first converting the arguments to the correct format. The line jldsave then saves the result of repeated_euler to a file 'outfile' as a dictionary. The way 'outfile' is defined means that the results from the different lines in the parameter file can be saved in different files, all in the directory  $SAVE. 

To load the results from a file, lets say its filename.jld2, we use the FileIO method 'load':
```julia
data = load(filename)
```
To load just the time evolution for each trajectory of the spin values, ie the results of DTWA.repeated_euler, we can do
```julia
collective_spin = load(filename, "collective_spin")
```
To load all the data from the jobarray, so in this case all 7 files saved in the directory $Base/results, we use the function
```julia
load_data(file_name, number, intro::String= "")
```
file_name is a string refering to the file directory containing all the files and must have a slash at the end, so in this case file_name is "$Base/results/". THe files in $Base/results/ must be saved as "(intro)(number).jld2", so in this example the files are saved as "take1_1.jld2". number is a range, refering to the number of files you want to plot from the directory "$Base/results/" and intro is an optional string. In this example, intro is "take1_" and number is range(1,7) (if I wanted to plot all 7 lines). The function load_data then loops through the range of numbers and calls load("(intro)(number).jld2") for all the integers in the number range. It saves these results into an array the same length as the range numbers, with each value being a dictionary of results. 

To plot the magnetisation of these trajectories, we have created a function:
```julia
plot_magnetisation_loading_data(file_name, number, intro::String= "")
```
where file_name, number and intro are the arguements that go into load_data. Inside this function, the function DTWA.magnetisation_3D() is then called for each dictionary in the list, with the 'average' argument being dict["average"] and the 'axis' argument being dict["axis"]. This will then return two arguments, the first a list of plots, the second a list of labels for those plots that indicates the different $\omega$ values. To use this, one could do:
```julia
begin
    plots, labels = plot_magnetisation_loading_data(file_name, number, intro)[1]
    plot(plots, label = labels)
end
```

To plot the spin squeezing parameter, we have created a very similar function:
```julia
plot_spin_sqeezing_loading_data(file_name, number, intro::String= "")
```
which similarly uses file_name, number and intro as the arguments to load_data and then calls DTWA.spin_sqeezing_param_3D(collective_spin, average_spin, number_spins) with collective_spin being dict["collective_spin"], the average_spin being dict["average] and number_spins being dim[1]*dim[2]*dim[3] where dims is dict["dim"]. Again, two arguments are returned, the plots of the spin squeezing parameter and the labels indicating this time the different $\alpha$ values. To use this, one could do:
```julia
begin
    plots, labels = plot_spin_sqeezing_loading_data(file_name, number, intro)[1]
    plot(plots, label = labels)
end
```

