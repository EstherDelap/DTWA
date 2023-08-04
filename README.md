# Realistic simulations of open quantum systems using DDTWA
[![Runtests](https://github.com/EstherDelap/DTWA/actions/workflows/Runtests.yml/badge.svg)][def]

[def]: https://github.com/EstherDelap/DTWA/actions/workflows/Runtests.yml

## Description
Julia scripts for simulating open, interacting quantum systems using the Disipative Discrete Truncated Wigner Approximation. This is based on Peter Rabl's paper "Realistic simulations of spin squeezing and cooperative coupling effects in large ensembles
of interacting two-level systems". Interactions of the system with the environment are modeled as local decay and dephasing. This is achieved by replacing the mean field equations used in the Discrete Truncated Wigner Approximation by stochastic trajectories. The two interaction hamiltonians that can currently be used in this script are the transverse Ising Hamiltonian and the dissipative XYZ Hamiltonian. The time evolution of the average magnetisation can be calculated and plotted, as can the time evolution of the spin squeezing parameter. Of particular interest is simulating the XYZ model with different Jx, Jy and Jz parameters.

## Instalation

Download [julia](https://julialang.org/) (at least version 1.9) and clone my repository. There is also a julia extension to [VS Code](https://code.visualstudio.com/docs/languages/julia) (for which you still need julia downloaded). Then in the terminal (once in the desired file directory) clone the directory by writing "git clone (http or ssh code)".

## Usage

To generate the time evolution of the collective spin for a time interval time_interval, a total number of timesteps N and an initial spin array S_0:
```julia
DTWA.euler_3D(N, time_interval, S_0, Γ_deph, Γ_decay, Ω, Jx, Jy, Jz) 
```
The collective spin is $$[\sum_i^n s_i^x/2, \sum_i^n s_i^y/2, \sum_i^n s_i^z/2]$$ so this function returns a vector of length N, each element being  vector of length 3. Γ_deph and Γ_decay are the dephasing and decay rates respectively, and Ω is the driving field term. S_0 needs to be a 3 dimensional array where each element is a vector of length 3, containing [s_i^x, s_i^y, s_i^z] for that position in the array. To produce a 3D array of spins, each aligned in direction 'dir' along 'axis':
```julia
DTWA.spin_array_3D(dim, axis, dir)
```
dim must be a vector of integers of length 3, axis is either 1, 2 r 3 for x, y or z respectively and dir is +1 or -1. Jx, Jy, Jz are the interaction matrices. To produce interaction matrices for the transverse Ising model, use:
```julia
DTWA.Jx_Ising(dim)
DTWA.Jz_Ising(dim,α)
```
The first function is used for both Jx and Jy and produces an array of zeros, since spins don't interact along those directions. The second produces interactions between spins that depend on their distance, according to the equation $$J_{ij}=\frac{1}{N} \frac{1}{|r_i - r_j|^α}$$
To produce interaction matrices for the XYZ model, where spins interact along the x, y and z directions with strength Jx, Jy and Jz respectively:
```julia
DTWA.J_XYZ(dims,j)
```
where j is either jx, jy or jz.

To generate the time evolution of the collective spin for multiple trajectories:
```julia
DTWA.repeated_euler(dim, N,number_repeats,Γ_deph, Γ_decay,Ω, α, method, axis, dir)
```
This calls DTWA.euler_3D a number_repeats number of times. Every time it calls euler_3D it produces a new initial spin array S_0 using DTWA.spin_array_3D where the spins are aligned along axis in direction dir. 'method' is either "Ising" or "XYZ" so that the function produces the correct interaction matrices to go as an argument in DTWA.euler_3D. α is a vector of length 3. If Ising is the model used, it is 3 of the same integers, and if XYZ is the model used, it is 3 floating point numbers, [jx, jy. jz]. 

You can then use DTWA.repeated_euler in a HPC facilitiy. This can be very useful if you want to plot multiple trajectories on one graph. If this is done, then save the results in a folder, and inside this folder have the first result saved in a .jld2 file as "pre_name_1.jld2". "pre_name" can be any string (eg "take1_") as long as it's the same for all files you wish to plot. Save the second result as "pre_name_2.jld2" the third as "pre_name_3.jld2" etc. This will then allow you to use the function; 
```julia
DTWA.load_data(file_name, number, intro)
```
to load the data as a vector of dictionaries. This function makes use of the FileIO function load(file_name) which loads the content of the file at the path file_name as a dictionary. This is achieved by saving the information from DTWA.repeated_euler using:
```julia
jldsave(outfile;dim=dim, α=[a,b,c], Γ_deph=Γ_deph, Γ_decay=Γ_decay, Ω = Ω, axis = axis, dir = dir, collective_spin = rv1, average = rv2)
```
done inside the HPC. So this means that load("home\Documents\...\take1_1.jld2")["collective_spin"] will give me the collective spin. Thus, the function DTWA.load_data is an array the length of the number of spin systems you wish to plot, each element being a dictionary containing informtation about the time evolution of one of your spin systems.


