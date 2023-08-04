function magnetisation(average, axis)
	#axis refers to 1 for x, 2 for y or 3 for z, depending along where you want to measure magnetisation
	return getindex.(average[:],axis)
end