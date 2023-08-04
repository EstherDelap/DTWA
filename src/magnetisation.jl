function magnetisation(average, axis)
	#Finds the time evolution of the magnetisation ie average sx, sy or sz. Axis refers to 1, 2 or 3 for x, y or z, respectively.
	return getindex.(average[:],axis)
end