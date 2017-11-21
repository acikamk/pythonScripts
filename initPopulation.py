import math
import numpy
import sys
import os
'''
Script that generates popSize random vectors
that can be used as initial population for simtools
optimization. 
Inputs:
	- popSize - size of population, number of random vectors
	- nK - number of variables in the cost function
	- nKCD - number of variable parameters in the model 
			= Variable K dim: in the model joblog.txt
'''
def main(popSize, nK, nKCD):

	f = open("population_human.txt", "w")
	boundsK = (-1,1)
	boundsKCD = (1.e-4, 1)
	for i in range(int(popSize)):
		vectorKCD = boundsKCD[0] + numpy.random.rand(int(nKCD))	\
			*(boundsKCD[1] - boundsKCD[0]) 	
		vectorK = boundsK[0] + numpy.random.rand(int(nK)) \
			*(boundsK[1] - boundsK[0])  
		f.write("VECTOR: {}\n".format(numpy.concatenate((vectorKCD, vectorK)).tolist()))
	f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
    #notch 100 187 4 *Variable K dim: 187.
    #horst 300 17874 1 *Variable K dim: 17874.
    #human 300 4821 13