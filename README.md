# Multi-Scale Modeling of Electrokinetic Porous Networks
The developed pipeline allows efficient modeling of electrokinetic phenomena in porous networks.
A porous network is modeled as a random network of micro-scale and nano-scale pores.

# Running Code
The modeling package is implemented in C++. You need C++ compiler (such as gcc) to compile the code and build the executable file.
Make sure that you have already installed lapack and blas library on your laptop/PC.

Once you have all the requirements ready, compile codes by running the Makefile:
make clean  
make

If everything goes through, it means all things are installed properly and now you are ready to run the model!

# Input Files
To run the code you have to provide Network.in located at ./input_files folder. Please refer to the guideline on how to prepare the input files at ./input_files folder.
The developed model uses a number of look-up tables to complete computations. These tables are already generated and located at the default folder ./tables.
If you put them i na different location, you just need to specify the path in the bash script in ./bin folder. 
Once the input files are ready, a bash script is used to run the code, and it's located at bin folder.

# Toy Problems
A few example problems were provided along with source code. You can find the corresponding input files 
under ./examples folder. 

[//]: # (the input file is designed for a network that involves a loop, as shown below:)

![alt text](https://github.com/shimaalizadeh/Electrokinetic_porous_networks/blob/master/toy_problem.png)


