# Numerical Modeling of Electrokinetic Porous Networks
The developed pipeline allows numerical modeling of electrokinetic phenomena in porous networks.
A porous network is modeled as a random network of micro-scale and nano-scale pores.

# Running Code
Before running code, make sure that you have already installed lapack and blas library on your laptop/PC.
Needless to say, you need to have GCC compiler to compile the model that is developed in c++.

Once you have all the requirements ready, compile codes by running the Makefile:
make clean
make

If everything goes through, it means all things are installed properly and now you are ready to run the model!

# Input Files
To run the code you have to provide Network.in located at input_files folder. The guidline of how to prepare the file is in input_files folder.
The developed model uses a number of look-up tables to complete computations. These tables are already generated. You just need to put them in a directory and specify the path in the bash script. 
Once the input files are ready, a bash script is used to run the code, and it's located at bin folder.

# Toy Problem
As a starting point, the input file is designed for a network that involves a loop, as shown below:

