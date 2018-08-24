This is the file that runs the whole pipeline

codedir: the location of the code. DO NOT forget to set this with the right path to the code!
program: the executive that is run by this script (DO NOT CHANGE!)

Please set the rest of the parameters based on the network you want to simulate:
num_pores: number of pores you have in your network
num_reservoirs: number of reservoirs (nodes) you have in the netwrok. 
h_pores: set to 1 (DO NOT CHANGE, it will be deprecated later)
dt: the time step used by the model
tmax: maximum time of the simulation 
restart: a boolean (0 or 1) 0 means simulation starts from t=0, 1 means simulation is continued from where it was stopped,
and restart files will be used instead of the original Network.in file.
table_path: the path to the look-up tables
network_input_file: put the correct name to the input file located at "input_files" folder. If restart=0, the model will look for
"network_input_file". 


After you enter the right inputs, save and close the bash file. Then simply run:
bash run.sh (any_name_you_put.sh)

This should start the simulation. As an output for now model prints the exit current of the system towards the selective membrane.
