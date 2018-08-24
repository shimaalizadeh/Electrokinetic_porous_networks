This is the file that runs the whole pipeline<br/>

**codedir**: the location of the code. **DO NOT forget to set this with the right path to the code!**<br/>
**program**: the executive that is run by this script (**DO NOT CHANGE!**)<br/>

Please set the rest of the parameters based on the network you want to simulate:<br/>

**num_pores**: number of pores you have in your network<br/>
**num_reservoirs**: number of reservoirs (nodes) you have in the netwrok.<br/>
**h_pores**: set to 1 (**DO NOT CHANGE, it will be deprecated later**)<br/>
**dt**: the time step used by the model<br/>
**tmax**: maximum time of the simulation <br/>
**restart**: a boolean (0 or 1) **0** means simulation starts from t=0, **1** means simulation is continued from where it was stopped,
and restart files will be used instead of the original Network.in file.<br/>
**table_path**: the path to the look-up tables. **Please change this to where you put the tables**.<br/>
**network_input_file**: put the correct name to the input file located at "input_files" folder. If restart=0, the model will look for "network_input_file".<br/>

After you enter the right inputs, save and close the bash file. Then simply run:<br/>
bash run.sh (any_name_you_put.sh)<br/>

This should start the simulation. As an output for now model prints the exit current of the system towards the selective membrane.<br/>
