- The first part if for pores' data. Note that each pore has an id from 0 to num_pores-1
There are 10 columns that are respectively:

Column1: type of pore (slit(0) or circular(1))<br/>
Column2: if the pore has dead-end(1) or not (0)<br/>
Column3: pore input reservoir id<br/>
Column4: pore output reservoir id<br/>
Column5: pore length (non-dimensional)<br/>
Column6: number of mesh points inside pore<br/>
Column7: sigma_star<br/>
Column8: lambda_D<br/>
Column9: pore area cross section (non-dimensional)<br/>
Column10: pore concentration (1 at t=0)<br/>

Example:<br/>
0 0 0 1 0.6 1600 -30 0.017 1.00 1.00<br/>
0 0 1 3 0.6 1600 -30 0.017 1.00 1.00<br/>
0 1 3 5 0.6 1600 -30 0.017 1.00 1.00<br/>
0 0 1 2 0.6 1600 -15 0.034 0.5 1.0<br/>
0 0 2 4 0.6 1600 -15 0.034 0.5 1.00<br/>
0 0 3 4 0.6 1600 -15 0.034 0.5 1.00<br/>

- The second part includes reservoirs data. Each reservoir has an id from 0 to num_nodes-1
There are 9 columns that are respectively:<br/>
Column1: reservoir type (slit(0) or circular(1))<br/>
Column2: If the pressure is known in the reservoir (1 or 0)<br/>
Column3: If the electrochemical potential is known in the reservoir (1 or 0)<br/>
Column4: value of pressure in the reservoir (put 0 if it's unknown)<br/>
Column5: value of electrochemical potential in the reservoir (put 0 if it's unknown)<br/>
Column6: value of concentration (1 at t=0)<br/>
Column7: sigma_star of reservoir<br/>
Column8: lambda_D of reservoir<br/>
Column9: volume of reservoir (0 if pressure or potential are both known, which means it's an end reservoir)<br/>

Example:    
0 0 0 0.0 40.0 1.00 -30 0.017 0.00<br/>
0 1 1 0.0 0.0 1.00 -30 0.017 8.068e-05<br/>
0 1 1 0.0 0.0 1.00 -15 0.034 8.068e-05<br/>
0 1 1 0.0 0.0 1.00 -30 0.017 5.072e-05<br/>
0 1 1 0.0 0.0 1.00 -15 0.034 7.602e-05<br/>
0 1 0 0.0 0.0 0.10 -30 0.017 0.00<br/>

- The third part of Network.in describe the network connectivity matrix, which says:
Reservoir i is connected to reservoir j by pore k. So there are three columns:<br/>
Column1: reservoir i<br/>
Column2: reservoir j<br/>
Column3: pore k<br/>

Example:  
0 1 0<br/>
1 3 1<br/>
1 2 3<br/>
2 4 4<br/>
3 4 5<br/>
3 5 5<br/>

- The last part of the file is no longer used.
1600 6 4 0.62 0.28 0.03 0.0 0.0 0.0 0.0<br/>
1600 6 4 0.62 0.19 0.03 0.0 -0.6 0.0 0.0<br/>
1600 6 4 0.62 0.25 0.03 0.0 -1.2 0.0 0.0<br/>
1600 6 4 0.62 0.26 0.03 0.0 -5.0 0.0 0.0<br/>
1600 6 4 0.62 0.33 0.03 0.0 -5.6 0.0 0.0<br/>
2 3 1 0.12 0.33 0.00 0.0 0.0 0.0 0.0<br/>
2 3 1 0.12 0.33 0.00 0.0 0.0 0.0 0.0<br/>
2 3 1 0.12 0.33 0.00 0.0 0.0 0.0 0.0<br/>
2 3 1 0.12 0.33 0.00 0.0 0.0 0.0 0.0<br/>
2 3 1 0.12 0.33 0.00 0.0 0.0 0.0 0.0<br/>
2 3 1 0.12 0.33 0.00 0.0 0.0 0.0 0.0<br/>
