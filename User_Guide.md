# User Manual for Matlab-based 'Coral Bleaching Model'. 
This is a Matlab version of the coral bleaching model contained within the C-based 
CSIRO Environmental Modelling Suite 
[(GitHub - csiro-coasts/EMS: Environmental Modelling Suite)](https://github.com/csiro-coasts/EMS).

In order to run this model the entire repository should be downloaded from GitHub:
[(GitHub - gbrrestoration/coral-bleaching-model)](https://github.com/gbrrestoration/coral-bleaching-model). The repository contains two Matlab scripts and a Matlab function:
- **Run_Bleach** - script that specifices the configuration of the model (initial conditions, forcings and parameters) and manages the run.
- **Bleach** - function that contains the model equations. It is called by Run_Bleach.
- **Plot_Bleach** - script that plots common model outputs after the model integration is complete. This script is not necessary but is included as an aid for initial users.

The 'Bleach' function can be added to the Matlab path by 'running' the function and accepting 'add to path' this then allows the 'Run_Bleach' script to call this function as necessary. 

Time in the model is specified in days and the limits are set in the variable 'tspan'.

Inputs to the model are split into three categories: initial conditions, forcings, parameters.

**Initial Conditions:**
Initial conditions for all are both the input and output variables of the model.

**Forcings:**
Primary forcing inputs can be input in the following manner:
Temperature is entered as a timeseries of temperature anamoly (degC over the Maximum Monthly Mean). This is entered as an array with the 1st row being time in days and 2nd row being the temperature. The temperature at any given time is calulated using a linear interpolation function within 'Bleach'. Light is a calculated timeseries using the variable 'Peak' as the maximum daily sunlight and then fitting a  curve around midday. This timeseries is calculated hourly for the duration of 'tspan'. The light at any given time is calulated using a linear interpolation function within 'Bleach'. The remainder of the forcing inputs are entered as individual constants which remain consistent throughout the  model (i.e. are not timeseries).

An option exists in which the forcings data, specifically temperature, light, DIN and DIP can be read in from  files in order to replicate real world data such as from tank experimenets. In order to enable this option  the user must specify file_path, SDate and MMM. These variables are included in 'RunBleach' but are commented
out for regular usage of the model.
'file_path' = The local directory in which the relevant data is stored.
'SDate' = The date which will correspond to time 0 in the model (note that data should exist before this date ideally)
'MMM' = the Maximum Monthly Mean of particular corals being modelled (this is used to calculate the temperature anamoly)
Data must be organised into seperate csv files for each forcing variable being added. Additionally the name (including capitalization) of the file must match the name of the variable to be used (i.e. Temperature.csv, DIN_w.csv, DIP_w.csv, Light.csv).
Note that these variables should NOT be removed or commented out from earlier in the script.

**Parameters:**
Input parameters are all input by the user as singular constants. These values do not need to be changed however have been included in the 'Run_Bleach' script in order to allow for experimentation. 
