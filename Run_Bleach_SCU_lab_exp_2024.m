%% Run bleach model

% Developed by Mark Baird and Julie Terp Jorgensen (2021).

% Description: 
% This process contains a coral model with autotrophic and heterotrophic growth, zooxanthellae physiology, xanthophyll cycle, reaction centre dynamics and reative oxygen build-up.

% This configuration of Run_Bleach was used to run the model for a laboratory experiment
  
close all 
clear all

tspan = [0 23]; % Days (integer format), adjust as necessary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter values (Table A.7 in Baird et al. (2018))

mN = 4.514805e-14*16*1000*14.01;  					% Nitrogen content of zooxanthellae cells (mg N cell-1)
rCS = 5e-6;                        					% Radius of zooxanthellae cells (m)
CSvol = 4*pi/3*rCS^3;    						% Volume of the cell (m3)
uCH = 0.05;                    						% Maximum growth rate of coral host (d-1)
Spart = 3.0;                    					% Rate coefficient of particle capture (m d-1)
uCS = 0.4;                      					% Maximum growth rate of zooxanthellae (d-1)
zetaCH = 0.0; %0.01;                  					% Quadratic mortality coefficient of polyps (d-1 (g N m-2)-1)
zetaCS = 0.0;                  						% Linear mortality of zooxanthellae (d-1) 
CHremin = 0.5;                   					% Remineralised fraction of coral mortality 
OCH = 2.0/1000;                      					% Nitrogen-specific host area coefficient of polyps (m2 mg N-1)
C2Chlmin = 20;                  					% max Chl:C ratio at which synthesis stops    

Chlmax = 2.0*2.09e7 * ((1.0e18*CSvol)^-0.310);
Xanth_tau = 72; 							% 1 per 20 minutes %%72 times per day
photon2rcii = 0.1e-6;
ROSthreshold = 1.418e-14;
photon2ros = 3500; 							% Initial condition for the number of photons that lead to the generation of one ROS was reduced in this configuration to represent a thermally sensitive coral species with a less efficient photosynthetic apparatus
chla2rcii = 0.002/893.49;
CSmaxbleachrate = 1.0; 							% (d-1)
repair_coefficient = 268;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial conditions:

% Where adjusted for laboratory experiment see comment
% See Run_Bleach.m for original initial conditions

% Initial conditions are passed to Bleach_SCU_lab_exp_2024.m and then re-calculated at each timepoint

CS_N0 = 1;								% Initial condition of symbiont biomass (mg N m-2) was increased to represent a healthy starting state
CH_N0 = 10;								% Initial condition of coral host biomass (g N m-2) was reduced to be more representative of a laboratory experiment
CS_Chl = CS_N0*5.6786/30/20; % 1/20th of maximum 			% Initial condition of symbiont chlorophyll a concentration (mg m-2) was increased to represent a healthy starting state

y0(1) = CS_N0; 
y0(2) = CS_N0*0.5;
y0(3) = CS_N0*(1*((1/16)*(30.97/14.01)))*0.5;
y0(4) = CS_N0*(((106/16)*(12.01/14.01)))*0.5;
y0(5) = CH_N0;

% The 15 entries where the first 4 are multiplied with CS_N0 
% Units: y(1): mg N, y(2): dimensionless, y(3): dimensionless, y(4): dimensionless

y0(6) = CS_Chl;
y0(7) = CS_Chl*0.2448*0.33;
y0(8) = CS_Chl*0.2448*0.67;

chla2rcii = 0.002/893.49;

y0(9) = 1.0607E-07;							% Initial condition for oxidised reaction centre concentration (mg m-2) obtained from a pre-experimental simulation 
y0(10) = 7.1695E-09;							% Initial condition for reduced reaction centre concentration (mg m-2) obtained from a pre-experimental simulation 
y0(11) = 1.0108E-07;							% Initial condition for inhibited reaction centre concentration (mg m-2) obtained from a pre-experimental simulation 
y0(12) = ROSthreshold*CS_N0/mN*0.5;					% Intitial condition of ROS concentration was set to half the ROS threshold (Suggett, 2008), to reduce the model spin-up time to a stable starting state

disp('Initial conditions: [CS_N Rn RP RC CH_N Chl Xp Xh Qox Qred Qin ROS]');

y0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input forcings here:

% Adjust as necessary
% NOTE forcings are overwritten below if file_path is specified

% Temperature timeseries
% 1st row is times (in days) 2nd is temperature anomaly:
Temperature = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24;
                    0 0 0 0 0 0 0.3 1.4 2 2.3 3.3 3.6 4.3 3.9 4 4 4.2 3.8 3.1 4 4 4 3.95 3.9 1.9]; 

% Make Light into a timeseries:
Peak = 300; 								% Peak daily sunlight (W/m2); average value supplied from experimental data
Light(1,:) = tspan(1):1/24:tspan(2);
Light(2,:) = Peak*max(0,sin((Light(1,:)-0.25-floor(Light(1,:)-0.25))*2*pi)); % W m-2 

% Nutrients:
DIN_w = 41.855; 							% (mg N m-3); average value supplied from experimental data
DIP_w = 17.9*(1*((1/16)*(30.97/14.01))); 				% (mg P m-3); average value supplied from experimental data       

POM = 0.0; 								% Particulate organic matter (mg m-3)
tau = 0.3; 								% Sheer stress at the bottom (N m-2) - values between 0-1
rho = 1026; 								% Density of the water (kg m-3)
Salinity = 35;								% Note that Salinity is not currently used in this model.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input options for reading forcing data from files (overwrites above forcing variables):

% This requires each recorded variable to be saved in individual .csv files in the filepath folder, with no other csv files.
% Comment out for normal model usage 

% file_path = 'C:/Users/Example/Bleaching_submodel/Test_Data/';		% Light, temperature and nutrient forcing data recorded in the laboratory experiment was supplied to the model

% SDate = datetime(2021,11,28,0,0,0); 					% Starting datetime [format = (Year,Month,Day,Hour,Minute,Second)], ensure forcings have a value before this datetime
% MMM = 28.6; 								% Maximum Monthly Mean of the corals in laboratory experiment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If file_path is specified then read in forcings from a folder of csv files:

if exist('file_path','var') == 1
    files = {dir([file_path,'*.csv'])};
    filenames = join([{files{1}.folder}; {files{1}.name}],"\",1);
    for i = 1:length(filenames)
        [pth, nme, ext] = fileparts(filenames{1,i}); 			% Extract name of variable
        dat = readtable(filenames{1,i});
        time = datenum(table2array(dat(:,1))) - datenum(SDate);
        data = table2array(dat(:,2));
        if exist(nme,'var') == 1
            disp(strcat("Finished reading in ",nme,", overwriting original variable"))
        else
            disp(strcat("Error: variable name from file does not match any existing variables: ",nme))
            return
        end
        assignin('base',sprintf(nme),[time.'; data.']);
        clear dat time data nme
    end
    Temperature(2,:) = Temperature(2,:) - MMM;
    Temperature(2,Temperature(2,:)<0) = 0;
end
clear files filenames pth ext i


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forcings are collated into forcings variable for passing to Bleach_SCU_lab_exp_2024.m function:

forcings = cell(8,1); 							% This must be the right number of forcings
forcings{1} = DIN_w;
forcings{2} = DIP_w;
forcings{3} = Temperature;
forcings{4} = POM;
forcings{5} = Salinity;
forcings{6} = tau;
forcings{7} = rho;
forcings{8} = Light;

% Check for non-finite and repeated values in forcings that will break 'Bleach_SCU_lab_exp_2024.m' function:

for i=1:length(forcings)
    if sum(sum(isfinite(forcings{i})<1)) > 0
        index = find(isfinite(forcings{i})<1);
        disp(strcat("Error: Forcings variable '",num2str(i),"' has a non finite value"))
            return
    end
    if length(forcings{i}(1,:)) ~= length(unique(forcings{i}(1,:)))	% If length of datetime is not equal to unique data points
        [v, w] = unique( forcings{i}(1,:), 'stable' );			% Extract the indices for duplicated points found
        duplicate_indices = setdiff( 1:numel(forcings{i}(1,:)), w );
        disp(strcat("Error: Forcings variable '",num2str(i),"' has a repeated time point ",datestr(datetime((forcings{i}(1,duplicate_indices)+datenum(SDate)),'ConvertFrom','datenum')))) % Returns an error message and stops the code, display replicated datetime 
        return
    end
end

% Parameters are collated into forcings variable for passing to Bleach_SCU_lab_exp_2024.m function:

params(1) = mN;
params(2) = rCS;
params(3) = CSvol;
params(4) = uCH;
params(5) = Spart;
params(6) = uCS;
params(7) = zetaCH;
params(8) = zetaCS;
params(9) = CHremin;
params(10) = OCH;
params(11) = C2Chlmin;
params(12) = Chlmax;
params(13) = Xanth_tau;
params(14) = photon2rcii;
params(15) = ROSthreshold;
params(16) = photon2ros;
params(17) = chla2rcii;
params(18) = CSmaxbleachrate;
params(19) = repair_coefficient; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display some initial conditions before the run:

% Some calculations:
Aeff = 1 - exp(-OCH * y0(5));						% Aeff is the effective projected area fraction of the coral community

normalisedprojectedarea =  CS_N0 / mN * pi * rCS * rCS/Aeff; 		% Normalised projected area of the symbiont (m2)

% Projected area of 1 cell (pi r2) and multiplied it by the number of cells / Aeff of the host
% Normalised means to the projected area of the polyp

% Output some settings:

disp(['Projected Area [Norm] = ',num2str(normalisedprojectedarea)])
disp(['Photons to ROS = ',num2str(photon2ros)])
disp(['Starting Zoothanthallae Chlorophyll = ',num2str(CS_Chl)])
disp(['Zoothanthallae Repair Coefficient = ',num2str(repair_coefficient)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',600/86400);		% Max step increased 
y00 = y0;
tstore = [];
ystore = [];
for tt = 1:tspan(2)
    [t,y] = ode45(@Bleach_SCU_lab_exp_2024, [tt-1 tt], y00, opts, forcings, params);
    y00 = y(end,:);
    tstore = [tstore; t];
    ystore = [ystore; y];
    disp(['completed day = ',num2str(tt),' with ',num2str(length(t)),' steps']);
end

t = tstore;
y = ystore;

% Calculations for plotting:

% Absorption cross-section (value between 0 and pi r^2)

rCS = 5e-6;mN = 4.514805e-14*16*1000*14.01; 
cells = y(:,1) /mN;CSvol = 4*pi/3*rCS^3;
C2Chlmin = 20;
ROSthreshold = 1.418e-14;
CSmaxbleachrate = 1.0; 							% (d-1)
				  
cellXp = y(:,7)./(cells*CSvol);    					% (mg pig m-3)
cellXh_= y(:,8)./(cells*CSvol);
cellChl = y(:,6)./(cells*CSvol);
				  
absorb =  cellXp * 0.06 + cellChl * 0.06;
temp = 2.0 * absorb * rCS;
aA = pi * rCS * rCS * (1.0 - 2.0 * (1.0 - (1.0 + temp) .* exp(-temp)) ./ (temp .* temp));

nChl2C = y(:,6)./y(:,1)/5.6786*C2Chlmin;

ROSpercell = y(:,12)./cells;
expel = CSmaxbleachrate*max(0.0,(ROSpercell - ROSthreshold)./ROSthreshold)*ROSthreshold;

rctotal = sum(y(:,9:11),2);

temperatures = interp1(Temperature(1,:),Temperature(2,:),t);

a_rub = min(1.0,max(0,(1-exp(-(2-temperatures)))/(1-exp(-2))));

Aeff = 1 - exp(-OCH * y(5));

normalisedprojectedarea =  cells * pi * rCS * rCS/Aeff; 		% (m2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
