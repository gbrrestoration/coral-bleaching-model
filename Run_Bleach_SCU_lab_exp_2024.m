%% Run bleach model

% developed by Mark Baird and Julie Terp Jorgensen (2021).

% Description: 
% This process contains a coral model with autotrophic and heterotrophic growth, zooxanthellae physiology, 
% xanthophyll cycle, reaction centre dynamics and reative oxygen build-up.

% Differences from coral bleaching model in eReefs 

% 1. Downwelling irradiance is not spectrally-resolved.
% 2. No temperature dependence of physiological rates (except hardwired through a_rub)
% 3. Host dynamics not working yet.
    
close all 
clear all

tspan = [0 26]; % Days (integer format)

% initial conditions:

CS_N0 = 0.4;
CH_N0 = 1000;
CS_Chl = CS_N0*5.6786/30/20; % 1/20th of maximum 

y0(1) = CS_N0; 
y0(2) = CS_N0*0.5;
y0(3) = CS_N0*(1*((1/16)*(30.97/14.01)))*0.5;
y0(4) = CS_N0*(((106/16)*(12.01/14.01)))*0.5;
y0(5) = CH_N0;

% The 15 entries where the first 4 are multiplied with CS_N0
% units: y(1): mg N, y(2): dimensionless, y(3): dimensionless, y(4): dimensionless

y0(6) = CS_Chl;
y0(7) = CS_Chl*0.2448*0.33;
y0(8) = CS_Chl*0.2448*0.67;

chla2rcii = 0.002/893.49;

y0(9) = CS_Chl * chla2rcii*0.05;
y0(10) = CS_Chl * chla2rcii*0.1;
y0(11) = CS_Chl * chla2rcii*0.85;
y0(12) = 0.0;

disp('Initial conditions: [CS_N Rn RP RC CH_N Chl Xp Xh Qox Qred Qin ROS]');

y0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input forcings here
% NOTE forcings are overwritten below if file_path is specified
% temperature timeseries, 1st row is times (in days) 2nd is temperature anamoly
Temperature = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26;
                    0 0 0 0.7 1.8 2.9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 2.79 1.2 0 0 0 0 0 ]; 

% Makes Light into a timeseries
Peak = 1000; %Peak daily sunlight (W/m2)
Light(1,:) = [tspan(1):1/24:tspan(2)];
Light(2,:) = Peak*max(0,sin((Light(1,:)-0.25-floor(Light(1,:)-0.25))*2*pi)); % W m-2 

DIN_w = 40.71428571; % (mg N m-3);
DIP_w = 13.57142857*(1*((1/16)*(30.97/14.01))); % (mg P m-3)       
POM = 0.0; % Particulate organic matter.
tau = 0.3; % sheer stress at the bottom (N m-2) - values between 0-1
rho = 1026; % density of the water (kg m-3)
Salinity = 35;% Note that Salinity is not currently used in this model.

% Input options for reading forcing data from files (overwrites above forcing variables)
% Comment out for normal model usage 
%file_path = 'C:/Users/Example/Bleaching_submodel/Test_Data/';
%SDate = datetime(2021,11,28,0,0,0); % Starting datetime [format = (Year,Month,Day,Hour,Minute,Second)]
%MMM = 26; %Maximum Monthly Mean for the corals you are modelling
% This requires each recorded variable to be saved in individual csv files
% in the filepath folder, with no other csv files.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters here
% Parameter values table A.7 in Baird et al. 2018
mN = 4.514805e-14*16*1000*14.01;  % Nitrogen content of zooxanthellae cells           (mg N cell-1)
rCS = 5e-6;                        % Radius of zooxanthellae cells                     (m)
CSvol = 4*pi/3*rCS^3;    % volume of the cell (m3)
uCH = 0.05;                     % Maximum growth rate of coral host                     (d-1)
Spart = 3.0;                    % Rate coefficient of particle capture              (m d-1)
uCS = 0.4;                      % Maximum growth rate of zooxanthellae              (d-1)
zetaCH = 0.0; %0.01;                  % Quadratic mortality coefficient of polyps         (d-1 (g N m-2)-1)
zetaCS = 0.0;                  % Linear mortality of zooxanthellae                 (d-1) 
CHremin = 0.5;                   % Remineralised fraction of coral mortality 
OCH = 2.0/1000;                      % Nitrogen-specific host area coefficient of polyps (m2 mg N-1)
C2Chlmin = 20;                  % max Chl:C ratio at which synthesis stops.    

Chlmax = 2.0*2.09e7 * ((1.0e18*CSvol)^-0.310);
Xanth_tau = 72; % 1 per 20 minutes %%72 times per day
photon2rcii = 0.1e-6;
ROSthreshold = 1.418e-14;
photon2ros = 7000.0;
chla2rcii = 0.002/893.49;
CSmaxbleachrate = 1.0; % d-1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If file_path is specified then read in forcings from a folder of csv files
if exist('file_path','var') == 1
    files = {dir([file_path,'*.csv'])};
    filenames = join([{files{1}.folder}; {files{1}.name}],"\",1);
    for i = 1:length(filenames)
        [pth, nme, ext] = fileparts(filenames{1,i}); % Extract name of variable
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
% Forcings are collated into forcings variable for passing to Bleach function
forcings = cell(8,1); % This must be the right number, if more forcing variables are added
forcings{1} = DIN_w;
forcings{2} = DIP_w;
forcings{3} = Temperature;
forcings{4} = POM;
forcings{5} = Salinity;
forcings{6} = tau;
forcings{7} = rho;
forcings{8} = Light;

% Check for non-finite and repeated values in forcings that will break 'Bleach' function
for i=1:length(forcings)
    if sum(sum(isfinite(forcings{i})<1)) > 0
        index = find(isfinite(forcings{i})<1);
        disp(strcat("Error: Forcings variable '",num2str(i),"' has a non finite value"))
            return
    end
    if length(forcings{i}(1,:)) ~= length(unique(forcings{i}(1,:)))
        [v, w] = unique( forcings{i}(1,:), 'stable' );
        duplicate_indices = setdiff( 1:numel(forcings{i}(1,:)), w );
        disp(strcat("Error: Forcings variable '",num2str(i),"' has a repeated time point ",datestr(datetime((forcings{i}(1,duplicate_indices)+datenum(SDate)),'ConvertFrom','datenum'))))
        return
    end
end

% Parameters are collated into forcings variable for passing to Bleach function
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',60/86400);
y00 = y0;
tstore = [];
ystore = [];
for tt = 1:tspan(2)
    [t,y] = ode45(@Bleach, [tt-1 tt], y00, opts, forcings, params);
    y00 = y(end,:);
    tstore = [tstore; t];
    ystore = [ystore; y];
    disp(['completed day = ',num2str(tt),' with ',num2str(length(t)),' steps']);
end

t = tstore;
y = ystore;

% calculations for plotting:

% absorption cross-section (value between 0 and pi r^2)

rCS = 5e-6;mN = 4.514805e-14*16*1000*14.01; 
cells = y(:,1) /mN;CSvol = 4*pi/3*rCS^3;
C2Chlmin = 20;
ROSthreshold = 1.418e-14;
CSmaxbleachrate = 1.0; % d-1
				  
cellXp = y(:,7)./(cells*CSvol);    % mg pig m-3
cellXh_= y(:,8)./(cells*CSvol);
cellChl = y(:,6)./(cells*CSvol);
				  
absorb =  cellXp * 0.06 + cellChl * 0.06;
temp = 2.0 * absorb * rCS;
aA = pi * rCS * rCS * (1.0 - 2.0 * (1.0 - (1.0 + temp) .* exp(-temp)) ./ (temp .* temp));

nChl2C = y(:,6)./y(:,1)/5.6786*C2Chlmin;

ROSpercell = y(:,12)./cells;
expel = CSmaxbleachrate*max(0.0,(ROSpercell - ROSthreshold)./ROSthreshold)*ROSthreshold;

rctotal = sum(y(:,9:11),2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
