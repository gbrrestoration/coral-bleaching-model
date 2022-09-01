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

tspan = [0 26];

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
% Enter in inputs here (lukes code fixers) 

DIN_w = 40.71428571; % (mg N m-3);
DIP_w = 13.57142857*(1*((1/16)*(30.97/14.01))); % (mg P m-3)

% temperature timeseries, 1st row is times (in days) 2nd is temperature
% anamoly at those times (this can be altered to come from a file)
temp_timeseries = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26;
                    0 0 0.7 1.8 2.9 4 4 4 4 4 4 4 4 4 4 4 4 4 4 2.79 1.2 0 0 0 0 0 ]; 

POM = 0.0; % Particulate organic matter.
salinity = 35;
tau = 0.3; % sheer stress at the bottom (N m-2) - values between 0-1
rho = 1026; % density of the water (kg m-3)
  
Peak = 500; %Peak daily sunlight (W/m2)



inputs = cell(10,1); % This must be the right number, if more variables are added
inputs{1} = DIN_w;
inputs{2} = DIP_w;
inputs{3} = temp_timeseries;
inputs{4} = POM;
inputs{5} = salinity;
inputs{6} = tau;
inputs{7} = rho;
inputs{8} = Peak;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'MaxStep',60/86400);

[t,y] = ode45(@Bleach, [tspan], y0, opts, inputs);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Edmax = Peak;
Ed = Edmax*max(0,sin((t-0.25-floor(t-0.25))*2*pi));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 rctotal = sum(y(:,9:11),2);
 
 
 figure(1)
 subplot(231)
 pp = plot(t,y(:,[1]), 'LineWidth', 2);
 ylabel('Zooxanthallae Biomass (mg N m-2)','FontSize', 8);
 xlim([0 26])
 title ('24 h shade','FontSize',8)
 %text(0.04,0.99,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12)
 %xticklabels({})
 
 
 
 figure(1)
 subplot(232)
 pp = plot(t,y(:,[2:4 1]), 'LineWidth', 2); % y contains 15 variables 
 legend(pp, 'RN', 'RP', 'RC','CS_N');
 title('Symbiont biomass','FontSize', 16);
 ylabel('Reserves','FontSize', 16);
 xlabel('Time [days]', 'FontSize', 16);
 
 subplot(233)
 pp2 = plot(t, [y(:,2)./y(:,1) y(:,3)./y(:,1)/(((1/16)*(30.97/14.01))) y(:,4)./y(:,1)/(((106/16)*(12.01/14.01)))], 'LineWidth', 2);
 lgnd = legend(pp2, 'N [mg N]', 'P [dimensionless]', 'C [dimensionless]', 'Location','northeast');
 lgnd.FontSize = 6;
 xlim([0 26])
 ylim([0 1.2])
 ylabel('Normalised reserves','FontSize', 8);
 text(0.04,0.99,'D','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12)
 xticklabels({})
 
 
 
 subplot(234)
 pp = plot(t,y(:,6:8), 'LineWidth', 2); % y contains 15 variables 
 legend(pp, 'Chl', 'Xp','Xh');
 title('Pigments', 'FontSize', 16);
 ylabel('Pigments','FontSize', 16);
 xlabel('Time [days]', 'FontSize', 16);
 
 subplot(235)
 pp = plot(t,y(:,9)./rctotal,t,y(:,10)./rctotal,t,y(:,11)./rctotal, 'LineWidth', 2); % y contains 15 variables 
 ll = legend(pp, 'Qox', 'Qred', 'Qin');
 title('Reaction centres', 'FontSize', 16);
 ylabel('Fraction in each state','FontSize', 16);
 xlabel('Time [days]', 'FontSize', 16);
 
 subplot(236)
 					  
pp = plot(t,ROSpercell,t,expel,t,ROSthreshold*ones(size(t)),'LineWidth',2); % y contains 15 variables 
ylabel('ROS per cell, Expulsion','FontSize', 8);
xlabel('Time (d)', 'FontSize', 8);
lgd = legend(pp,'ROS','norm expulsion','Threshold', 'Location', 'northwest');
lgd.FontSize = 6;
set(gca,'ylim',[0 4*1.418e-14]);
set(gca,'xlim',[0 26]);
xlim = get(gca,'xlim')
ylim = get(gca,'ylim');
text(xlim(1) + 0.4*diff(xlim),ylim(1) + 0.53*diff(ylim),['Thres = ',num2str(ROSthreshold)]);
text(0.04,0.99,'G','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',12)

 			 
 					  
 subplot(236) % opaqueness of the cell
 pp = plot(t,aA/rCS/rCS/pi, t, nChl2C, t,Ed/Peak,'LineWidth', 2); % y contains 15 variables 
 title('Opaqueness and C:Chl', 'FontSize', 16);
 ylabel('[\alpha / PA],[norm. C:Chl,[norm. Ed]');
 xlabel('Time [days]', 'FontSize', 16);
 legend(pp,'opaq','norm Chl:C','norm Ed');
 xlim = get(gca,'xlim');ylim = get(gca,'ylim');
 text(xlim(1) + 0.2*diff(xlim),ylim(1) + 0.1*diff(ylim),['C:Chlmin = ',num2str(C2Chlmin)]);
 
 set(gcf,'Position',[720 500 1005 598]);
 eval(['print -dpng RB_3_12_',datestr(now,'dd-mmm-yyyy'),'.fig'])
 			  
