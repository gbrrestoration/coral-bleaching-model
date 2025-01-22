function dydt = bleach(t, y, forcings, params)

% Developed by Mark Baird and Julie Terp Jorgensen (2021).

% Water column environmental conditions:
  
% Extract 'forcings' variables created in Run_Bleach.m:

DIN_w = forcings{1};
DIP_w = forcings{2};
Temperature = forcings{3};
POM = forcings{4};
Salinity = forcings{5};
tau = forcings{6};
rho = forcings{7};
Light = forcings{8};

% If timeseries data is available, interpolate across timeseries:

Fnames = {'DIN_w', 'DIP_w','Light','Temperature'};
for f = 1:length(Fnames)
    if length(eval(Fnames{f})) > 1
        dat = eval(Fnames{f});
        eval(sprintf('%s = %g;',Fnames{f},interp1(dat(1,:),dat(2,:),t)));
    end
end
clear Fnames dat f

% Rename key variables:

Ed = Light;
Tanomaly = Temperature;

% Extract parameters from Run_Bleach.m:

mN = params(1);  			% Nitrogen content of zooxanthellae cells (mg N cell-1)
rCS = params(2);                	% Radius of zooxanthellae cells (m)
CSvol = params(3);    			% Volume of the cell (m3)
uCH = params(4);                	% Maximum growth rate of coral host (d-1)
Spart = params(5);              	% Rate coefficient of particle capture (m d-1)
uCS = params(6);               		% Maximum growth rate of zooxanthellae (d-1)
zetaCH = params(7);             	% Quadratic mortality coefficient of polyps (d-1 (g N m-2)-1)
zetaCS = params(8);             	% Linear mortality of zooxanthellae (d-1) 
CHremin = params(9);            	% Remineralised fraction of coral mortality 
OCH = params(10);               	% Nitrogen-specific host area coefficient of polyps (m2 mg N-1)
C2Chlmin = params(11);          	% max Chl:C ratio at which synthesis stops    

Chlmax = params(12);
Xanth_tau = params(13); 		% 1 per 20 minutes %%72 times per day
photon2rcii = params(14);
ROSthreshold = params(15);
photon2ros = params(16);
chla2rcii = params(17);
CSmaxbleachrate = params(18);
repair_coefficient = params(19);

a_rub = min(1.0,max(0,(1-exp(-(2-Tanomaly)))/(1-exp(-2))));

% Parameters required for conversion to spectrally-resolved:

% c = 2.998e8;                  	% Speed of light (m s-1)
% h = 6.626e-34;                	% Planck constant (J s-1)
% AV = 6.02e23;                 	% Avagadro constant (mol-1)

dydt(1:length(y)) = 0; 			% This sets every derivative to 0 so you can progressively add to them.

% State variables (Table A.1 in Baird et al. (2018)):

CS_N = y(1);				% Zooxanthellae biomass                         CS     (mg N m-2)
RN   = y(2); 				% Reserves of nitrogen                          RN     (mg N m-2)
RP   = y(3);  				% Reserves of phosphorus                        RP     (mg P m-2)
RC   = y(4);  				% Reserves of carbon                            RC     (mg C m-2)
CH_N   = y(5);  			% Coral biomass                                 CH     (mg N m-2)
Chl  = y(6); 				% Zooxanthellae chlorophyll a concentration     Chl    (mg m-2)
Xp = y(7); 				% Zooxanthellae diadinoxanthin concentration    Xp     (mg m-2)
Xh = y(8); 				% Zooxanthellae diatoxanthin concentration      Xh     (mg m-2)
Qox = y(9);				% Oxidised reaction centre concentration        Qox    (mg m-2)
Qred = y(10);				% Reduced reaction centre concentration         Qred   (mg m-2)
Qin = y(11); 				% Inhibited reaction centre concentration       Qin    (mg m-2)
ROS = y(12);				% Reactive oxygen species concentration         [ROS]  (mg m-2]

Qtot = Qox + Qin + Qred;

% Constant values (Table A.7 in Baird et al. (2018)):

D_N = 17.5e-10;                		% Molecular diffusivity of NO3    (m2 s-1)
D_P = 8.46e-10;                		% Molecular diffusivity of PO4    (m2 s-1)
v = 1.05e-6;                   		% Kinematic viscosity of water    (m2 s-1)

% Conversion of W m-2 to mmol photon d-1 for light at 490 nm. 
  
W2mmol_490 = 8.359335857479461e-9*490*1000*86400; %  = (1e9 h c)-1 / Av

% 8.359335857479461e-9 == (1E+9 * h * c)-1 /Av

% phi = 0.1                    		% Fractional (of uCS) respiration rate

% Schmidt number  = diffusivity momentum / diffusivity of nutrient ions

Sc_N = v / D_N; 	       		% Ratio between the kinematic viscosity of water and the molecular diffusivity of NO3
Sc_P = v / D_P; 	      		% Ratio between the kinematic viscosity of water and the molecular diffusivity of PO4

% Mass transfer rate coefficient (height of water cleared of nutrient per unit of time by the water-coral exchange): 

Sx_N = 2850 * (2 * tau / rho)^0.38 * Sc_N^-0.6;  % (m d-1)
Sx_P = 2850 * (2 * tau / rho)^0.38 * Sc_P^-0.6; % (m d-1)

Aeff = 1 - exp(-OCH * CH_N);   		% Aeff is the effective projected area fraction of the coral community (m2) 

% Translocation (represents the one-way consumption of zooxanthellae organic matter produced through either 
% zooxanthellae growth or mortality):

% Derived calculations: 

cells = CS_N / mN;
  
Nquota = RN / CS_N; 			% Normalised reserves
kN_mass = DIN_w * Sx_N * Aeff;
Nuptake = kN_mass * (1.0 - Nquota) ; 	% Dissolved inorganic nitrogen (DIN)(mg N m-2 s-1) 

Pquota = RP / CS_N / ((1/16) * (30.97/14.01));
kP_mass = DIP_w * Sx_P * Aeff;
Puptake = kP_mass * (1.0 - Pquota); 	% Dissolved inorganic phosphorus (DIP)(mg P m-2 s-1)

% Pigment synthesis - do calculations at 490 nm assuming it is PAR.

cellXp = Xp/(cells*CSvol);    		% mg pig m-3
cellXh_= Xh/(cells*CSvol);
cellChl = Chl/(cells*CSvol);

yC =  cellXp * 0.06 + cellChl * 0.06; % cellXp * bio->yC_diadinoxanthin[w] + cellChl * bio->yC_symbiodinium[w];

% Calculate the absorption cross-section (value between 0 and pi r^2):

temp = 2.0 * yC * rCS;
aA = 0.0;
if (temp > 3e-4)
  aA = pi * rCS * rCS * (1.0 - 2.0 * (1.0 - (1.0 + temp) * exp(-temp)) / (temp * temp));
end				  

% Calculate fraction benefit of increasing pigment concentration [number between 0 and 4/3]:
  
temp = yC*rCS;

Chlsynfactor = (1.0-exp(-2.0*temp)*(2.0*temp*temp+2.0*temp+1.0))/(temp*temp*temp);

Iquota = RC / CS_N / ((106/16) * (12.01/14.01));
									 
dChldt_syn = Chlmax * uCS * min(1.33,max(0.0,Chlsynfactor * (1.0 - Iquota) * (1.0 - (Qin/Qtot))));

if (CS_N * 5.6786 < C2Chlmin * Chl ) 	% Don't synthesise if can't fit in anymore. 
  dChldt_syn = 0.0;
end
if Ed == 0
  dChldt_syn = 0.0;
end
dydt(6) = dydt(6) + dChldt_syn * CSvol * cells;
dydt(7) = dydt(7) +  0.2448 * dChldt_syn * CSvol * cells ;

% aA - m2 cell-1; cells ind./m2 
								 
Iuptake = Ed * aA * Qox/Qtot * a_rub * W2mmol_490 * (1.0 - Iquota) * cells;  	% mmol photons cell-1 d-1

symb_grow_IN = uCS * Nquota * Pquota * Iquota;

CStoCHfrac = (cells * pi * rCS * rCS) / (2 * CH_N * OCH);

% Eq. 7 in Baird et al. (2018) where growth rate of zooxanthellae is given by the maximum growth rate, uCS, 
% multiplied by the normalised reserves, Nquota, Pquota and Iquota. Iquota is the energy/carbon reserve. 

% Calculations used in the derivatives:

% Host feeding

if Ed > 0
 Spart = 0.0;
end
    
pot_grazing = Spart * Aeff* POM; 		% Change units from mg to g.

total_grazing = min([uCH*CH_N, pot_grazing]); 	% g N m-2 s-1

polypmort = zetaCH * CH_N;

symb_grow = min([uCS,symb_grow_IN + CHremin * polypmort * CH_N / CS_N]);

mucus = max([0.0,symb_grow_IN * CS_N + CHremin * polypmort * CH_N - uCS * CS_N]);

mucus_attNR = max([0.0,mucus - CHremin * polypmort * CH_N]);

% Translocate (+ve symbiont to host), including growth translocates mortality of symbionts and 
% mortality of host (-ve).

translocate = symb_grow * CS_N * CStoCHfrac + zetaCS * CS_N;

translocate_attNR = symb_grow * CS_N * CStoCHfrac;

host_growth = total_grazing + translocate;


% Mucus production in g N: 

mucus2 = max(0.0,host_growth - uCH*CH_N);
host_growth = host_growth - mucus2;
mucus = mucus + mucus2;
 
% Symbiont and host growth:  

dydt(1) = symb_grow * CS_N - translocate - polypmort * CS_N;
dydt(5) = host_growth - polypmort * CH_N + total_grazing;

% Changes in internal reserves:

dydt(2) = dydt(2) + Nuptake - (symb_grow_IN * CS_N - translocate_attNR * Nquota); 
dydt(3) = dydt(3) + Puptake - (symb_grow_IN * CS_N + translocate_attNR * Pquota) * ((1/16) * (30.97/14.01)) ; 		% Redfield ratio 106C:16N:1P
dydt(4) = dydt(4) + Iuptake - (symb_grow_IN * CS_N +  translocate_attNR * Iquota)* ((1060/16) * (12.01/14.01)); 	% Redfield ratio 1060 photon:16N:1P
% The units for dydt(4) have been corrected since this in in mmol photon cell-1 d-1

% Xanthophyll switching:

xans = Xp+Xh;

metric = Qin/Qtot;
parabole = 1.0;
if (metric > 0.5)
  if (Xh >= Xp)
    parabole = 1.0-4.0*(Xp/xans-0.5)^2;
  end
else									 
  if (Xp >= Xh)
    parabole = 1.0-4.0*(Xh/xans-0.5)^2;     
  end
end

% Factor of 8 = 2^3 so that metric-0.5 ranges between -1 and 1.

swwitch = 8.0*(metric-0.5)^3 * Xanth_tau * parabole * xans;									 
      
dydt(7) = dydt(7) - swwitch;  % Xp
dydt(8) = dydt(8) + swwitch;  % Xh

% Symbiont mortality:

dydt(1) = dydt(1) - zetaCS * y(1);
dydt(2) = dydt(2) - zetaCS * y(2);
dydt(3) = dydt(3) - zetaCS * y(3);
dydt(4) = dydt(4) - zetaCS * y(4);
dydt(6) = dydt(6) - zetaCS * y(6);
dydt(7) = dydt(7) - zetaCS * y(7);
dydt(8) = dydt(8) - zetaCS * y(8);
dydt(9) = dydt(9) - zetaCS * y(9);
dydt(10) = dydt(10) - zetaCS * y(10);
dydt(11) = dydt(11) - zetaCS * y(11);
dydt(12) = dydt(12) - zetaCS * y(12);

% Reaction centre dynamics:
							 
Qi2Qox = repair_coefficient * photon2rcii * Qin * 86400; % This repairs 10 mol ph m-2 d-1 hitting 1 mmol of RCII.

% Now need to expel:

ROSpercell = y(12)/cells;

% Rate of ROS detoxification:

% Detoxification rate based on ROS above ROSthreshold*0.5
% Detoxification = 0 when ROS is below ROSthreshold*0.5
% Detoxification brings ROS back to ROSthreshold*0.5 when ROS is above ROSthreshold*0.5

ARO = 0;

if ROSpercell <= ROSthreshold/2
    ARO=0;
elseif ROSpercell > ROSthreshold/2 && ROSpercell <= ROSthreshold
    ARO = uCS * Nquota * Iquota * Pquota * (ROSpercell - ROSthreshold/2)*cells;
elseif ROSpercell > ROSthreshold
    ARO = uCS * Nquota * Iquota * Pquota * max((ROSthreshold - ROSthreshold/2)*cells,0);
end

absorb = Ed * aA * W2mmol_490 * cells * photon2rcii; % Units now mmol rcii m-2 d-1 
C_fix = absorb * (Qox/Qtot) * a_rub * (1.0 - Iquota);
C_notfixed = absorb * (Qox/Qtot) - C_fix;									 

dydt(9) = dydt(9) + dChldt_syn * CSvol * cells * chla2rcii;   					% Ox								 
dydt(9) = dydt(9) - C_notfixed + Qi2Qox;                      					% Ox
dydt(10) = dydt(10) + C_notfixed - absorb * (Qred / Qtot);    					% Red
dydt(11) = dydt(11) - Qi2Qox + absorb * (Qred / Qtot);        					% In
dydt(12) = dydt(12) - ARO + 0.5 * (Qin/Qtot) * absorb / photon2rcii / photon2ros / 86400;   	% 50% diffusion of intracellular ROS concentration to extracellular ROS concentration
						 
expulsionrate = CSmaxbleachrate*min(1.0,max(0.0,(ROSpercell - ROSthreshold)/ROSthreshold));

dydt(1) = dydt(1) - expulsionrate * y(1);
dydt(2) = dydt(2) - expulsionrate * y(2);
dydt(3) = dydt(3) - expulsionrate * y(3);
dydt(4) = dydt(4) - expulsionrate * y(4);
dydt(6) = dydt(6) - expulsionrate * y(6);
dydt(7) = dydt(7) - expulsionrate * y(7);
dydt(8) = dydt(8) - expulsionrate * y(8);
dydt(9) = dydt(9) - expulsionrate * y(9);
dydt(10) = dydt(10) - expulsionrate * y(10);
dydt(11) = dydt(11) - expulsionrate * y(11);
dydt(12) = dydt(12) - expulsionrate * y(12);
									 
dydt = dydt';
% Need to turn off derivatives if CS_N is less than 1000 times the ODE tolerance.

% if (y(1) < 1e-5)
% dydt = dydt*0;
% end







