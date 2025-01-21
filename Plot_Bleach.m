figure(1)
 
subplot(231) %Symbiont reserves of N, P, C and zooxanthallae biomass 
plot(t,y(:,[2:4 1]), 'LineWidth', 2); % y contains 15 variables 
ylabel('Symbiont reserves','FontSize', 12);
xlabel('Time [d]', 'FontSize', 12);
ylim = get(gca,'ylim');
set(gca,'xlim',tspan); %change this to timeframe
text(0.04,0.99,'A','Units', 'normalized', 'VerticalAlignment', 'Top','FontSize',12)
title('Symbiont biomass','FontSize', 14);
legend('RN (mg N m^-^2)', 'RP (mg P m^-^2)', 'RC (mg C m^-^2)','CS_N (mg N m^-^2)','Position',[0.17 0.83 0.03 0.06]);


subplot(232) %Normalised reserves of N, P, C
plot(t, [y(:,2)./y(:,1) y(:,3)./y(:,1)/(((1/16)*(30.97/14.01))) y(:,4)./y(:,1)/(((106/16)*(12.01/14.01)))], 'LineWidth', 2);
ylabel('Normalised internal reserves','FontSize', 12);
xlabel('Time (d)', 'FontSize', 12);
set(gca,'ylim',[0 1.4]);
set(gca,'xlim',tspan); %change this to timeframe
text(0.04,0.99,'B','Units', 'normalized', 'VerticalAlignment', 'Top','FontSize',12)
title('Normalised Reserves','FontSize', 14);
legend('N [dimensionless]', 'P [dimensionless]', 'C [dimensionless]', 'Position',[0.48 0.85 0.03 0.06]);


subplot(233) %Pigments
plot(t,y(:,6:8), 'LineWidth', 2); % y contains 15 variables 
ylabel('Pigment concentration','FontSize', 12);
xlabel('Time (d)', 'FontSize', 12);
set(gca,'ylim',[0 0.08]);
set(gca,'xlim',tspan); %change this to timeframe
text(0.04,0.99,'C','Units', 'normalized', 'VerticalAlignment', 'Top','FontSize',12)
title('Pigments', 'FontSize', 14);
legend('Chla (mg m^-^2)', 'Xp (mg m^-^2)','Xh (mg m^-^2)','Position',[0.73 0.85 0.03 0.06]);

 
subplot(234) %Reaction centre fraction state
plot(t,y(:,9)./rctotal,t,y(:,10)./rctotal,t,y(:,11)./rctotal, 'LineWidth', 2); % y contains 15 variables 
ylabel('Normalised reaction centre state','FontSize', 12);
xlabel('Time (d)', 'FontSize', 12);
ylim = get(gca,'ylim');
set(gca,'xlim',tspan); %change this to timeframe
text(0.04,0.99,'D','Units', 'normalized', 'VerticalAlignment', 'Top','FontSize',12)
title('Reaction centres', 'FontSize', 14);
legend('Qox (dimensionless)', 'Qred (dimensionless)', 'Qin (dimensionless)','Position',[0.17 0.375 0.03 0.06]);


subplot(235) %ROS per cell build-up, ROS expulsion
plot(t,ROSpercell,t,expel,t,ROSthreshold*ones(size(t)),'LineWidth',2); % y contains 15 variables 
ylabel('ROS concentration, Expulsion','FontSize', 12);
xlabel('Time (d)', 'FontSize', 12);
set(gca,'ylim',[0 4*1.418e-14]);
set(gca,'xlim',tspan);
ylim = get(gca,'ylim');
%set(gca,'xlim',[0 26]); %change this to timeframe
%text(xlim(1) + 0.1*diff(xlim),ylim(1) + 0.33*diff(ylim),['Thres = ',num2str(ROSthreshold)]);
text(0.04,0.99,'E','Units', 'normalized', 'VerticalAlignment', 'Top','FontSize',12)
title('ROS, Expulsion', 'FontSize', 14);
legend('ROS concentration (mg O cell^-^1)','Symbiont cell expulsion rate (d-1)','Threshold', 'Position',[0.47 0.375 0.03 0.06]);

 			 
subplot(236) % opaqueness of the cell
plot(t,aA/rCS/rCS/pi, t, nChl2C, Light(1,:),Light(2,:)/Peak,'LineWidth', 2); % y contains 15 variables 
ylabel('[\alpha / PA],[norm. C:Chl,[norm. Ed]','FontSize',12);
xlabel('Time (d)', 'FontSize', 12);
ylim = get(gca,'ylim');
set(gca,'xlim',tspan); %change this to timeframe
%text(xlim(1) + 0.3*diff(xlim),ylim(2) + 0.1*diff(ylim),['C:Chlmin = ',num2str(C2Chlmin)]);
text(0.04,0.99,'F','Units', 'normalized', 'VerticalAlignment', 'Top','FontSize',12)
title('Opaqueness and C:Chl', 'FontSize', 14);
legend('opaq','norm Chl:C','norm Ed','Position',[0.74 0.375 0.03 0.06]);


set(gcf,'Position',[720 500 1005 598]);
% insert name of scipt below (%%)
% eval(['print -dpng %%',datestr(now,'dd-mmm-yyyy'),'.fig'])
 			  
