figure(1)

subplot(2,4,1)
%subplot(231) %Symbiont reserves of N, P, C and zooxanthallae biomass 
plot(t,[y(:,[2:4 1])], 'LineWidth', 2); % y contains 15 variables 
ylabel('Reserves','FontSize', 12,'FontName','Times');
xlabel('Time (d)', 'FontSize', 12,'FontName','Times');
ylim = get(gca,'ylim');
set(gca,'xlim',tspan); %change this to timeframe
text(0.04,0.99,'A','Units', 'normalized', 'VerticalAlignment', 'Top','FontSize',12,'FontName','Times')
title('Symbiont biomass','FontSize', 14,'FontName','Times');
legend('RN', 'RP', 'RC','CS_N','Position',[0.2252 0.8099 0.056 0.1046],'FontName','Times');

subplot(2,4,2)
%subplot(232) %Normalised projected area and RuBisCo activity
plot(t,[normalisedprojectedarea,a_rub],'LineWidth',2)
ylabel('Normalised','FontSize', 12,'FontName','Times');
xlabel('Time (d)', 'FontSize', 12,'FontName','Times');
set(gca,'ylim',[0 7]);
set(gca,'xlim',tspan); %change this to timeframe
text(0.04,0.99,'B','Units', 'normalized', 'VerticalAlignment', 'Top','FontSize',12,'FontName','Times')
title('Projected Area, RuBisCo Activity','FontSize', 14,'FontName','Times');
legend('Area [dimensionless]', 'a rub [dimensionless]', 'Position',[0.3428 0.693 0.1169 0.06],'FontName','Times');

subplot(2,4,3)
%subplot(233) %Normalised reserves of N, P, C
plot(t, [y(:,2)./y(:,1) y(:,3)./y(:,1)/(((1/16)*(30.97/14.01))) y(:,4)./y(:,1)/(((106/16)*(12.01/14.01)))], 'LineWidth', 2);
ylabel('Normalised reserves','FontSize', 12,'FontName','Times');
xlabel('Time (d)', 'FontSize', 12,'FontName','Times');
set(gca,'ylim',[0 7]);
set(gca,'xlim',tspan); %change this to timeframe
text(0.04,0.99,'C','Units', 'normalized', 'VerticalAlignment', 'Top','FontSize',12,'FontName','Times')
title('Normalised Reserves','FontSize', 14,'FontName','Times');
legend('N [mg N]', 'P [dimensionless]', 'C [dimensionless]', 'Position',[0.5902 0.8397 0.1044 0.0741],'FontName','Times');

%subplot(234) %Pigments
subplot(2,4,4)
plot(t,[y(:,6) y(:,7)*2 y(:,8)*2], 'LineWidth', 2); % y contains 15 variables 
ylabel('Pigments','FontSize', 12,'FontName','Times');
xlabel('Time (d)', 'FontSize', 12,'FontName','Times');
set(gca,'ylim',[0 0.4]);
set(gca,'xlim',tspan); %change this to timeframe
text(0.04,0.99,'D','Units', 'normalized', 'VerticalAlignment', 'Top','FontSize',12,'FontName','Times')
title('Pigments', 'FontSize', 14,'FontName','Times');
legend('Chl', 'Xp*2','Xh*2','Position',[0.8402 0.8408 0.0596 0.0741],'FontName','Times');

%subplot(235) %Reaction centre fraction state
subplot(2,4,5)
plot(t,(y(:,9)./rctotal),t,y(:,10)./rctotal,t,y(:,11)./rctotal, 'LineWidth', 2); % y contains 15 variables 
ylabel('Fraction in each state','FontSize', 12,'FontName','Times');
xlabel('Time (d)', 'FontSize', 12,'FontName','Times');
ylim = get(gca,'ylim');
set(gca,'xlim',tspan); %change this to timeframe
text(0.04,0.99,'E','Units', 'normalized', 'VerticalAlignment', 'Top','FontSize',12,'FontName','Times')
title('Reaction centres', 'FontSize', 14,'FontName','Times');
legend('Qox', 'Qred', 'Qin','Position',[0.0612 0.414 0.0581 0.0741],'FontName','Times');

%plot(t,1-(y(:,9)./rctotal))


%subplot(236) %ROS per cell build-up, ROS expulsion
subplot(2,4,6)
plot(t,ROSpercell,t,expel,t,ROSthreshold*ones(size(t)),'LineWidth',2); % y contains 15 variables 
ylabel('ROS per cell, Expulsion','FontSize', 12,'FontName','Times');
xlabel('Time (d)', 'FontSize', 12,'FontName','Times');
set(gca,'ylim',[0 4*1.418e-14]);
set(gca,'xlim',tspan);
ylim = get(gca,'ylim');
%set(gca,'xlim',[0 26]); %change this to timeframe
text(0.0265, 0.3464, ['Thres = ',num2str(ROSthreshold)],'Units', 'normalized', 'VerticalAlignment', 'Top', 'FontName','Times')
text(0.04,0.99,'F','Units', 'normalized', 'VerticalAlignment', 'Top','FontSize',12,'FontName','Times')
title('ROS', 'FontSize', 14,'FontName','Times');
legend('ROS','Expulsion [norm]','Threshold','Position',[0.3919 0.3669 0.0956 0.0741],'FontName','Times');

%subplot(237) % opaqueness of the cell
subplot(2,4,7)
plot(t,aA/rCS/rCS/pi, t, nChl2C, Light(1,:),Light(2,:)/Peak,'LineWidth', 2); % y contains 15 variables 
ylabel('[\alpha / PA],[norm. C:Chl,[norm. Ed]','FontSize',12,'FontName','Times');
xlabel('Time (d)', 'FontSize', 12,'FontName','Times');
ylim = get(gca,'ylim');
set(gca,'xlim',tspan); %change this to timeframe
text(0.5766, 0.9861, ['C:Chlmin = ',num2str(C2Chlmin)],'Units', 'normalized', 'VerticalAlignment', 'Top', 'FontName','Times')
text(0.04,0.99,'G','Units', 'normalized', 'VerticalAlignment', 'Top','FontSize',12,'FontName','Times')
title('Opaqueness and C:Chl', 'FontSize', 14,'FontName','Times');
legend('opaq','norm Chl:C','norm Ed','Position',[0.7033 0.3773 0.0826 0.0741],'FontName','Times');


%set(gcf,'Position',[106.3333   44.3333  996.0000  586.6667]);
% insert name of scipt below (%%)
eval(['print -dpng Run_Bleach_7_',datestr(now,'dd-mmm-yyyy'),'.fig'])
		  
