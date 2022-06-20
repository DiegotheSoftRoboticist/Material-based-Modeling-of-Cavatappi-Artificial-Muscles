%E1 Precursor Structure by Diego Higueras-Ruiz
 
clear all
close all
clc

% % Retrieve Data Sample 1
%     filename = 'E2 Sample 1.xlsx';
%     sheet = 'Sheet1';
%     E2Sample1 = xlsread(filename, sheet,'A:C');
%         timeE2_1 = E2Sample1(:,1);  
%         ForceE2_1 = E2Sample1(:,2);
%         GapE2_1 = E2Sample1(:,3);
%      
%         A = 4.5*4.5; %20.25mm
%         
%         StressE2_1 = -ForceE2_1/19; %MPa
%         StrainE2_1 = -(GapE2_1-GapE2_1(1)) /GapE2_1(1) %Strain
%         
%   % Retrieve Data Sample 2
%     filename = 'E2 Sample 2.xlsx';
%     sheet = 'Sheet1';
%     E2Sample2 = xlsread(filename, sheet,'A:C');
%         timeE2_2 = E2Sample2(:,1);  
%         ForceE2_2 = E2Sample2(:,2);
%         GapE2_2 = E2Sample2(:,3);
%      
%         A = 4.5*4.5; %20.25mm
%         
%         StressE2_2 = -ForceE2_2/19; %MPa
%         StrainE2_2 = -(GapE2_2-GapE2_2(1)) /GapE2_2(1) %Strain      
% 
% % Retrieve Data Sample 3
%     filename = 'E2 Sample 3.xlsx';
%     sheet = 'Sheet1';
%     E2Sample3 = xlsread(filename, sheet,'A:C');
%         timeE2_3 = E2Sample3(:,1);  
%         ForceE2_3 = E2Sample3(:,2);
%         GapE2_3 = E2Sample3(:,3);
%      
%         A = 4.5*4.5; %20.25mm
%         
%         StressE2_3 = -ForceE2_3/19; %MPa
%         StrainE2_3 = -(GapE2_3-GapE2_3(1)) /GapE2_3(1) %Strain

% Retrieve Data Sample 4
    filename = 'E2 Sample 4.xlsx';
    sheet = 'Sheet1';
    E2Sample4 = xlsread(filename, sheet,'A:C');
        timeE2_4 = E2Sample4(:,1);  
        ForceE2_4 = E2Sample4(:,2);
        GapE2_4 = E2Sample4(:,3);
     
        A = (4*4)*1.9; %20.25mm
        
        StressE2_4 = -ForceE2_4/A; %MPa
        StrainE2_4 = -(GapE2_4-GapE2_4(1)) /GapE2_4(1) %Strain
        
%   % Retrieve Data Sample 5
  filename = 'E2 Sample 5.xlsx';
    sheet = 'Sheet1';
    E2Sample5 = xlsread(filename, sheet,'A:C');
        timeE2_5 = E2Sample5(:,1);  
        ForceE2_5 = E2Sample5(:,2);
        GapE2_5 = E2Sample5(:,3);
     
        A = (5.5*5.5)*1.7; %20.25mm
        
        StressE2_5 = -ForceE2_5/A; %MPa
        StrainE2_5 = -(GapE2_5-GapE2_5(1)) /GapE2_5(1) %Strain
%            
%         
%    % Retrieve Data Sample 6
  filename = 'E2 Sample 6.xlsx';
    sheet = 'Sheet1';
    E2Sample6 = xlsread(filename, sheet,'A:C');
        timeE2_6 = E2Sample6(:,1);  
        ForceE2_6 = E2Sample6(:,2);
        GapE2_6 = E2Sample6(:,3);
     
        A = (5*5)*1.8; %20.25mm
        
        StressE2_6 = -ForceE2_6/A; %MPa
        StrainE2_6 = -(GapE2_6-GapE2_6(1)) /GapE2_6(1) %Strain    
%         
%     
%    % Retrieve Data Sample 4
%     filename = 'E1 Sample 4.xlsx';
%     sheet = 'Sheet1';
%     E1Sample4 = xlsread(filename, sheet,'A:C');
%         timeE1_4 = E1Sample4(:,1);  
%         ForceE1_4 = E1Sample4(:,2);
%         GapE1_4 = E1Sample4(:,3);
%         
%         StressE1_4 = -ForceE1_4/A; %MPa
%         StrainE1_4 = (GapE1_4-GapE1_4(1)) /GapE1_4(1) %Strain    
        
        
        
        
set(groot, 'DefaultTextInterpreter', 'tex', ...
           'DefaultAxesTickLabelInterpreter', 'tex', ...
           'DefaultAxesFontName', 'tex', ...
           'DefaultLegendInterpreter', 'tex', ...
           'defaultFigureColor','w');
        
        
     %E1 Sample 1 Temporal plot
     
fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);

yyaxis right;
H1 = plot(timeE2_4,StrainE2_4,':','LineWidth',1,'Color',[0.5,0.5,0.5]); %plot(time4003(1:(end)),-Strain400_filt3(1:(end)),':','Linewidth',1.5)
ylabel('Strain (mm/mm)')
set(gca,'ycolor',[0.5,0.5,0.5])
% ylim([0 60])
yyaxis left;
H2 =  plot(timeE2_4,-StressE2_4,':','Linewidth',1,'Color',[0, 0, 0]);% plot(time4003(1:(end)),p400_filt3(1:(end)),':','Linewidth',1.5)
ylabel('Stress (MPa)')
set(gca,'ycolor',[0, 0, 0])
xlabel('Time (s)')
% xlim([0,8]);
% set(gca,'XTick',[0:1:8]);
% title('400 grams')
% legend('SAMPLE 1','SAMPLE 1','Location','Northwest')
grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color

% print(gcf,'MaxStrain.png','-dpng','-r700'); 


% 
% %Stress vs. Strain
% 
fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);

%   plot(StrainE2_1(1004:1175)-StrainE2_1(1004), StressE2_1(1004:1175)-StressE2_1(1004),'Color',[0.4 0 0],'Linewidth',1)      
%   plot(StrainE2_2(1004:1175)-StrainE2_2(1004), StressE2_2(1004:1175)-StressE2_2(1004),'Color',[0.1 0.9 0],'Linewidth',1)      
%   plot(StrainE2_3(1004:1175)-StrainE2_3(1004), StressE2_3(1004:1175)-StressE2_3(1004),'Color',[0.5 0 1],'Linewidth',1)    
  plot(StrainE2_4(1195:1228)-StrainE2_4(1195), StressE2_4(1195:1228)-StressE2_4(1195),'Color',[0 0 0],'Linewidth',1)      
  plot(StrainE2_5(1195:1230)-StrainE2_5(1195), StressE2_5(1195:1230)-StressE2_5(1195),'Color',[0.1 0.5 0],'Linewidth',1)      
  plot(StrainE2_6(1195:1226)-StrainE2_6(1195), StressE2_6(1195:1226)-StressE2_6(1195),'Color',[0.5 0 0],'Linewidth',1)      

%              ylim([-0.5 1.5]);set(gca,'YTick',[0:0.5:1.5]); 
       ylabel('Stress (MPa)')
       xlabel('Strain (mm/mm)')
%         leg = legend('ɛ_{Axial} Nylon Monofilament',...
%      'ɛ_{Radial} Nylon Monofilament','Location','northwest') 
%  leg.ItemTokenSize = [5,20];
%    grid on 
% set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 
% % print(gcf,'Strains_DrawnandNONDrawn.png','-dpng','-r700');
% 
% 
%STATISTIC PLOT
StrainE2_4_Stat = StrainE2_4(1195:1228)-StrainE2_4(1195);
StressE2_4_Stat = StressE2_4(1195:1228)-StressE2_4(1195);
StrainE2_5_Stat = StrainE2_5(1195:1230)-StrainE2_5(1195);
StressE2_5_Stat = StressE2_5(1195:1230)-StressE2_5(1195);
StrainE2_6_Stat = StrainE2_6(1195:1226)-StrainE2_6(1195);
StressE2_6_Stat = StressE2_6(1195:1226)-StressE2_6(1195);

 StrainE2_Stat = [0:0.005:0.05];
% 
n = 1;
for i =1:length(StrainE2_Stat)
     if StrainE2_Stat(i) == 0
   StressE2_4_New(n) = 0;
   StressE2_5_New(n) = 0;
   StressE2_6_New(n) = 0;
     else
   StressE2_4_New(n) = interp1q(StrainE2_4_Stat,StressE2_4_Stat,StrainE2_Stat(i));
   StressE2_5_New(n) = interp1q(StrainE2_5_Stat,StressE2_5_Stat,StrainE2_Stat(i));
   StressE2_6_New(n) = interp1q(StrainE2_6_Stat,StressE2_6_Stat,StrainE2_Stat(i));
 
     end
 n = n+1;
end
 



n =1;
for i = 1:length(StressE2_5_New)
    Stress_calulate = [StressE2_4_New(i),StressE2_5_New(i),StressE2_6_New(i)];
   Stress_mean(n) = mean(Stress_calulate) % Mean Of All Experiments At Each Value Of ‘x’
   Stress_mean(i) = Stress_mean(n);
   SEM_Stress(n) = std(Stress_calulate)/sqrt(length(Stress_calulate)); % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
   n = n+1;
end
CI95 = tinv([0.025 0.975], length(Stress_calulate)-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, SEM_Stress, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

p = polyfit(StrainE2_Stat(1:11),Stress_mean(1:11),1);
fit = polyval(p,StrainE2_Stat(1:11));

p_manual = [-9,0];
fit_manual = polyval(p_manual,StrainE2_Stat(1:11));

close all
 
fig=figure('units','inch','position',[0,0,3.6,3.2]); hold on; grid on; set(gca,'FontSize',10);
 
plot(StrainE2_Stat(1:11),Stress_mean(1:11),'LineWidth',1.5,'Color',[0.75,0.75,0.75]);

% plot(StrainE2_Stat(1:11),fit(1:11),'LineWidth',0.4,'Color',[0.8,0,0]);

% plot(StrainE2_Stat(1:11),-4.5*StrainE2_Stat(1:11),'LineWidth',0.4,'Color',[0.8,0,0]);

 a = errorbar(StrainE2_Stat,Stress_mean,yCI95(2,:),'s','MarkerSize',4,...
    'MarkerEdgeColor',[0.75,0.75,0.75],...
    'MarkerFaceColor',[0, 0.75, 0.75],'LineWidth',0.5,'Color',[0, 0.75, 0.75])

plot(StrainE2_Stat(1:11),fit_manual(1:11),'LineWidth',0.85,'Color',[0.8,0,0]);

   ylim([-0.85 4]); %set(gca,'YTick',[0:0.5:2.1]);
% % 
 xlim([0 0.052]); %set(gca,'xTick',[0:50:250]);


%  legend([a b c d e],'m = 0.1 kg','m = 0.2 kg','m = 0.3 kg','m = 0.4 kg','m = 0.5 kg')


ylabel('Stress, \sigma (MPa)')
xlabel('Strain, ɛ')
 
grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color
print(gcf,'E1&E2.png','-dpng','-r700');   


% fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);
% 
%     
%   plot(StrainE2_4(1135:1225+240)-StrainE2_4(1135), StressE2_4(1135:1225+240)-StressE2_4(1135),'Color',[0 0 0],'Linewidth',1)      
%      
% 
% %              ylim([-0.5 1.5]);set(gca,'YTick',[0:0.5:1.5]); 
%        ylabel('Stress (MPa)')
%        xlabel('Strain (mm/mm)')

