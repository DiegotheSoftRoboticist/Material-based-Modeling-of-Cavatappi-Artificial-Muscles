
%E1 Precursor Structure by Diego Higueras-Ruiz
 
clear all
close all
clc

% Retrieve Data Sample 1
    filename = 'E1 Sample 1.xlsx';
    sheet = 'Sheet1';
    E1Sample1 = xlsread(filename, sheet,'A:C');
        timeE1_1 = E1Sample1(:,1);  
        ForceE1_1 = E1Sample1(:,2);
        GapE1_1 = E1Sample1(:,3);
        
        
        A = 1.12;
%         do = (1.8); %mm
%         di = (0.80); %mm
%         A = (do^2-di^2)*pi/4
        
        StressE1_1 = -ForceE1_1/A; %MPa
        StrainE1_1 = (GapE1_1-GapE1_1(1)) /GapE1_1(1) %Strain
        
  % Retrieve Data Sample 2
    filename = 'E1 Sample 2.xlsx';
    sheet = 'Sheet1';
    E1Sample2 = xlsread(filename, sheet,'A:C');
        timeE1_2 = E1Sample2(:,1);  
        ForceE1_2 = E1Sample2(:,2);
        GapE1_2 = E1Sample2(:,3);
        
        StressE1_2 = -ForceE1_2/A; %MPa
        StrainE1_2 = (GapE1_2-GapE1_2(1)) /GapE1_2(1); %Strain
           
        
   % Retrieve Data Sample 3
    filename = 'E1 Sample 3.xlsx';
    sheet = 'Sheet1';
    E1Sample3 = xlsread(filename, sheet,'A:C');
        timeE1_3 = E1Sample3(:,1);  
        ForceE1_3 = E1Sample3(:,2);
        GapE1_3 = E1Sample3(:,3);
        
        StressE1_3 = -ForceE1_3/A; %MPa
        StrainE1_3 = (GapE1_3-GapE1_3(1)) /GapE1_3(1) %Strain     
        
    
   % Retrieve Data Sample 4
    filename = 'E1 Sample 4.xlsx';
    sheet = 'Sheet1';
    E1Sample4 = xlsread(filename, sheet,'A:C');
        timeE1_4 = E1Sample4(:,1);  
        ForceE1_4 = E1Sample4(:,2);
        GapE1_4 = E1Sample4(:,3);
        
        StressE1_4 = -ForceE1_4/A; %MPa
        StrainE1_4 = (GapE1_4-GapE1_4(1)) /GapE1_4(1); %Strain    
        
        
        
        
set(groot, 'DefaultTextInterpreter', 'tex', ...
           'DefaultAxesTickLabelInterpreter', 'tex', ...
           'DefaultAxesFontName', 'tex', ...
           'DefaultLegendInterpreter', 'tex', ...
           'defaultFigureColor','w');
        
        
     %E1 Sample 1 Temporal plot
     
fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);

yyaxis right;
H1 = plot(timeE1_1,StrainE1_1,':','LineWidth',1,'Color',[0.5,0.5,0.5]); %plot(time4003(1:(end)),-Strain400_filt3(1:(end)),':','Linewidth',1.5)
ylabel('Strain (mm/mm)')
set(gca,'ycolor',[0.5,0.5,0.5])
% ylim([0 60])
yyaxis left;
H2 =  plot(timeE1_1,StressE1_1,'Linewidth',1,'Color',[0, 0, 0]);% plot(time4003(1:(end)),p400_filt3(1:(end)),':','Linewidth',1.5)
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

     %E1 Sample 2 Temporal plot

fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);

yyaxis right;
H1 = plot(timeE1_2,StrainE1_2,':','LineWidth',1,'Color',[0.5,0.5,0.5]); %plot(time4003(1:(end)),-Strain400_filt3(1:(end)),':','Linewidth',1.5)
ylabel('Strain (mm/mm)')
set(gca,'ycolor',[0.5,0.5,0.5])
% ylim([0 60])
yyaxis left;
H2 =  plot(timeE1_2,StressE1_2,'Linewidth',1,'Color',[0, 0, 0]);% plot(time4003(1:(end)),p400_filt3(1:(end)),':','Linewidth',1.5)
ylabel('Stress (MPa)')
set(gca,'ycolor',[0, 0, 0])
xlabel('Time (s)')
% xlim([0,8]);
% set(gca,'XTick',[0:1:8]);
% title('400 grams')
% legend('SAMPLE 1','SAMPLE 1','Location','Northwest')
grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color

%E1 Sample 3 Temporal plot

fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);

yyaxis right;
H1 = plot(timeE1_3,StrainE1_3,':','LineWidth',1,'Color',[0.5,0.5,0.5]); %plot(time4003(1:(end)),-Strain400_filt3(1:(end)),':','Linewidth',1.5)
ylabel('Strain (mm/mm)')
set(gca,'ycolor',[0.5,0.5,0.5])
% ylim([0 60])
yyaxis left;
H2 =  plot(timeE1_3,StressE1_3,'Linewidth',1,'Color',[0, 0, 0]);% plot(time4003(1:(end)),p400_filt3(1:(end)),':','Linewidth',1.5)
ylabel('Stress (MPa)')
set(gca,'ycolor',[0, 0, 0])
xlabel('Time (s)')
% xlim([0,8]);
% set(gca,'XTick',[0:1:8]);
% title('400 grams')
% legend('SAMPLE 1','SAMPLE 1','Location','Northwest')
grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color

%Stress vs. Strain

fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);

%   plot(StrainE1_1(975:1099), StressE1_1(975:1099),'Color',[0 0 0],'Linewidth',1)      
  plot(StrainE1_2(971:1098)-0.005, StressE1_2(971:1098),'Color',[0 0 0],'Linewidth',1)      
  plot(StrainE1_3(982:1098)-0.005, StressE1_3(982:1098),'Color',[0.1 0.5 0],'Linewidth',1)      
  plot(StrainE1_4(991:1098)-0.007, StressE1_4(991:1098),'Color',[0.5 0 0],'Linewidth',1)      

%              ylim([-0.5 1.5]);set(gca,'YTick',[0:0.5:1.5]); 
       ylabel('Stress (MPa)')
       xlabel('Strain (mm/mm)')
%         leg = legend('ɛ_{Axial} Nylon Monofilament',...
%      'ɛ_{Radial} Nylon Monofilament','Location','northwest') 
%  leg.ItemTokenSize = [5,20];
%    grid on 
% set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 
% print(gcf,'Strains_DrawnandNONDrawn.png','-dpng','-r700');


%STATISTIC PLOT
StrainE1_2_Stat = StrainE1_2(971:1098)-0.005;
StressE1_2_Stat = StressE1_2(971:1098);
StrainE1_3_Stat = StrainE1_3(982:1098)-0.005;
StressE1_3_Stat = StressE1_3(982:1098);
StrainE1_4_Stat = StrainE1_4(991:1098)-0.007;
StressE1_4_Stat = StressE1_4(991:1098);

StrainE1_Stat = [0:0.005:0.05]

n = 1;
for i =1:length(StrainE1_Stat)
     if StrainE1_Stat(i) == 0
   StressE1_2_New(n) = 0;
   StressE1_3_New(n) = 0;
   StressE1_4_New(n) = 0;
     else
   StressE1_2_New(n) = interp1q(StrainE1_2_Stat,StressE1_2_Stat,StrainE1_Stat(i));
   StressE1_3_New(n) = interp1q(StrainE1_3_Stat,StressE1_3_Stat,StrainE1_Stat(i));
   StressE1_4_New(n) = interp1q(StrainE1_4_Stat,StressE1_4_Stat,StrainE1_Stat(i));
 
     end
 n = n+1;
 end
figure 
plot(StrainE1_Stat,StressE1_2_New); hold on;
plot(StrainE1_Stat,StressE1_3_New); hold on;
plot(StrainE1_Stat,StressE1_4_New)


n =1;
for i = 1:length(StressE1_2_New)
    Stress_calulate = [StressE1_2_New(i),StressE1_3_New(i),StressE1_4_New(i)];
   Stress_mean(n) = mean(Stress_calulate) % Mean Of All Experiments At Each Value Of ‘x’
   Stress_mean(i) = Stress_mean(n);
   SEM_Stress(n) = std(Stress_calulate)/sqrt(length(Stress_calulate)); % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
   n = n+1;
end
CI95 = tinv([0.05 0.95], length(Stress_calulate)-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, SEM_Stress, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

p = polyfit(StrainE1_Stat,Stress_mean,1);
fit = polyval(p,StrainE1_Stat);

p_manual = [65 0]
fit_manual = polyval(p_manual,StrainE1_Stat);

 close all
fig=figure('units','inch','position',[0,0,3.6,2.5]); hold on; grid on; set(gca,'FontSize',10);

plot(StrainE1_Stat,Stress_mean,'LineWidth',1.5,'Color',[0.75,0.75,0.75]);

a = errorbar(StrainE1_Stat,Stress_mean,yCI95(2,:),'o','MarkerSize',4,...
    'MarkerEdgeColor',[0.75,0.75,0.75],...
    'MarkerFaceColor',[0, 0.4470, 0.7410],'LineWidth',0.5,'Color',[0, 0.4470, 0.7410])

plot(StrainE1_Stat,fit_manual,'LineWidth',0.85,'Color',[0.8,0,0]);

 ylim([0 4]); set(gca,'YTick',[0:1:4]);
% 
xlim([0 0.052]); set(gca,'XTick',[0:0.01:0.05]);%set(gca,'xTick',[0:50:250]);


%  legend([a b c d e],'m = 0.1 kg','m = 0.2 kg','m = 0.3 kg','m = 0.4 kg','m = 0.5 kg')


ylabel('Stress, \sigma (MPa)')
xlabel('Strain, ɛ')
 
grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color
print(gcf,'E1.png','-dpng','-r700');   
