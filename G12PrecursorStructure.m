%G12 Precursor Structure by Diego Higueras-Ruiz
 
clear all
close all
clc

do = (1.8); %mm
        di = (0.85); %mm
        J = ((do^4-di^4)*pi)/32;
        
% Retrieve Data Sample 1
    filename = 'Shear Sample 1 Clean.xlsx';
    sheet = 'Sheet1';
    G12Sample1 = xlsread(filename, sheet,'A:C');
        timeG12_1 = G12Sample1(:,1);  
        TorqueG12_1 = G12Sample1(:,2);
        ThetaG12_1 = G12Sample1(:,3);

        L1 = 9; %mm
        SStressG12_1 = (TorqueG12_1*(do/2))/(J*1000); %MPa
        SStrainG12_1 = ((ThetaG12_1-ThetaG12_1(1))*(do/2))/L1; %Strain
        
  % Retrieve Data Sample 2
      filename = 'Shear Sample 2 Clean.xlsx';
    sheet = 'Sheet1';
    G12Sample2 = xlsread(filename, sheet,'A:C');
        timeG12_2 = G12Sample2(:,1);  
        TorqueG12_2 = G12Sample2(:,2);
        ThetaG12_2 = G12Sample2(:,3);
        
        L2 = 10; %mm
        SStressG12_2 = (TorqueG12_2*(do/2))/(J*1000); %MPa
        SStrainG12_2 = ((ThetaG12_2-ThetaG12_2(1))*(do/2))/L1; %Strain
        
   % Retrieve Data Sample 3
    filename = 'Shear Sample 3 Clean.xlsx';
    sheet = 'Sheet1';
    G12Sample3 = xlsread(filename, sheet,'A:C');
        timeG12_3 = G12Sample3(:,1);  
        TorqueG12_3 = G12Sample3(:,2);
        ThetaG12_3 = G12Sample3(:,3);
        
        L3 = 8.5; %mm
        SStressG12_3 = (TorqueG12_3*(do/2))/(J*1000); %MPa
        SStrainG12_3 = ((ThetaG12_3-ThetaG12_3(1))*(do/2))/L1; %Strain    
        
   % Retrieve Data Sample 4
  filename = 'Shear Sample 4 Clean.xlsx';
    sheet = 'Sheet1';
    G12Sample4 = xlsread(filename, sheet,'A:C');
        timeG12_4 = G12Sample4(:,1);  
        TorqueG12_4 = G12Sample4(:,2);
        ThetaG12_4 = G12Sample4(:,3);
        
        L4 = 9; %mm
        SStressG12_4 = (TorqueG12_4*(do/2))/(J*1000); %MPa
        SStrainG12_4 = ((ThetaG12_4-ThetaG12_4(1))*(do/2))/L1; %Strain 
      
        % Retrieve Data Sample 5
  filename = 'Shear Sample 5 Clean.xlsx';
    sheet = 'Sheet1';
    G12Sample5 = xlsread(filename, sheet,'A:C');
        timeG12_5 = G12Sample5(:,1);  
        TorqueG12_5 = G12Sample5(:,2);
        ThetaG12_5 = G12Sample5(:,3);
        
        L5 = 9; %mm
        SStressG12_5 = (TorqueG12_5*(do/2))/(J*1000); %MPa
        SStrainG12_5 = ((ThetaG12_5-ThetaG12_5(1))*(do/2))/L1; %Strain 
        
        
set(groot, 'DefaultTextInterpreter', 'tex', ...
           'DefaultAxesTickLabelInterpreter', 'tex', ...
           'DefaultAxesFontName', 'tex', ...
           'DefaultLegendInterpreter', 'tex', ...
           'defaultFigureColor','w');
        
        
     %E1 Sample 1 Temporal plot
     
fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);

yyaxis right;
H1 = plot(timeG12_3,SStrainG12_3,':','LineWidth',1,'Color',[0.5,0.5,0.5]); %plot(time4003(1:(end)),-Strain400_filt3(1:(end)),':','Linewidth',1.5)
ylabel('Strain, \gamma_{12} (\circ/\circ)')
set(gca,'ycolor',[0.5,0.5,0.5])
% ylim([0 60])
yyaxis left;
H2 =  plot(timeG12_3,SStressG12_3,'Linewidth',1,'Color',[0, 0, 0]);% plot(time4003(1:(end)),p400_filt3(1:(end)),':','Linewidth',1.5)
ylabel('Stress, \tau_{12} (MPa)')
set(gca,'ycolor',[0, 0, 0])
xlabel('Time (s)')
% xlim([0,8]);
% set(gca,'XTick',[0:1:8]);
% title('400 grams')
% legend('SAMPLE 1','SAMPLE 1','Location','Northwest')
grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color

% print(gcf,'MaxStrain.png','-dpng','-r700'); 



%Stress vs. Strain

fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);

  plot(SStrainG12_1(955:1021)-SStrainG12_1(955), SStressG12_1(955:1021)-SStressG12_1(955),'Color',[0 0 0],'Linewidth',1)      
  plot(SStrainG12_2(955:1021)-SStrainG12_2(955), SStressG12_2(955:1021)-SStressG12_2(955),'Color',[0.1 0.5 0],'Linewidth',1)      
  plot(SStrainG12_3(955:1021)-SStrainG12_3(955), SStressG12_3(955:1021)-SStressG12_3(955),'Color',[0.5 0 0],'Linewidth',1)      
  plot(SStrainG12_4(955:1018)-SStrainG12_4(955), SStressG12_4(955:1018)-SStressG12_4(955),'Color',[0.5 1 0],'Linewidth',1)      
  plot(SStrainG12_5(955:1021)-SStrainG12_5(955), SStressG12_5(955:1021)-SStressG12_5(955),'Color',[0.5 0 1],'Linewidth',1)      

%              ylim([-0.5 1.5]);set(gca,'YTick',[0:0.5:1.5]); 
       ylabel('Stress, \tau_{12} (MPa)')
       xlabel('Strain, \gamma_{12} (\circ/\circ)')
%         leg = legend('ɛ_{Axial} Nylon Monofilament',...
%      'ɛ_{Radial} Nylon Monofilament','Location','northwest') 
%  leg.ItemTokenSize = [5,20];
%    grid on 
% set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 
% print(gcf,'Strains_DrawnandNONDrawn.png','-dpng','-r700');


%STATISTIC PLOT
StrainG12_2_Stat = SStrainG12_2(908:1021)-SStrainG12_2(908);
StressG12_2_Stat = SStressG12_2(908:1021)-SStressG12_2(908);
StrainG12_3_Stat = SStrainG12_3(908:1021)-SStrainG12_3(908);
StressG12_3_Stat = SStressG12_3(908:1021)-SStressG12_3(908);
StrainG12_4_Stat = SStrainG12_4(908:1021)-SStrainG12_4(908);
StressG12_4_Stat = SStressG12_4(908:1021)-SStressG12_4(908);

StrainG12_Stat = [0:0.005:0.05];

n = 1;
for i =1:length(StrainG12_Stat)
     if StrainG12_Stat(i) == 0
   StressG12_2_New(n) = 0;
   StressG12_3_New(n) = 0;
   StressG12_4_New(n) = 0;
     else
   StressG12_2_New(n) = interp1q(StrainG12_2_Stat,StressG12_2_Stat,StrainG12_Stat(i));
   StressG12_3_New(n) = interp1q(StrainG12_3_Stat,StressG12_3_Stat,StrainG12_Stat(i));
   StressG12_4_New(n) = interp1q(StrainG12_4_Stat,StressG12_4_Stat,StrainG12_Stat(i));
 
     end
 n = n+1;
 end
figure 
plot(StrainG12_Stat,StressG12_2_New); hold on;
plot(StrainG12_Stat,StressG12_3_New); hold on;
plot(StrainG12_Stat,StressG12_4_New)


n =1;
for i = 1:length(StressG12_3_New)
    Stress_calulate = [StressG12_2_New(i),StressG12_3_New(i),StressG12_4_New(i)];
   Stress_mean(n) = mean(Stress_calulate) % Mean Of All Experiments At Each Value Of ‘x’
   Stress_mean(i) = Stress_mean(n);
   SEM_Stress(n) = std(Stress_calulate)/sqrt(length(Stress_calulate)); % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
   n = n+1;
end
CI95 = tinv([0.05 0.95], length(Stress_calulate)-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, SEM_Stress, CI95(:)); % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’

p = polyfit(StrainG12_Stat,Stress_mean,1);
fit = polyval(p,StrainG12_Stat);
p_manual = [9.5 0]
fit_manual = polyval(p_manual,StrainG12_Stat);
p_23 = [3 0]
fit_23 = polyval(p_23,StrainG12_Stat); 

close all
 
fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',10);
 
plot(StrainG12_Stat,Stress_mean,'LineWidth',1.5,'Color',[0.75,0.75,0.75]);


 a = errorbar(StrainG12_Stat,Stress_mean,yCI95(2,:),'o','MarkerSize',4,...
    'MarkerEdgeColor',[0.75,0.75,0.75],...
    'MarkerFaceColor',[0, 0.4470, 0.7410],'LineWidth',0.5,'Color',[0, 0.4470, 0.7410])

plot(StrainG12_Stat,fit_manual,'LineWidth',0.85,'Color',[0.8,0,0]);
plot(StrainG12_Stat,fit_23,':','LineWidth',0.85,'Color',[0.8,0,0]);

ylim([0 0.6]); set(gca,'YTick',[0:0.1:0.6]);
 
xlim([0 0.052]); set(gca,'xTick',[0:0.01:0.052]);


%  legend([a b c d e],'m = 0.1 kg','m = 0.2 kg','m = 0.3 kg','m = 0.4 kg','m = 0.5 kg')



       ylabel('Stress, \tau (MPa)')
       xlabel('Strain, \gamma')
 
grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color
print(gcf,'G12.png','-dpng','-r700');   
