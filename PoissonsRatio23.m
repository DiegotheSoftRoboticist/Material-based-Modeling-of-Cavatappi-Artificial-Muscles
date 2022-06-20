%Poisson's Ratio%
clear all 
close all
clc

% Retrieve Data Sample 1
% % % % % % % %     filename = 'Poissons ratio_23. Increments of 0.025ml.xlsx';
% % % % % % % %     sheet = 'Sheet1';
% % % % % % % %     FT_1 = xlsread(filename, sheet,'C:E');
% % % % % % % %         timeFT_1 = FT_1(:,1);  
% % % % % % % %         PressureFT_1 = FT_1(:,2);
        
        L = 77; %mm
        
        R = 2/2; %mm Outter Radius tube
        ri = 0.8/2;  %mm Inner Radius tube
        
E2 = 9; %MPa

Vol_ml = [0:0.02:0.08];  
Vol_0 = 2*ri*L*pi;
Vol_f = Vol_0+Vol_ml*1000;
rf = Vol_f/(2*L*pi);   
u_ri = rf-ri;
p = [13 65 135 200 255]/145.038; %MPa
    
% % % fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);
% % % plot(timeFT_1,PressureFT_1,':','Linewidth',1,'Color',[0.5,0.5,0.5])
% % % % xlim([0 200])
% % % % ylim([-2 40]);set(gca,'YTick',[0:10:40]); 
% % %  xlabel('Time, (s)')
% % %  ylabel('Pressure, (psi)')


 grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 

    n =1;    
      for i = 1:length(p) 
syms u_23
        
eqn = u_ri(i) == (((ri^3)*(p(i)))/(E2*(R^2-ri^2)))*((1-u_23)+(1+u_23)*(R^2/ri^2));
 
S(i,:) = solve(eqn,u_23); 
u_23(n) = u_23;
n = n+1;

      end   
      
u_23 = S(:,1)


%https://www.makeitfrom.com/material-properties/Plasticized-Flexible-Polyvinyl-Chloride-PVC-P
%https://www.shimadzu.com/an/industries/engineering-materials/film/poisson/index.html
