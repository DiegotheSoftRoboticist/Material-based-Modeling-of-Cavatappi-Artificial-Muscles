%Comparing thick vs. thin-walled thickness vessel

clear all
clc
close all

p = [0:0.14:2]; %MPa which is equal to 200 psi
po = 0; %Mpa
r = (2/2); %mm
ri = (0.9/2); %mm
r_ave = (r+ri)/2; %mm
t = 0.0005; %m 

%Mechanical Properties
u_12 = 0.19;
u_23 = 0.55;
E1 = 60; %MPa
E2 = 8.5; %MPa %very sensitive parameter
G12 = 6; %MPa
u_21 = (E2*u_12)/E1;

% %Stresses in a thin-walled thickness vessel
% Sigma_z_thin = (p*r)/(2*t);
% Sigma_theta_thin = (p*r)/t;
% Sigma_radial_thin = 0;

% %Stresses in a thick-walled thickness vessel
% Sigma_z_thick = (p*ri^2-po*r^2)/(r^2-ri^2);
%    Term2 =((p-po)*r^2*ri^2)/((r^2-ri^2)*r^2);
% Sigma_theta_thick = Sigma_z_thick + Term2;
% Sigma_radial_thick = Sigma_z_thick - Term2;

%     Term2_ave =((p-po)*r^2*ri^2)/((r^2-ri^2)*r_ave^2);
% Sigma_theta_thick_ave = Sigma_z_thick + Term2_ave;
% Sigma_radial_thick_ave = Sigma_z_thick - Term2_ave;
% 
%     Term2_i =((p-po)*r^2*ri^2)/((r^2-ri^2)*ri^2);
% Sigma_theta_thick_i = Sigma_z_thick + Term2_i;
% Sigma_radial_thick_i = Sigma_z_thick - Term2_i;
% 
% 
% Sigma_z = (p*ri^2-po*r^2)/(r^2-ri^2)
% Sigma_theta = Sigma_z + ((p-po)*r^2*ri^2)/((r^2-ri^2)*r^2)
% Sigma_radial = Sigma_z - ((p-po)*r^2*ri^2)/((r^2-ri^2)*r^2)
% 
% %Strains in a thin-walled thickness vessel
% e_z_thin = (1/E1)*Sigma_z_thin - (u_12/E1)*Sigma_theta_thin - (u_12/E1)*Sigma_radial_thin;
% e_theta_thin = -(u_12/E1)*Sigma_z_thin + (1/E2)*Sigma_theta_thin - (u_23/E2)*Sigma_radial_thin;
% e_rad_thin = -(u_12/E1)*Sigma_z_thin - (u_23/E2)*Sigma_theta_thin + (1/E2)*Sigma_radial_thin;
% 
% %Strains in a thick-walled thickness vessel
% e_z_thick = (1/E1)*Sigma_z_thick - (u_12/E1)*Sigma_theta_thick - (u_12/E1)*Sigma_radial_thick;
% e_theta_thick = -(u_12/E1)*Sigma_z_thick + (1/E2)*Sigma_theta_thick - (u_23/E2)*Sigma_radial_thick;
% e_rad_thick = -(u_12/E1)*Sigma_z_thick - (u_23/E2)*Sigma_theta_thick + (1/E2)*Sigma_radial_thick;
% 
% e_z_thick_ave = (1/E1)*Sigma_z_thick - (u_12/E1)*Sigma_theta_thick_ave - (u_12/E1)*Sigma_radial_thick_ave;
% e_theta_thick_ave = -(u_12/E1)*Sigma_z_thick + (1/E2)*Sigma_theta_thick_ave - (u_23/E2)*Sigma_radial_thick_ave;
% e_rad_thick_ave = -(u_12/E1)*Sigma_z_thick - (u_23/E2)*Sigma_theta_thick_ave + (1/E2)*Sigma_radial_thick_ave;
% 
% e_z_thick_i = (1/E1)*Sigma_z_thick - (u_12/E1)*Sigma_theta_thick_i - (u_12/E1)*Sigma_radial_thick_i;
% e_theta_thick_i = -(u_12/E1)*Sigma_z_thick + (1/E2)*Sigma_theta_thick_i - (u_23/E2)*Sigma_radial_thick_i;
% e_rad_thick_i = -(u_12/E1)*Sigma_z_thick - (u_23/E2)*Sigma_theta_thick_i + (1/E2)*Sigma_radial_thick_i;
% 
% 
% 
% set(groot, 'DefaultTextInterpreter', 'tex', ...
%            'DefaultAxesTickLabelInterpreter', 'tex', ...
%            'DefaultAxesFontName', 'tex', ...
%            'DefaultLegendInterpreter', 'tex', ...
%            'defaultFigureColor','w');
% 
% e_z_R_Experimental =[0 -0.027 -0.24 -0.55 -0.82 -1.12 -1.27 -1.29 -1.24 -1.17 -1.13];
% 
% e_theta_R_Experimental = [0 0.50 4.90 10.40 13.90 18.90 22.32 25.65 28.20 30.950 34];
% 
% p_Experimental_radial = [0 9.04 37 58 73 102 121 142 160 182 207];
% 
% p_Experimental_axial = [0 9 38 59.92 74.29 102.3 123.4 145 163.5 183 192.5];
% 
% 
% fig1=figure('units','inch','position',[0,0,4.5,3.75]); hold on; grid on; set(gca,'FontSize',10);
% plot(p*140,e_z_thin,':','Linewidth',1,'Color',[0.5,0.5,0.5])
% plot(p*140,e_theta_thin,'-','Linewidth',1,'Color',[1,0.5,0.5])
% plot(p*140,e_rad_thin,'-','Linewidth',1,'Color',[0,0.25,0])
% 
% plot(p_Experimental_radial,e_z_R_Experimental/100,'+','Linewidth',1,'Color',[0,0,0])
% plot(p_Experimental_axial,e_theta_R_Experimental/100,'*','Linewidth',1,'Color',[0,0,0])
% % xlim([0 200])
% % ylim([-2 40]);set(gca,'YTick',[0:10:40]); 
% 
% legend('e_z thin','e_{theta} thin','e_{rad} thin','e_{theta} Experimental r = R_{outer}','e_{z} Experimental')
% xlabel('Pressure, (MPa)')
% ylabel('Strains')
% 
% 
%  grid on 
% set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 
% 
% e_z_R_Experimental =[0 -0.027 -0.24 -0.55 -0.82 -1.12 -1.27 -1.29 -1.24 -1.17 -1.13];
% 
% e_theta_R_Experimental = [0 0.4 3 8.3 10.42 14.20 16.53 19.23 21.15 24.05 27];
% 
% p_Experimental_radial = [0 9.04 37 58 73 102 121 142 160 182 207];
% 
% p_Experimental_axial = [0 9 38 59.92 74.29 102.3 123.4 145 163.5 183 192.5];
% 
%     
% 
% fig2=figure('units','inch','position',[0,0,4.5,3.75]); hold on; grid on; set(gca,'FontSize',10); 
% plot(p*140,e_z_thick,':','Linewidth',1,'Color',[0.5,0.3,0.5])
% plot(p*140,e_theta_thick,'-.','Linewidth',1,'Color',[1,0.1,0.5])
% plot(p*140,e_rad_thick,'-','Linewidth',1,'Color',[0,0.25,0])
% 
% % plot(p*140,e_z_thick_ave,':','Linewidth',1,'Color',[0.5,0.4,0.5])
% plot(p*140,e_theta_thick_ave,'-.','Linewidth',1,'Color',[1,0.4,0.5])
% plot(p*140,e_rad_thick_ave,'-','Linewidth',1,'Color',[0,0.7,0])
% 
% % plot(p*140,e_z_thick_i,':','Linewidth',1,'Color',[0.5,0.5,0.5])
% plot(p*140,e_theta_thick_i,'-.','Linewidth',1,'Color',[1,0.9,0.5])
% plot(p*140,e_rad_thick_i,'-','Linewidth',1,'Color',[0,1,0])
% 
% plot(p_Experimental_radial+14,e_z_R_Experimental/100,'+','Linewidth',1,'Color',[0,0,0])
% plot(p_Experimental_axial+14,e_theta_R_Experimental/100,'*','Linewidth',1,'Color',[0,0,0])
% 
%  xlim([0 200])
% % ylim([-2 40]);set(gca,'YTick',[0:10:40]); 
% legend('e_z thick','e_{theta} thick at r = R_{outer}','e_{rad} thick at r = R_{outer}',...
%    'e_{theta} thick at r = R_{average}','e_{rad} thick r = R_{average}',...
%    'e_{theta} thick r = R_{inner}','e_{rad} thick r = R_{inner}'...
%    ,'e_{theta} Experimental r = R_{outer}','e_{z} Experimental')
%  xlabel('Pressure, (psi)')
%  ylabel('Strains')
% 
% 
%  grid on 
% set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 
% 
% print(gcf,'Notwistedsample Strains.png','-dpng','-r700');
% % 
u = 2*(((ri^2)*p*r)/(E2*(r^2-ri^2)));
fig3=figure('units','inch','position',[0,0,4.5,3.75]); hold on; grid on; set(gca,'FontSize',10);
 
plot(p*140,u,'Linewidth',1,'Color',[0.5,0.5,0.5])
% plot(p_Experimental_radial(1:9)+14,(e_theta_R_Experimental(1:9)*1.8/100)/2,'+','Linewidth',1,'Color',[0,0,0])



 xlim([0 250])
% ylim([-2 40]);set(gca,'YTick',[0:10:40]); 
legend('u (mm)','e_{theta} thick','e_{rad} thick')
 xlabel('Pressure, (MPa)')
 ylabel('Stresses')


 grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 


