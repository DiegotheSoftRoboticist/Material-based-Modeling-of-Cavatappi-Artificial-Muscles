close all
clear all
clc

% Retrieve Data Sample 1
    filename = '15deg FreeTorsion Sample 1.xlsx';
    sheet = 'Sheet1';
    FT_1 = xlsread(filename, sheet,'C:E');
        timeFT_1 = FT_1(:,1);  
        PressureFT_1 = FT_1(:,2);
        DispFT_1 = FT_1(:,3);
        timeFT_1_p = timeFT_1;
        L1 = 8; %mm
        r_spool = 3; %mm Spool radius
        r_tube = 2.30/2; %mm Outter Radius tube
        
        ThetaFT_1 = (DispFT_1)/(r_spool); %rad
        ThetaFT_1_degree = (ThetaFT_1-ThetaFT_1(1))*(180/pi);
        GammaFT_1 = (ThetaFT_1*r_tube)/(L1); %rad/rad
% 
% % Retrieve Data Sample 2
    filename = '15deg FreeTorsion Sample 2.xlsx';
    sheet = 'Sheet1';
    FT_2 = xlsread(filename, sheet,'C:E');
        timeFT_2 = FT_2(:,1);  
        PressureFT_2 = FT_2(:,2);
        DispFT_2 = FT_2(:,3);
        
        L2 = 13; %mm
        
        ThetaFT_2 = (DispFT_2)/(r_spool); %rad
        GammaFT_2 = (ThetaFT_2*r_tube)/(L2); %rad/rad
        
% PressureFT_1 = PressureFT_1+386.35;
% PressureFT_2 = PressureFT_2+386.35;

fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);

yyaxis right;
H1 = plot(timeFT_2,PressureFT_2,':','LineWidth',1,'Color',[0.5,0.5,0.5]); %plot(time4003(1:(end)),-Strain400_filt3(1:(end)),':','Linewidth',1.5)
ylabel('Pressure, (MPa)')
set(gca,'ycolor',[0.5,0.5,0.5])
% ylim([0 60])
yyaxis left;
H2 =  plot(timeFT_2,-GammaFT_2,'Linewidth',1,'Color',[0, 0, 0]);% plot(time4003(1:(end)),p400_filt3(1:(end)),':','Linewidth',1.5)
ylabel('Shear strain, \gamma_{z\theta}')
set(gca,'ycolor',[0, 0, 0])
xlabel('Time (s)')
% xlim([0,8]);
% set(gca,'XTick',[0:1:8]);
% title('400 grams')
% legend('SAMPLE 1','SAMPLE 1','Location','Northwest')
grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color

Pi = [0:0.01:1.4]; %MPa which is equal to 200 psi\
P0 = 0; %Mpa

r0 = (2/2); %mm
ri = (0.8/2); %mm
r = ri:0.0001:r0; %mm
d0 = r0*2; %mm
di =ri*2; %mm
R = 1.5; %mm
a = 10; %deg

%MATERIAL PROPERTIES OF THE STRAIGHT DRAWN TUBE

u_12 = 0.19;
u_23 = 0.55;
E2 = 9.1; %MPa
G12 = 10; %MPa
G23 = 3.2; %MPa
E1 = 61.1; %MPa
u_21 = (E2*u_12)/E1;

%%%%%TRANSVERSILY ISOTROPIC MODEL%%%

%COMPLIANCE MATERIX
S = [1/E1 -(u_12/E1) -(u_12/E1) 0 0 0;
    -(u_12/E1) 1/E2 -(u_23/E2)  0 0 0;
    -(u_12/E1) -(u_23/E2) 1/E2 0 0 0;
    0 0 0 1/G23 0 0;
    0 0 0 0 1/G12 0;
    0 0 0 0 0 1/G12];

S_bar(1,1) = (cosd(a)^4)*S(1,1) + (cosd(a)^2)*(sind(a)^2)*(2*S(1,2)*S(6,6)) + (sind(a)^4)*S(2,2);
S_bar(1,2) = (cosd(a)^2)*(sind(a)^2)*(S(1,1)+S(2,2)-S(6,6)) + ((cosd(a)^4)+(sind(a)^4))*S(1,2);
S_bar(1,3) = (cosd(a)^2)*S(1,3) + (sind(a)^2)*S(2,3); 
S_bar(1,6) = (cosd(a))*(sind(a))*((cosd(a)^2)-(sind(a)^2))*(2*S(1,2)+S(6,6)) + 2*(cosd(a))*(sind(a))*((sind(a)^2)*S(2,2)-(cosd(a)^2)*S(1,1));
S_bar(2,2) = (cosd(a)^4)*S(2,2) + (cosd(a)^2)*(sind(a)^2)*(2*S(1,2)*S(6,6)) + (sind(a)^4)*S(1,1);
S_bar(3,3) = S(3,3);
S_bar(2,3) = (cosd(a)^2)*S(2,3) + (sind(a)^2)*S(1,3); 
S_bar(2,6) = (cosd(a))*(sind(a))*((sind(a)^2)-(cosd(a)^2))*(2*S(1,2)+S(6,6)) + 2*(cosd(a))*(sind(a))*((cosd(a)^2)*S(2,2)-(sind(a)^2)*S(1,1));
S_bar(6,6) = 4*(cosd(a)^2)*(sind(a)^2)*(S(1,1)+S(2,2)-2*S(1,2)) + (((cosd(a)^2)-(sind(a)^2))^2)*S(6,6);
S_bar(3,6) = 2*(cosd(a))*(sind(a))*(S(2,3)-S(1,3));
 
for i = 1:length(Pi)
%%%%%THICK WALL VESSEL. THEORY%%%
r = ri:0.001:r0; %mm
u_r = (((ri^2)*(Pi(i))*r)./(E2*(r0^2-ri^2))).*((1-u_23)+(1+u_23)*(r0^2./r.^2));

Sigma_z = (Pi(i)*(ri+u_r(1))^2-P0*(r0+u_r(end))^2)/((r0+u_r(end))^2-(ri+u_r(1))^2); 
Term2 =((Pi(i)-P0)*(r0+u_r(end))^2*(ri+u_r(1))^2)./(((r0+u_r(end))^2-(ri+u_r(1))^2)*(u_r+r).^2);
Sigma_theta = Sigma_z + Term2;
Sigma_rad = Sigma_z - Term2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Plot CHECK%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig=figure('units','inch','position',[0,0,4.5,3.75]); hold on; grid on; set(gca,'FontSize',10);
% Sigma_z(1:length(Term2))=Sigma_z;
% plot(r,Sigma_z,'-.','Linewidth',1,'Color',[0.5,0.5,0.5]);
% plot(r,Sigma_theta,'Linewidth',1,'Color',[0,0,0]);
% plot(r,Sigma_rad,'.','Linewidth',1,'Color',[1,0.5,0.5]);
% plot(r,Term2,'Linewidth',1,'Color',[0.5,0.5,1]);
% legend('Sigma_z','Sigma_theta','Sigma_rad',...
%     'Term2')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = (pi*(((r0+u_r(end)).^4)-((ri+u_r(1)).^4)))/2;
A = (pi*(((r0+u_r(end)).^2)-((ri+u_r(1)).^2)));

%SHEAR STRESS PROFILE


F_iso = 0; % N. Load of the mass

tau_torsion = ((F_iso*(r_spool)*(r0+u_r(end)))./J);
tau_direct = 0;

tau(i) = tau_torsion+tau_direct;


gammaZ_oi(i) = S_bar(1,6)*Sigma_z(1);
gammaTheta_oi(i) = S_bar(2,6)*Sigma_theta(end);
gammaR_oi(i) = S_bar(3,6)*Sigma_rad(end);
gammaTau_oi(i) = S_bar(6,6)*tau_torsion;%-S_bar(6,6)*Tau(1);
gammaTau_SUB_oi(i) = gammaTau_oi(i)-gammaTau_oi(1);
Gamma_z_theta_oi(i) = gammaZ_oi(i) + gammaTheta_oi(i) + gammaR_oi(i) - abs(gammaTau_SUB_oi(i));



% Gamma_z_theta_avg(i) = Gamma_z_theta;
end

set(groot, 'DefaultTextInterpreter', 'tex', ...
           'DefaultAxesTickLabelInterpreter', 'tex', ...
           'DefaultAxesFontName', 'tex', ...
           'DefaultLegendInterpreter', 'tex', ...
           'defaultFigureColor','w');


  PressureFT_1MPa =PressureFT_1/145.038;
 PressureFT_2MPa =PressureFT_2/145.038;

fig=figure('units','inch','position',[0,0,3.6,3.2]); hold on; grid on; set(gca,'FontSize',10);
 
% plot(p,GammExperimental,'-','Linewidth',1,'Color',[1,0.5,0.5])
% plot(PressureFT_1MPa(39075:40678),(-GammaFT_1(39075:40678)+GammaFT_1(39075)))
 plot(PressureFT_1MPa(27298-200:100:30364-200+4500)-PressureFT_1MPa(27298-200),(-GammaFT_1(27298:100:30364+4500)+GammaFT_1(27298)),'-','Linewidth',1,'Color',[0 0.4 0.9])
 plot(PressureFT_2MPa(14700-100:40:16100-100+1000)-PressureFT_2MPa(14700-100),(-GammaFT_2(14700:40:16100+1000)+GammaFT_2(14700)),'-','Linewidth',1,'Color',[0.2 0.8 0.8])
% plot(p,Epsilon_z,':','Linewidth',1,'Color',[0.5,1,0.25])
plot(Pi(1:end),Gamma_z_theta_oi(1:end)-Gamma_z_theta_oi(1),':','Linewidth',2.3,'Color',[0,0,0])

 xlim([0 1.5])
 ylim([-0.01 0.07])
 
%  set(gca,'YTick',[0:10:40]); 
legend('Sample 1',...
    'Sample 2','Model')
 xlabel('Pressure, (MPa)')
 ylabel('Shear strain, \gamma_{zθ}')


grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 
 print(gcf,'TorsionalCavatappi_15deg_FreeTorsion','-dpng','-r700');

 
 fig1=figure('units','inch','position',[0,0,3.6,3.2]); hold on; grid on; set(gca,'FontSize',10);
 
% plot(p,GammExperimental,'-','Linewidth',1,'Color',[1,0.5,0.5])
% plot(PressureFT_1MPa(39075:40678),(-GammaFT_1(39075:40678)+GammaFT_1(39075)))
 plot(PressureFT_1MPa(27298-200:100:30364-200)-PressureFT_1MPa(27298-200),(-GammaFT_1(27298:100:30364)+GammaFT_1(27298)),'-','Linewidth',1,'Color',[0 0.4 0.9])
 plot(PressureFT_2MPa(14700-100:40:16100-100+200)-PressureFT_2MPa(14700-100),(-GammaFT_2(14700:40:16100+200)+GammaFT_2(14700)),'-','Linewidth',1,'Color',[0.2 0.8 0.8])
% plot(p,Epsilon_z,':','Linewidth',1,'Color',[0.5,1,0.25])
plot(Pi(1:end),Gamma_z_theta_oi(1:end)-Gamma_z_theta_oi(1),':','Linewidth',2.3,'Color',[0,0,0])

 xlim([0 1.5])
 ylim([-0.01 0.07])
 
%  set(gca,'YTick',[0:10:40]); 
legend('Sample 1',...
    'Sample 2','Model')
 xlabel('Pressure, (MPa)')
 ylabel('Shear strain, \gamma_{zθ}')

%\alpha = 12^\circ
grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 
 print(gcf,'TorsionalCavatappi_15deg_FreeTorsionLOADING','-dpng','-r700');
 
 
% fig2=figure('units','inch','position',[0,0,4.5,3.75]); hold on; grid on; set(gca,'FontSize',10);
%  
% plot(R*1000,Sigma_z_R*ones(1,11),':','Linewidth',1,'Color',[0.5,0.5,0.5])
% plot(R*1000,Sigma_theta_R,'-','Linewidth',1,'Color',[1,0.5,0.5])
% plot(R*1000,Sigma_rad_R,':','Linewidth',1,'Color',[0.5,1,0.25])
% 
%  xlim([(0.85/2) 1.85/2])
% % ylim([-2 40]);set(gca,'YTick',[0:10:40]); 
% legend('\sigma_{z}','\sigma_{\theta}','\sigma_{r}')
%  xlabel('Radius, (mm)')
%  ylabel('Stresses \sigma')
% 
% 
%  grid on 
% set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 
% % % print(gcf,'optimalalpha0.01deg.png','-dpng','-r700');
