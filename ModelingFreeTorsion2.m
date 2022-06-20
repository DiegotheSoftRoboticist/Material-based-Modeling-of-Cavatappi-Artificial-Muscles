close all
clear all
clc

% Retrieve Data Sample 1
    filename = '42deg FreeTorsion Sample 1.xlsx';
    sheet = 'Sheet1';
    FT_1 = xlsread(filename, sheet,'C:E');
        timeFT_1 = FT_1(:,1);  
        PressureFT_1 = FT_1(:,2);
        DispFT_1 = FT_1(:,3);
        timeFT_1_p = timeFT_1+0.1;
        L1 = 16; %mm
        r = 3; %mm Spool radius
        r_tube = 1.85/2; %mm Outter Radius tube
        
        ThetaFT_1 = (DispFT_1)/(r); %rad
        ThetaFT_1_degree = (ThetaFT_1-ThetaFT_1(1))*(180/pi);
        GammaFT_1 = (ThetaFT_1*r_tube)/(L1); %rad/rad

% Retrieve Data Sample 2
    filename = '42deg FreeTorsion Sample 2.xlsx';
    sheet = 'Sheet1';
    FT_2 = xlsread(filename, sheet,'C:E');
        timeFT_2 = FT_2(:,1);  
        PressureFT_2 = FT_2(:,2);
        DispFT_2 = FT_2(:,3);
        
        L2 = 12; %mm
        
        ThetaFT_2 = (DispFT_2)/(r); %rad
        GammaFT_2 = (ThetaFT_2*r_tube)/(L2); %rad/rad
        
PressureFT_1 = PressureFT_1+386.35;
PressureFT_2 = PressureFT_2+386.35;

fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);

yyaxis right;
H1 = plot(timeFT_1_p,PressureFT_1,':','LineWidth',1,'Color',[0.5,0.5,0.5]); %plot(time4003(1:(end)),-Strain400_filt3(1:(end)),':','Linewidth',1.5)
ylabel('Pressure, (MPa)')
set(gca,'ycolor',[0.5,0.5,0.5])
% ylim([0 60])
yyaxis left;
H2 =  plot(timeFT_1,ThetaFT_1_degree,'Linewidth',1,'Color',[0, 0, 0]);% plot(time4003(1:(end)),p400_filt3(1:(end)),':','Linewidth',1.5)
ylabel('Shear strain, \gamma_{z\theta}')
set(gca,'ycolor',[0, 0, 0])
xlabel('Time (s)')
% xlim([0,8]);
% set(gca,'XTick',[0:1:8]);
% title('400 grams')
% legend('SAMPLE 1','SAMPLE 1','Location','Northwest')
grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color

p = [0:0.1:1.4]; %MPa which is equal to 200 psi
po = 0.1; %Mpa
P = 1.4;
r = (1.85/2)/1000; %m
ri = (0.85/2)/1000; %m
r_ave = (r+ri)/2;
t = 0.0005; %m 
a = 42; %deg

u_12 = 0.19;
u_23 = 0.32;
%  E1 = 50; %MPa
E2 = 3; %MPa
G12 = 6; %MPa
G23 = 6; %MPa
E1 = 31.72; %MPa
% E2 = 2.09; %MPa
% G12 = 2.78; %MPa
% G12 = E1/(2*(1+u_12)); %MPa
u_21 = (E2*u_12)/E1;

S = [1/E1 -(u_12/E1) -(u_12/E1) 0 0 0;
    -(u_12/E1) 1/E2 -(u_23/E2)  0 0 0;
    -(u_12/E1) -(u_23/E2) 1/E2 0 0 0;
    0 0 0 1/G23 0 0;
    0 0 0 0 1/G12 0;
    0 0 0 0 0 1/G12];

T = [cosd(a)^2 sind(a)^2 0 0 0 2*sind(a)*cosd(a);
    sind(a)^2 cosd(a)^2 0 0 0 -2*sind(a)*cosd(a);
    0 0 1 0 0 0;
    0 0 0 cosd(a) -sind(a) 0;
    0 0 0 cosd(a) sind(a) 0;
    -sind(a)*cosd(a) sind(a)*cosd(a) 0 0 0 (cosd(a)^2)-(sind(a)^2)];

 S_bar = transpose(T)*S*T;
 
%Epsilon_z is assumed to be 0.

Sigma_z = (p*ri^2)/(r^2-ri^2);
Term2 =((p)*r^2*ri^2)/((r^2-ri^2)*r^2);
Sigma_theta = Sigma_z + Term2;
Sigma_rad = Sigma_z - Term2;
Tau_z_theta = 0;

R = [(0.85/2)/1000:0.00005:(1.85/2)/1000];

Sigma_z_R = (P*ri^2)/(r^2-ri^2);
Term2_R =(P*r^2*ri^2)./((R.^2)*(r^2-ri^2));
Sigma_theta_R = Sigma_z_R + Term2_R;
Sigma_rad_R = Sigma_z_R - Term2_R;
Tau_z_theta_R = 0;


Epsilon_z = S_bar(1,1)*Sigma_z + S_bar(1,2)*Sigma_theta + S_bar(1,3)*Sigma_rad + S_bar(1,6)*Tau_z_theta;

Epsilon_theta = S_bar(2,1)*Sigma_z + S_bar(2,2)*Sigma_theta + S_bar(2,3)*Sigma_rad + S_bar(2,6)*Tau_z_theta;

Epsilon_rad = S_bar(3,1)*Sigma_z + S_bar(3,2)*Sigma_theta + S_bar(3,3)*Sigma_rad + S_bar(3,6)*Tau_z_theta;

Gamma_z_theta = S_bar(1,6)*Sigma_z + S_bar(2,6)*Sigma_theta + S_bar(3,6)*Sigma_rad + S_bar(6,6)*Tau_z_theta;



ThetaExperimental = [0 2.5 4.5 7 10 15 18.5 22.5 27 32.5];
GammExperimental = ThetaExperimental*0.185*(2/3)*(pi/180);

set(groot, 'DefaultTextInterpreter', 'tex', ...
           'DefaultAxesTickLabelInterpreter', 'tex', ...
           'DefaultAxesFontName', 'tex', ...
           'DefaultLegendInterpreter', 'tex', ...
           'defaultFigureColor','w');


  PressureFT_1MPa =PressureFT_1/145.038;
  PressureFT_2MPa =PressureFT_2/145.038;

fig=figure('units','inch','position',[0,0,4.5,3.75]); hold on; grid on; set(gca,'FontSize',10);
 
plot(p(2:end),Gamma_z_theta(2:end)-Gamma_z_theta(2),':','Linewidth',1,'Color',[0.5,0.5,0.5])
% plot(p,GammExperimental,'-','Linewidth',1,'Color',[1,0.5,0.5])
% plot(PressureFT_1MPa(39075:40678),(-GammaFT_1(39075:40678)+GammaFT_1(39075)))
plot(PressureFT_1MPa(16706:17500),(-GammaFT_1(16706:17500)+GammaFT_1(16706)))
plot(PressureFT_2MPa(38974:40025),(-GammaFT_2(38974:40025)+GammaFT_2(38974)))
% plot(p,Epsilon_z,':','Linewidth',1,'Color',[0.5,1,0.25])

% xlim([0 200])
% ylim([-2 40]);set(gca,'YTick',[0:10:40]); 
legend('\gamma_{z\theta} Model','\gamma_{z\theta} Experimental SAMPLE 1','\gamma_{z\theta} Experimental SAMPLE 2')
 xlabel('Pressure, (MPa)')
 ylabel('Shear strain, \gamma_{z\theta}')


grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 
% print(gcf,'optimalalpha0.01deg.png','-dpng','-r700');

fig2=figure('units','inch','position',[0,0,4.5,3.75]); hold on; grid on; set(gca,'FontSize',10);
 
plot(R*1000,Sigma_z_R*ones(1,11),':','Linewidth',1,'Color',[0.5,0.5,0.5])
plot(R*1000,Sigma_theta_R,'-','Linewidth',1,'Color',[1,0.5,0.5])
plot(R*1000,Sigma_rad_R,':','Linewidth',1,'Color',[0.5,1,0.25])

 xlim([(0.85/2) 1.85/2])
% ylim([-2 40]);set(gca,'YTick',[0:10:40]); 
legend('\sigma_{z}','\sigma_{\theta}','\sigma_{r}')
 xlabel('Radius, (mm)')
 ylabel('Stresses \sigma')


 grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 
% print(gcf,'optimalalpha0.01deg.png','-dpng','-r700');
