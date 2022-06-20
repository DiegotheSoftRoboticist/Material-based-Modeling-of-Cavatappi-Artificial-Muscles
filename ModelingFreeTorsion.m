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
        
        L1 = 15; %mm
        r = 3; %mm Spool radius
        r_tube = 1.85/2; %mm Outter Radius tube
        
        ThetaFT_1 = (DispFT_1)/(r); %rad
        GammaFT_1 = (ThetaFT_1*r_tube)/(L1); %rad/rad

% Retrieve Data Sample 2
    filename = '42deg FreeTorsion Sample 2.xlsx';
    sheet = 'Sheet1';
    FT_2 = xlsread(filename, sheet,'C:E');
        timeFT_2 = FT_2(:,1);  
        PressureFT_2 = FT_2(:,2);
        DispFT_2 = FT_2(:,3);
        
        L2 = 15; %mm
        
        ThetaFT_2 = (DispFT_2)/(r); %rad
        GammaFT_2 = (ThetaFT_2*r_tube)/(L2); %rad/rad
        
PressureFT_1 = PressureFT_1+386.35;
PressureFT_2 = PressureFT_2+386.35;

fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);

yyaxis right;
H1 = plot(timeFT_1,PressureFT_1,':','LineWidth',1,'Color',[0.5,0.5,0.5]); %plot(time4003(1:(end)),-Strain400_filt3(1:(end)),':','Linewidth',1.5)
ylabel('Pressure, (MPa)')
set(gca,'ycolor',[0.5,0.5,0.5])
% ylim([0 60])
yyaxis left;
H2 =  plot(timeFT_1,DispFT_1,'Linewidth',1,'Color',[0, 0, 0]);% plot(time4003(1:(end)),p400_filt3(1:(end)),':','Linewidth',1.5)
ylabel('Shear strain, \gamma_{z\theta}')
set(gca,'ycolor',[0, 0, 0])
xlabel('Time (s)')
% xlim([0,8]);
% set(gca,'XTick',[0:1:8]);
% title('400 grams')
% legend('SAMPLE 1','SAMPLE 1','Location','Northwest')
grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color

p = [0.1:0.14:1.4]; %MPa which is equal to 200 psi
po = 0.1; %Mpa
r = (1.85/2)/1000; %m
ri = (0.85/2)/1000; %m
r_ave = (r+ri)/2;
t = 0.0005; %m 
a = 32; %deg

u_12 = 0.19;
%  E1 = 50; %MPa
E2 = 2; %MPa
% G12 = 2.78; %MPa
G12 = 6; %MPa

E1 = 31.72; %MPa
% E2 = 2.09; %MPa
% G12 = E1/(2*(1+u_12)); %MPa
u_21 = (E2*u_12)/E1;

Q_11 = E1/(1-u_12*u_21);
Q_12 = (u_12*E2)/(1-u_12*u_21);
Q_22 = E2/(1-u_12*u_21);
Q_66 = G12;


Q_11_bar = Q_11*(cosd(a))^4 + Q_22*(sind(a))^4 + 2*(Q_12 + 2*Q_66)*(sind(a))^2*(cosd(a))^2;
Q_12_bar = (Q_11 + Q_22 - 4*Q_66)*(((sind(a))^2)*(cosd(a))^2) + Q_12*((cosd(a))^4 + (sind(a))^4); 
Q_22_bar = Q_11*(sind(a))^4 + Q_22*(cosd(a))^4 + 2*(Q_12 + 2*Q_66)*(sind(a))^2*(cosd(a))^2;
Q_16_bar = (Q_11 - Q_12 - 2*Q_66)*((cosd(a)))^3*(sind(a)) - (Q_22 - Q_12 - 2*Q_66)*(cosd(a))*(sind(a))^3;
Q_26_bar = (Q_11 - Q_12 - 2*Q_66)*(cosd(a))*(sind(a))^3 - (Q_22 - Q_12 - 2*Q_66)*((cosd(a))^3)*(sind(a));
Q_66_bar = (Q_11 + Q_22 - 2*Q_12 - 2*Q_66)*((sind(a))^2)*(cosd(a))^2 + Q_66*(((sind(a))^4)+(cosd(a))^4);



%Epsilon_z is assumed to be 0.



syms Epsilon_z Epsilon_theta Gamma_z_theta 
 for i = 1:length(p) 
     
%      Sigma_z(i) = (p(i)*r)/(2*t)
%      Sigma_theta(i) = (p(i)*r)/t


Sigma_z(i) = (p(i)*ri^2-po*r^2)/(r^2-ri^2)
Sigma_theta(i) = Sigma_z(i) + ((p(i)-po)*r^2*ri^2)/((r^2-ri^2)*r_ave^2)
Sigma_rad(i) = Sigma_z(i) - ((p(i)-po)*r^2*ri^2)/((r^2-ri^2)*r_ave^2)

equ1 = Sigma_z(i) == Q_11_bar*Epsilon_z + Q_12_bar*Epsilon_theta + Q_16_bar*Gamma_z_theta;

equ2 = Sigma_theta(i) == Q_12_bar*Epsilon_z + Q_22_bar*Epsilon_theta + Q_26_bar*Gamma_z_theta;

equ3 =    0 == Q_16_bar*Epsilon_z + Q_26_bar*Epsilon_theta + Q_66_bar*Gamma_z_theta;

[A,B] = equationsToMatrix([equ1, equ2, equ3], [Epsilon_z, Epsilon_theta, Gamma_z_theta]);
X(i,:) = linsolve(A,B);
% E_thetaANDG_z_theta(i,:) = B
% 
 end
% 
%  m =2
% 
% E1 = Sigma_z(m) == Q_11_bar*Epsilon_z + Q_12_bar*Epsilon_theta + Q_16_bar*Gamma_z_theta/2;
% 
% E2 = Sigma_theta(m) == Q_12_bar*Epsilon_z + Q_22_bar*Epsilon_theta + Q_26_bar*Gamma_z_theta/2;
% 
% E3 =    0 == Q_16_bar*Epsilon_z + Q_26_bar*Epsilon_theta + Q_66_bar*Gamma_z_theta/2;
% 
% result = solve(E1,E2,E3)
% 
% vpa(result.Epsilon_z)
% vpa(result.Epsilon_theta)
% vpa(result.Gamma_z_theta)
% 
% 
Gamma = X(:,3);
Epsilon_zz = X(:,1);
Epsilon_tt = X(:,2);

 
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
 
plot(p,(-Gamma+Gamma(1)),':','Linewidth',1,'Color',[0.5,0.5,0.5])
% plot(p,GammExperimental,'-','Linewidth',1,'Color',[1,0.5,0.5])
% plot(PressureFT_1MPa(39075:40678),(-GammaFT_1(39075:40678)+GammaFT_1(39075)))
plot(PressureFT_1MPa(16706:17500),(-GammaFT_1(16706:17500)+GammaFT_1(16706)))
plot(PressureFT_2MPa(38974:40025),(-GammaFT_2(38974:40025)+GammaFT_2(38974)))

% xlim([0 200])
% ylim([-2 40]);set(gca,'YTick',[0:10:40]); 
legend('\gamma_{z\theta} Model','\gamma_{z\theta} Experimental SAMPLE 1','\gamma_{z\theta} Experimental SAMPLE 2')
 xlabel('Pressure, (MPa)')
 ylabel('Shear strain, \gamma_{z\theta}')


 grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 
% print(gcf,'optimalalpha0.01deg.png','-dpng','-r700');


ThetaExperimental = [2.5 4.5 7 10 15 18.5 22.5 27 32.5 35];
GammExperiemntal = ThetaExperimental*0.185;

