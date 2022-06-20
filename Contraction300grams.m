%Cavtappi Contraction 0.2 kg 


clear all 
close all
clc

L1 = (86+19)*0.9; %mm
L2 = (52+15)*0.9; %mm
L3 = (115+25)*0.9; %mm

% Retrieve Data Sample 1
    filename = '300grams Sample 1 Test 1.xlsx';
    sheet = 'Sheet1';
    
    FT_1 = xlsread(filename, sheet,'C:E');
        time_1 = FT_1(:,1);  
        Pressure_1 = FT_1(:,2);
        Disp_1 = FT_1(:,3);
        time_1_p = time_1;
    
        Strain1 = (Disp_1/L1)*100;
      
% Retrieve Data Sample 2
    filename = '300grams Sample 2 Test 1.xlsx';
    sheet = 'Sheet1';
    
    FT_2 = xlsread(filename, sheet,'C:E');
        time_2 = FT_2(:,1);  
        Pressure_2 = FT_2(:,2);
        Disp_2 = FT_2(:,3);
        time_2_p = time_1;
    
        Strain2 = (Disp_2/L2)*100;
        
% Retrieve Data Sample 3
    filename = '300grams Sample 3 Test 1.xlsx';
    sheet = 'Sheet1';
    
    FT_3 = xlsread(filename, sheet,'C:E');
        time_3 = FT_3(:,1);  
        Pressure_3 = FT_3(:,2);
        Disp_3 = FT_3(:,3);
        time_3_p = time_3;
    
        Strain3 = (Disp_3/L3)*100;    
        
        
fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);

yyaxis right;
H1 = plot(time_3,Pressure_3,':','LineWidth',1,'Color',[0.5,0.5,0.5]); 
ylabel('Pressure, (MPa)')
set(gca,'ycolor',[0.5,0.5,0.5])
% ylim([0 60])
yyaxis left;
H2 =  plot(time_3,Strain3-Strain3(1),'Linewidth',1,'Color',[0, 0, 0]);ylabel('Act. Strain')
set(gca,'ycolor',[0, 0, 0])
xlabel('Time (s)')
% xlim([0,8]);
% set(gca,'XTick',[0:1:8]);
% title('400 grams')
% legend('SAMPLE 1','SAMPLE 1','Location','Northwest')
grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color        
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = 3 %N
Pi = [0:0.01:1.9]; %MPa which is equal to 200 psi\
P0 = 0; %Mpa

r0 = (2/2); %mm
ri = (0.8/2); %mm
r = ri:0.0001:r0; %mm
d0 = r0*2; %mm
di =ri*2; %mm
R = 1.5; %mm
a = 55; %deg

%MATERIAL PROPERTIES OF THE STRAIGHT DRAWN TUBE

u_12 = 0.19;
u_23 = 0.6;
E2 = 9.1*1; %MPa
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
% u_r = (((r(1)^2)*(Pi(i))*r)./(E2*(r(end)^2-r(1)^2))).*((1-u_23)+(1+u_23)*(r(end)^2./r.^2));
u_r = (((ri^2)*(Pi(i))*r)./(E2*(r0^2-ri^2))).*((1-u_23)+(1+u_23)*(r0^2./r.^2));

% r = (ri+u_r(1)):0.025:(r0+u_r(end)); %mm

% Sigma_z = (Pi(i)*(ri)^2-P0*(r0)^2)/((r0)^2-(ri)^2); 
% Term2 =((Pi(i)-P0)*(r0)^2*(ri)^2)./(((r0)^2-(ri)^2)*(r).^2);
% Sigma_theta = Sigma_z + Term2;
% Sigma_rad = Sigma_z - Term2;

Sigma_z = (Pi(i)*(ri+u_r(1))^2-P0*(r0+u_r(end))^2)/((r0+u_r(end))^2-(ri+u_r(1))^2); 
Term2 =((Pi(i)-P0)*(r0+u_r(end))^2*(ri+u_r(1))^2)./(((r0+u_r(end))^2-(ri+u_r(1))^2)*(u_r+r).^2);
Sigma_theta = Sigma_z + Term2;
Sigma_rad = Sigma_z - Term2;

% Sigma_z = (Pi(i)*(r(1))^2-P0*(r(end))^2)/((r(end))^2-(r(1))^2); 
% Term2 =((Pi(i)-P0)*(r(end))^2*(r(1))^2)./(((r(end))^2-(r(1))^2)*(r).^2);
% Sigma_theta = Sigma_z + Term2;
% Sigma_rad = Sigma_z - Term2;

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

x1 = (r0+u_r(end))-((2*(r0+u_r(end)))^2/(16*R));
x2 = -(r0+u_r(end))-((2*(r0+u_r(end)))^2/(16*R));
x = x2:0.000005:x1;

tau_torsion = (F*R^2*x)./(J*(R-((2*(r0+u_r(end)))^2/(16*R))-x));
tau_direct = F/A;

tau = tau_torsion+tau_direct;

a = 0;
b = 0;
th = 0:pi/50:2*pi;
xunit0 = (r0+u_r(end))*cos(th) + a;
yunit0 = (r0+u_r(end)) * sin(th) + b;
xuniti = (ri+u_r(1))*cos(th) + a;
yuniti = (ri+u_r(1)) * sin(th) + b;
r_graph = -(r0+u_r(end)):0.000005:(r0+u_r(end));

%Cleaning up Shear Stresses at the hollowed section
for j = 1:length(r_graph)  
    if -(ri+u_r(1))<r_graph(j)&& r_graph(j)<(ri+u_r(1))
        tau(j) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Plot CHECK%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h = figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);
% axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
% 
% yyaxis right;
% H1 = plot(xunit0,yunit0,xuniti,yuniti,'-','LineWidth',1,'Color',[0.5,0.5,0.5]);
% ylim([-1.5 1.5])
% xlim([-1.5 1.5])
% ylabel('Length, (mm)')
% set(gca,'ycolor',[0.5,0.5,0.5])
% yyaxis left;
% H3 =  plot(r_graph,tau,':','Linewidth',1.8,'Color',[0, 0, 0]);
% ylabel('Shear Stress,\tau (MPa)')
% ylim([-6 6])
% set(gca,'ycolor',[0, 0, 0])
% xlabel('Length, (mm)')
% grid on 
% set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color
% drawnow 
%       % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if i == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end 
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ri = ri+u_r(1)
% r0 = r0+u_r(end)
% r = ri:0.001:r0
 r = ri+u_r(1):0.0001:r0+u_r(end);
%COMPLIANCE MATERIX. SOLVING STRAINS.
tau_in = interp1(r_graph(1:end),tau,r);
tau_out = abs(interp1(r_graph(1:end),tau,-(r)));
tau_avg = (abs(tau_in)+abs(tau_out))/2;



%%%%%THESE ARE THE INNER AND OUTER GAMMAS USING THE AVERAGE SHEAR STRESS
%%%%%PROFILE%%%%%%%%%%%%%%%%%%%%
% ri = (ri+u_r(1));
% r0 = (r0+u_r(end));
% gammaZ_inner(i) = S_bar(1,6)*Sigma_z;
% gammaTheta_inner(i) = S_bar(2,6)*Sigma_theta(10);
% gammaR_inner(i) = S_bar(3,6)*Sigma_rad(10);
% gammaTau_inner(i) = S_bar(6,6)*tau_avg(10);%-S_bar(6,6)*Tau(1);
% gammaTau_SUB_inner(i) = gammaTau_inner(i)-gammaTau_inner(1);
% 
% Gamma_z_theta_inner(i) = gammaZ_inner(i) + gammaTheta_inner(i) + gammaR_inner(i) + gammaTau_SUB_inner(i);
% 
% gammaZ_outer(i) = S_bar(1,6)*Sigma_z;
% gammaTheta_outer(i) = S_bar(2,6)*Sigma_theta(end);
% gammaR_outer(i) = S_bar(3,6)*Sigma_rad(end);
% gammaTau_outer(i) = S_bar(6,6)*tau_avg(end);%-S_bar(6,6)*Tau(1);
% gammaTau_SUB_outer(i) = gammaTau_outer(i)-gammaTau_outer(1);
% 
% Gamma_z_theta_outer(i) = gammaZ_outer(i) + gammaTheta_outer(i) + gammaR_outer(i) + gammaTau_SUB_outer(i);

%%%%%THESE ARE THE GAMMAS AT THE FOR CRITICAL LOCATIONS, gamma_ii,
%%%%%gamma_io, gamma_oi, gamma_oo.
gammaZ_oo(i) = S_bar(1,6)*Sigma_z(1);
gammaTheta_oo(i) = S_bar(2,6)*Sigma_theta(end);
gammaR_oo(i) = S_bar(3,6)*Sigma_rad(end);
gammaTau_oo(i) = S_bar(6,6)*tau_out(end);%-S_bar(6,6)*Tau(1);
gammaTau_SUB_oo(i) = gammaTau_oo(i)-gammaTau_oo(1);
Gamma_z_theta_oo(i) = gammaZ_oo(i) + gammaTheta_oo(i) + gammaR_oo(i) - abs(gammaTau_SUB_oo(i));

gammaZ_io(i) = S_bar(1,6)*Sigma_z(1);
gammaTheta_io(i) = S_bar(2,6)*Sigma_theta(2);
gammaR_io(i) = S_bar(3,6)*Sigma_rad(2);
gammaTau_io(i) = S_bar(6,6)*tau_out(2);%-S_bar(6,6)*Tau(1);
gammaTau_SUB_io(i) = gammaTau_io(i)-gammaTau_io(1);
Gamma_z_theta_io(i) = gammaZ_io(i) + gammaTheta_io(i) + gammaR_io(i) - abs(gammaTau_SUB_io(i));

gammaZ_oi(i) = S_bar(1,6)*Sigma_z(1);
gammaTheta_oi(i) = S_bar(2,6)*Sigma_theta(end);
gammaR_oi(i) = S_bar(3,6)*Sigma_rad(end);
gammaTau_oi(i) = S_bar(6,6)*tau_in(end);%-S_bar(6,6)*Tau(1);
gammaTau_SUB_oi(i) = gammaTau_oi(i)-gammaTau_oi(1);
Gamma_z_theta_oi(i) = gammaZ_oi(i) + gammaTheta_oi(i) + gammaR_oi(i) - abs(gammaTau_SUB_oi(i));

gammaZ_ii(i) = S_bar(1,6)*Sigma_z(1);
gammaTheta_ii(i) = S_bar(2,6)*Sigma_theta(4);
gammaR_ii(i) = S_bar(3,6)*Sigma_rad(4);
gammaTau_ii(i) = S_bar(6,6)*tau_in(4);%-S_bar(6,6)*Tau(1);
gammaTau_SUB_ii(i) = gammaTau_ii(i)-gammaTau_ii(1);
Gamma_z_theta_ii(i) = gammaZ_ii(i) + gammaTheta_ii(i) + gammaR_ii(i) - abs(gammaTau_SUB_ii(i));

gammaZ = S_bar(1,6)*Sigma_z(1);
gammaTheta = S_bar(2,6)*Sigma_theta;
gammaR = S_bar(3,6)*Sigma_rad;
gammaTau = S_bar(6,6)*tau_avg;%-S_bar(6,6)*Tau(1);

gammaZ_avg(i) = gammaZ;
gammaTheta_avg(i) = mean(gammaTheta(1:end));
gammaR_avg(i) = mean(gammaR(1:end));
gammaTau_avg(i) = mean(gammaTau(1:end));
gammaTau_avgSUB(i) = gammaTau_avg(i)-gammaTau_avg(1);

Gamma_z_theta_avg(i) = gammaZ_avg(i) + gammaTheta_avg(i) + gammaR_avg(i) - abs(gammaTau_avgSUB(i));

Delta(i) = (2*R^2*23*Gamma_z_theta_oi(i)*pi)/(r0+u_r(end));
Strain(i) = (Delta(i)/96)*100; %for .3 kg load


% Gamma_z_theta_avg(i) = Gamma_z_theta;
end


fig=figure('units','inch','position',[0,0,3.6,2.9]); hold on; grid on; set(gca,'FontSize',10);
% plot((Pressure_1(51017-70:50:55236-70)-Pressure_1(51017-70))/145.038,-(Strain1(51017:50:55236)-Strain1(51017)),'Linewidth',1,'Color',[0.5,0.5,0.5])
plot((Pressure_2(36005-30:60:38401-30)-Pressure_2(36005-30))/145.038,-(Strain2(36005:60:38401)-Strain2(36005)),'-','Linewidth',1,'Color',[0 0.4 0.9])
plot((Pressure_3(32901-50:60:38763-50)-Pressure_3(32901-50))/145.038,-(Strain3(32901:60:38763)-Strain3(32901)),'-','Linewidth',1,'Color',[0.2 0.8 0.8])
plot(Pi(1:end),Strain(1:end),':','Linewidth',2.3,'Color',[0,0,0])

% plot(145.038*PressureFT_2MPa(14587:60:19215)-145.038*PressureFT_2MPa(14587),(-GammaFT_2(14587:60:19215)+GammaFT_2(14587)),':','Linewidth',1.5)
% plot(145.038*PressureFT_2MPa(27467:60:33976)-145.038*PressureFT_2MPa(27467),(-GammaFT_2(27467:60:33976)+GammaFT_2(27467)),':','Linewidth',1.5)
  ylim([-5 35])
% -Gamma_z_theta(1)
  xlim([0 2])
% ylim([-2 40]);set(gca,'YTick',[0:10:40]); 
%  legend('Sample 1','Sample 2','Model')
 xlabel('Pressure, (MPa)')
 ylabel('Act. Strain, {\it ɛ} (%)')
%  ylim([-0.01 0.08])

grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 
print(gcf,'CavatappiArtificialMuscle_300grams.png','-dpng','-r700');   

fig=figure('units','inch','position',[0,0,3.6,2.9]); hold on; grid on; set(gca,'FontSize',10);
% plot((Pressure_1(51017-70:50:55236-70)-Pressure_1(51017-70))/145.038,-(Strain1(51017:50:55236)-Strain1(51017)),'Linewidth',1,'Color',[0.5,0.5,0.5])
plot((Pressure_2(36005-30:60:38401-30-1350)-Pressure_2(36005-30))/145.038,-(Strain2(36005:60:38401-1350)-Strain2(36005)),'-','Linewidth',1,'Color',[0 0.4 0.9])
plot((Pressure_3(32901-50:60:38763-50-3700)-Pressure_3(32901-50))/145.038,-(Strain3(32901:60:38763-3700)-Strain3(32901)),'-','Linewidth',1,'Color',[0.2 0.8 0.8])
plot(Pi(1:end),Strain(1:end),':','Linewidth',2.3,'Color',[0,0,0])

% plot(145.038*PressureFT_2MPa(14587:60:19215)-145.038*PressureFT_2MPa(14587),(-GammaFT_2(14587:60:19215)+GammaFT_2(14587)),':','Linewidth',1.5)
% plot(145.038*PressureFT_2MPa(27467:60:33976)-145.038*PressureFT_2MPa(27467),(-GammaFT_2(27467:60:33976)+GammaFT_2(27467)),':','Linewidth',1.5)
  ylim([-5 35])
% -Gamma_z_theta(1)
  xlim([0 2])
% ylim([-2 40]);set(gca,'YTick',[0:10:40]); 
%  legend('Sample 1','Sample 2','Model')
 xlabel('Pressure, (MPa)')
 ylabel('Act. Strain, {\it ɛ} (%)')
%  ylim([-0.01 0.08])

grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 
print(gcf,'CavatappiArtificialMuscle_300gramsLOADING.png','-dpng','-r700');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%                    TEMPORAL PLOT                                %%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p_temporal = (Pressure_3(27366:39098)-Pressure_3(27366))/145.038;
% po = 0; %Mpa
% P = 1.4;
% r = (2/2)/1000; %m
% ri = (0.85/2)/1000; %m
% r_ave = (1.45/2)/1000;
% t = 0.0005; %m 
% r_spool = 3; %mm Spool radius
% a = 60; %deg
% 
% u_12 = 0.19;
% u_23 = 0.85;
% %  E1 = 50; %MPa
% E2 = 4.5; %MPa
% G12 = 6; %MPa
% G23 = 9; %MPa
% E1 = 31.72; %MPa
% % E2 = 2.09; %MPa
% % G12 = 2.78; %MPa
% % G12 = E1/(2*(1+u_12)); %MPa
% u_21 = (E2*u_12)/E1;
% 
% S = [1/E1 -(u_12/E1) -(u_12/E1) 0 0 0;
%     -(u_12/E1) 1/E2 -(u_23/E2)  0 0 0;
%     -(u_12/E1) -(u_23/E2) 1/E2 0 0 0;
%     0 0 0 1/G23 0 0;
%     0 0 0 0 1/G12 0;
%     0 0 0 0 0 1/G12];
% 
% S_bar(1,1) = (cosd(a)^4)*S(1,1) + (cosd(a)^2)*(sind(a)^2)*(2*S(1,2)*S(6,6)) + (sind(a)^4)*S(2,2);
% S_bar(1,2) = (cosd(a)^2)*(sind(a)^2)*(S(1,1)+S(2,2)-S(6,6)) + ((cosd(a)^4)+(sind(a)^4))*S(1,2);
% S_bar(1,3) = (cosd(a)^2)*S(1,3) + (sind(a)^2)*S(2,3); 
% S_bar(1,6) = (cosd(a))*(sind(a))*((cosd(a)^2)-(sind(a)^2))*(2*S(1,2)+S(6,6)) + 2*(cosd(a))*(sind(a))*((sind(a)^2)*S(2,2)-(cosd(a)^2)*S(1,1));
% S_bar(2,2) = (cosd(a)^4)*S(2,2) + (cosd(a)^2)*(sind(a)^2)*(2*S(1,2)*S(6,6)) + (sind(a)^4)*S(1,1);
% S_bar(3,3) = S(3,3);
% S_bar(2,3) = (cosd(a)^2)*S(2,3) + (sind(a)^2)*S(1,3); 
% S_bar(2,6) = (cosd(a))*(sind(a))*((sind(a)^2)-(cosd(a)^2))*(2*S(1,2)+S(6,6)) + 2*(cosd(a))*(sind(a))*((cosd(a)^2)*S(2,2)-(sind(a)^2)*S(1,1));
% S_bar(6,6) = 4*(cosd(a)^2)*(sind(a)^2)*(S(1,1)+S(2,2)-2*S(1,2)) + (((cosd(a)^2)-(sind(a)^2))^2)*S(6,6);
% S_bar(3,6) = 2*(cosd(a))*(sind(a))*(S(2,3)-S(1,3));
%  
% %Epsilon_z is assumed to be 0.
% Sigma_z = (p_temporal*ri^2-po*r^2)/(r^2-ri^2);
% 
% % Sigma_z = ((p*ri^2)-(po*r^2))/(r^2-ri^2);
% Term2 =((p_temporal-po)*r^2*ri^2)/((r^2-ri^2)*r_ave^2);
% Sigma_theta = Sigma_z + Term2;
% Sigma_rad = Sigma_z - Term2;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F_iso = 1.4; % N. Load of the mass
% 
% u_r = (2*((ri^2)*(p_temporal)*r)/(E2*(r^2-ri^2)));
% u_ri = (((ri^3)*(p_temporal))/(E2*(r^2-ri^2)))*((1-u_23)+(1+u_23)*(r^2/ri^2));
% u_ave = (((ri^2)*(p_temporal)*r_ave)/(E2*(r^2-ri^2)))*((1-u_23)+(1+u_23)*(r^2/r_ave^2))
% 
%  J = (pi*(((r+u_r).^4)-((ri+u_ri).^4)))/2;
% % J = 1.85806490894224e-12;
% Tau_z_theta = zeros(length(p),1);
% Tau_z_theta = -((F_iso*(r_spool/1000)*(r_ave+u_ave))./J)/10^6;
% 
% % Tau_z_theta = 0;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % R = [(0.85/2)/1000:0.00005:(1.85/2)/1000];
% % Sigma_z_R = (P*ri^2)/(r^2-ri^2);
% % Term2_R =(P*r^2*ri^2)./((R.^2)*(r^2-ri^2));
% % Sigma_theta_R = Sigma_z_R + Term2_R;
% % Sigma_rad_R = Sigma_z_R - Term2_R;
% % Tau_z_theta_R = 0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% gammaTheta_temporal = S_bar(2,6)*Sigma_theta
% gammaR_temporal = S_bar(3,6)*Sigma_rad
% gammaTau_temporal = S_bar(6,6)*Tau_z_theta-S_bar(6,6)*Tau_z_theta(1)
% 
% Epsilon_z = S_bar(1,1)*Sigma_z + S_bar(1,2)*Sigma_theta + S_bar(1,3)*Sigma_rad + S_bar(1,6)*Tau_z_theta;
% 
% Epsilon_theta = S_bar(1,2)*Sigma_z + S_bar(2,2)*Sigma_theta + S_bar(2,3)*Sigma_rad + S_bar(2,6)*Tau_z_theta;
% 
% Epsilon_rad = S_bar(1,3)*Sigma_z + S_bar(2,3)*Sigma_theta + S_bar(3,3)*Sigma_rad + S_bar(3,6)*Tau_z_theta;
% 
% Gamma_z_theta_temporal = S_bar(1,6)*Sigma_z + S_bar(2,6)*Sigma_theta + S_bar(3,6)*Sigma_rad - gammaTau_temporal;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% n_c = 20; %# of inserted coils
% D_mean = 2.85;
% Z = n_c*D_mean*pi*Gamma_z_theta_temporal
% Strain_model_temporal = (Z/94.5)*100;
% 
% fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);
% 
% subplot(2,1,1);
% plot(time_3(27366:39098)-time_3(27366),p_temporal,':','LineWidth',1.5,'Color',[0,0,0]); 
% ylabel('Pressure, (MPa)')
% xlabel('Time (s)')
% 
%  ylim([-0.05 2.2])
% 
% subplot(2,1,2); 
% plot(time_3(27366:39098)-time_3(27366),-(Strain3(27366:39098)-Strain3(27366)),':','Linewidth',1.5,'Color',[0.5,0.5,0.5]); hold on
% plot(time_3(27366:39098)-time_3(27366),(Strain_model_temporal),'Linewidth',0.8,'Color',[0,0,0]);
% ylabel('Act. Strain, ɛ (%)')
% xlabel('Time (s)')
%  legend('Model','Sample 1')
% 
% set(gca,'ycolor',[0, 0, 0])
% xlabel('Time (s)')
% % xlim([0,8]);
% % set(gca,'XTick',[0:1:8]);
% % title('400 grams')
% % legend('SAMPLE 1','SAMPLE 1','Location','Northwest')
% grid on 
% set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color  
% % 
% % fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',8);
% % 
% % yyaxis right;
% % H1 = plot(time_3(27366:39098)-time_3(27366),p_temporal,':','LineWidth',1,'Color',[0.5,0.5,0.5]); 
% % ylabel('Pressure, (MPa)')
% % set(gca,'ycolor',[0.5,0.5,0.5])
% %  ylim([-0.05 2.2])
% % yyaxis left;
% % H2 =  plot(time_3(27366:39098)-time_3(27366),-(Strain3(27366:39098)-Strain3(27366)),'Linewidth',1,'Color',[0, 0, 0]);ylabel('Act. Strain')
% % H3 =  plot(time_3(27366:39098)-time_3(27366),(Strain_model_temporal),'Linewidth',1,'Color',[1, 0, 1]);
% % ylabel('Act. Strain')
% % 
% % set(gca,'ycolor',[0, 0, 0])
% % xlabel('Time (s)')
% % % xlim([0,8]);
% % % set(gca,'XTick',[0:1:8]);
% % % title('400 grams')
% % % legend('SAMPLE 1','SAMPLE 1','Location','Northwest')
% % grid on 
% % set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color  
% % 
% % 
