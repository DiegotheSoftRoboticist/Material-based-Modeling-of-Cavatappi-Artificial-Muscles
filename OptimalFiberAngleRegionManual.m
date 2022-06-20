
close all 
clc
clear all

F = 0.3 %N
Pi = [0:0.01:2]; %MPa which is equal to 200 psi\
P0 = 0; %Mpa

r0 = (2/2); %mm
ri = (0.8/2); %mm
r = ri:0.0001:r0; %mm
d0 = r0*2; %mm
di =ri*2; %mm
R = 1.5; %mm
a = 0; %deg

%MATERIAL PROPERTIES OF THE STRAIGHT DRAWN TUBE

u_12 = 0.19;
u_23 = 0.55;
E2 = 9; %MPa
G12 = 10; %MPa
 G23 = 0; %MPa
E1 = 60; %MPa
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
Strain(i) = (Delta(i)/110)*100; %for .4 kg load

% Gamma_z_theta_avg(i) = Gamma_z_theta;
end



fig=figure('units','inch','position',[0,0,3.6,2.9]); hold on; grid on; set(gca,'FontSize',10);
  plot(Pi(1:end),Strain(1:end),'-','Linewidth',1.2,'Color',[0 0 0])
%  plot(Pi(1:end),Strain(1:end),':','Linewidth',1.7,'Color',[0.1 0.6 0.7])
%plot(Pi(1:end),Strain(1:end),'-.','Linewidth',1.7,'Color',[0 0.4470 0.7410])
  ylim([-8 2])
% -Gamma_z_theta(1)
  xlim([0 2])
% ylim([-2 40]);set(gca,'YTick',[0:10:40]); 
%  legend('Sample 1','Sample 2','Model')
 xlabel('Pressure, (MPa)')
 ylabel('Act. Strain, {\it ɛ} (%)')
%  ylim([-0.01 0.08])

grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color 
print(gcf,'CavatappiArtificialMuscle_PolarMomentActuation.png','-dpng','-r700');   


