close all
clc
clear all

F = 2 %N
Pi = [0:0.01:2]; %MPa 
P0 = 0; %Mpa

r0 = (2/2); %mm
ri = (0.8/2); %mm
r = ri:0.0001:r0; %mm
d0 = r0*2; %mm
di =ri*2; %mm
R = 1.5; %mm
a = [0:1:90]; %deg

%MATERIAL PROPERTIES OF THE STRAIGHT DRAWN TUBE

u_12 = 0.19;
u_23 = 0.6;
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

  for i = 1:length(a)    
    for j = 1:length(Pi)

        j
        %%%%%THICK WALL VESSEL. THEORY%%%
        r = ri:0.001:r0; %mm
        u_r = (((ri^2)*(Pi(j))*r)./(E2*(r0^2-ri^2))).*((1-u_23)+(1+u_23)*(r0^2./r.^2));

        Sigma_z(j) = (Pi(j)*(ri+u_r(1))^2-P0*(r0+u_r(end))^2)/((r0+u_r(end))^2-(ri+u_r(1))^2); 
        Term2(j) =((Pi(j)-P0)*(r0+u_r(end))^2*(ri+u_r(1))^2)./(((r0+u_r(end))^2-(ri+u_r(1))^2)*(u_r(end)+r0).^2);
        Sigma_theta(j) = Sigma_z(j) + Term2(j);
        Sigma_rad(j) = Sigma_z(j) - Term2(j);

        J(j) = (pi*(((r0+u_r(end)).^4)-((ri+u_r(1)).^4)))/2;
        A(j) = (pi*(((r0+u_r(end)).^2)-((ri+u_r(1)).^2)));

        %SHEAR STRESS PROFILE

        x1 = (r0+u_r(end))-((2*(r0+u_r(end)))^2/(16*R));
        x2 = -(r0+u_r(end))-((2*(r0+u_r(end)))^2/(16*R));
%         x = x2:0.000005:x1;

        tau_torsion(j) = (F*R^2*x1)./(J(j)*(R-((2*(r0+u_r(end)))^2/(16*R))-x1));
        tau_direct(j) = F/A(j);

        tau_in(j) = tau_torsion(j)+tau_direct(j);


%         r_graph = -(r0+u_r(end)):0.000005:(r0+u_r(end));

            %Cleaning up Shear Stresses at the hollowed section
%         for k = 1:length(r_graph)  
%             if -(ri+u_r(1))<r_graph(k)&& r_graph(k)<(ri+u_r(1))
%                 tau(k) = 0;
%             end
%         end
%         r = ri+u_r(1):0.0001:r0+u_r(end);

        %COMPLIANCE MATERIX. SOLVING STRAINS.
%         tau_in = interp1(r_graph(1:end),tau,r);
        
    i
    S_bar(1,1) = (cosd(a(i)).^4)*S(1,1) + (cosd(a(i)).^2)*(sind(a(i)).^2)*(2*S(1,2)*S(6,6)) + (sind(a(i)).^4)*S(2,2);
    S_bar(1,2) = (cosd(a(i)).^2)*(sind(a(i)).^2)*(S(1,1)+S(2,2)-S(6,6)) + ((cosd(a(i)).^4)+(sind(a(i)).^4))*S(1,2);
    S_bar(1,3) = (cosd(a(i)).^2)*S(1,3) + (sind(a(i)).^2)*S(2,3); 
    S_bar(1,6) = (cosd(a(i)))*(sind(a(i)))*((cosd(a(i)).^2)-(sind(a(i)).^2))*(2*S(1,2)+S(6,6)) + ...
    2*(cosd(a(i)))*(sind(a(i)))*((sind(a(i)).^2)*S(2,2)-(cosd(a(i)).^2)*S(1,1));
    S_bar(2,2) = (cosd(a(i)).^4)*S(2,2) + (cosd(a(i)).^2)*(sind(a(i)).^2)*(2*S(1,2)*S(6,6)) + (sind(a(i)).^4)*S(1,1);
    S_bar(3,3) = S(3,3);
    S_bar(2,3) = (cosd(a(i)).^2)*S(2,3) + (sind(a(i)).^2)*S(1,3); 
    S_bar(2,6) = (cosd(a(i)))*(sind(a(i)))*((sind(a(i)).^2)-(cosd(a(i)).^2))*(2*S(1,2)+S(6,6)) + 2*(cosd(a(i)))*(sind(a(i)))*((cosd(a(i)).^2)*S(2,2)-(sind(a(i)).^2)*S(1,1));
    S_bar(6,6) = 4*(cosd(a(i)).^2)*(sind(a(i)).^2)*(S(1,1)+S(2,2)-2*S(1,2)) + (((cosd(a(i)).^2)-(sind(a(i))^2)).^2)*S(6,6);
    S_bar(3,6) = 2*(cosd(a(i)))*(sind(a(i)))*(S(2,3)-S(1,3));
    
%%%%%THESE ARE THE GAMMAS AT THE FOR CRITICAL LOCATIONS, gamma_ii,
%%%%%gamma_io, gamma_oi, gamma_oo.

    gammaZ_oi(i,j) = S_bar(1,6)*Sigma_z(j);
    gammaTheta_oi(i,j) = S_bar(2,6)*Sigma_theta(j);
    gammaR_oi(i,j) = S_bar(3,6)*Sigma_rad(j);
    gammaTau_oi(i,j) = S_bar(6,6)*tau_in(j);%-S_bar(6,6)*Tau(1);
    gammaTau_SUB_oi(i,j) = gammaTau_oi(i,j)-S_bar(6,6)*tau_in(1);
    Gamma_z_theta_oi(i,j) = gammaZ_oi(i,j) + gammaTheta_oi(i,j) + gammaR_oi(i,j) - abs(gammaTau_SUB_oi(i,j));

% %     %2 Newtons
    Delta(i,j) = (2*R^2*23*Gamma_z_theta_oi(i,j)*pi)/(r0+u_r(end));
    Strain(i,j) = (Delta(i,j)/80)*100; %for .2 kg load
    
%     %3 Newtons
%     Delta(i,j) = (2*R^2*23*Gamma_z_theta_oi(i,j)*pi)/(r0+u_r(end));
%     Strain(i,j) = (Delta(i,j)/96)*100; %for .3 kg load
%     
    %4 Newtons
%     Delta(i,j) = (2*R^2*23*Gamma_z_theta_oi(i,j)*pi)/(r0+u_r(end));
%     Strain(i,j) = (Delta(i,j)/110)*100; %for .4 kg load
    end
 
end

fig=figure('units','inch','position',[0,0,3.6*1.5,2.9*1.5]); hold on; grid on; set(gca,'FontSize',10);
s = surf(a,Pi,Strain')
% contour3(a,Pi,Strain',[25 20 15 10 5 0 -5],'-','Linewidth',1,'Color',[0.1 0.7 0.7])
 s.EdgeColor = 'none';
colormap('gray')
 colorbar
% plot(a,Strain,'-','Linewidth',1,'Color',[0.1 0.7 0.7])
% 
% 
 xlabel('Fiber Angle, \alpha ( ^\circ )'); xlim([0 90])
 ylabel('Pressure (MPa)'); ylim([0 2])  
 zlabel('Act. Strain, {\it ɛ} (%)'); zlim([-15 50])
% 
% grid on 
% set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color
print(gcf,'CavatappiArtificialMuscle_AlphaOptimal.png','-dpng','-r700');   


fig2 = figure('units','inch','position',[0,0,3.6,2.9]); hold on; grid on; set(gca,'FontSize',10);
  %plot(a,Strain(:,101),'-','Linewidth',1.7,'Color',[0.1 0.7 0.5])
plot(a,Strain(:,101),'-','Linewidth',1.7,'Color',[0.1 0.8 0.5])
 %plot(a,Strain(:,101),'-','Linewidth',1.7,'Color',[0.1 0.4 0.1])
 xlabel('Fiber Angle, \alpha ( ^\circ )'); xlim([0 90])
 ylabel('Act. Strain, {\it ɛ} (%)'); ylim([-8 12])
% 
grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color
print(gcf,'CavatappiArtificialMuscle_AlphaOptimalLOADS.png','-dpng','-r700');   



fig3 = figure('units','inch','position',[0,0,3.6,2.9]); hold on; grid on; set(gca,'FontSize',10);

 xlabel('Pressure (MPa)'); xlim([0 2])
 ylabel('Act. Strain, {\it ɛ} (%)'); ylim([-8 2])
% 
grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color
print(gcf,'CavatappiArtificialMuscle_AlphaOptimalLOADS.png','-dpng','-r700');   


