
poisson_12 = [0.195, 0.218, 0.201];
poisson_23 = [0.579, 0.507, 0.618];
label = [1,2,3];

mean_12 = mean(poisson_12)
mean_23 = mean(poisson_23)
mean_12 = mean_12*ones(3,1)
mean_23 = mean_23*ones(3,1)

v=500*ones(360,1)

fig=figure('units','inch','position',[0,0,3.5,2.5]); hold on; grid on; set(gca,'FontSize',10);
 
% plot(StrainG12_Stat,Stress_mean,'LineWidth',1.5,'Color',[0.75,0.75,0.75]);


plot(label,poisson_12,'o','MarkerSize',2,...
    'MarkerEdgeColor',[0.75,0.75,0.75],...
    'MarkerFaceColor',[0, 0.4470, 0.7410],'LineWidth',0.5,'Color',[0, 0.4470, 0.7410])

plot(label,poisson_23,'s','MarkerSize',2,...
    'MarkerEdgeColor',[0.75,0.75,0.75],...
    'MarkerFaceColor',[0, 0.75, 0.75],'LineWidth',0.5,'Color',[0, 0.75, 0.75])

plot(label,mean_12,'LineWidth',0.85,'Color',[0.8,0,0]);
plot(label,mean_23,'LineWidth',0.85,'Color',[0.8,0,0]);

 ylim([0 1]); set(gca,'YTick',[0:0.2:1]);
%  
 xlim([0.75 3.25]); set(gca,'xTick',[1:3]);


%  legend([a b c d e],'m = 0.1 kg','m = 0.2 kg','m = 0.3 kg','m = 0.4 kg','m = 0.5 kg')



       ylabel('Poisson Ratios, \nu')
%        xlabel('Sample, \gamma')
 
grid on 
set(gca,'GridColor',[0.1 0.2 0.9]) % a bluish color
print(gcf,'Poissons12and23.png','-dpng','-r700');   
