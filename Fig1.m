% Code to create Figure 1 of NCLIM-20051128B, "A Mechanism for Regional
% Variations in Snowpack Melt Under Rising Temperature" by Evan and
% Eisenman

% Load SNOTEL data. This can be acquired from:
% https://doi.pangaea.de/10.1594/PANGAEA.896396
swe = swe(:,d365,:);

g = ~isnan(mean(nanmean(swe,2),3)) & lat<50;
swe = swe(g,:,:);
sitenum = sitenum(g);
state = cellstr(state(g,:));
lat = lat(g);
lon = lon(g);

% SNOTEL Site: Virginia Lakes Ridge 
% State: California 
% Site Number: 846 s
% County: Mono 
% Latitude: 38 deg; 4 min N 
% Longitude: 119 deg; 14 min W 
% Elevation: 9400 feet 
% Reporting since: 1978-10-01 

y1 = swe(30,:,years==2017)';
y2 = swe(30,:,years==2014)';
x = (1:365)';

p1 = plot(x,y1,'linewidth',2);
grid on
box on
set(gca,'xlim',[30 270],'xtick',30:30:270)
set(gca,'ylim',[0 120],'ytick',0:15:120)
xlabel('$t$ (water year day)','interpreter','Latex')
ylabel('$S$ (cm)','interpreter','Latex')
set(gca,'fontsize',12,'TickLabelInterpreter','Latex')

m1 = round(mean(find(y1==max(y1))));
f1 = find(y1<0.1 & x<m1,1,'last');
z1 = find(x>m1 & y1<0.1,1);
Rc1 = max(y1)/(m1-f1);
Rm1 = max(y1)/(z1-m1);

hold on
p1c = plot([f1 m1],[0 max(y1)],'--','Color',p1.Color);
p1m = plot([m1 z1],[max(y1) 0],'-.','Color',p1.Color);
hold off


hold on
p2 = plot(x,y2,'r','linewidth',2);
m2 = round(mean(find(y2==max(y2))));
f2 = find(y2<0.1 & x<m2,1,'last');
z2 = find(x>m2 & y2<0.1,1);
Rc2 = max(y2)/(m2-f2);
Rm2 = max(y2)/(z2-m2);
p2c = plot([f2 m2],[0 max(y2)],'--','Color',p2.Color);
p2m = plot([m2 z2],[max(y2) 0],'-.','Color',p2.Color);
hold off

ttl = title( ...
    'Virginia Lakes Ridge SNOTEL Station (38$^\circ$N 119$^\circ$W)', ...
    'Interpreter','Latex','FontSize',12);


t1 = text(195.5,31,'$\zeta_{2017}=260$','Interpreter','Latex', ...
    'FontSize',12,'BackgroundColor','w');
hold on
q1 = quiver(211,27,46,-26,1,'k'); %,'MaxHeadSize',.5);

t2 = text(156.5,11,'$\zeta_{2014}=226$','Interpreter','Latex', ...
    'FontSize',12,'BackgroundColor','w');
q2 = quiver(176,7,46,-6,1,'k'); %,'MaxHeadSize',.5);

hold off


lg = legend([p1 p1c p1m p2 p2c p2m],{ ...
    '2017','$R_a=0.69$ cm d$^{-1}$','$R_m=1.80$ cm d$^{-1}$', ...
    '2014','$R_a=0.17$ cm d$^{-1}$','$R_m=0.56$ cm d$^{-1}$'}, ...
    'Interpreter','Latex','Location','NorthWest','Fontsize',12);

pos = get(gca,'Position');
set(gca,'Position',[pos(1) pos(2) pos(3)*3/4 pos(4)*3/4]) 






