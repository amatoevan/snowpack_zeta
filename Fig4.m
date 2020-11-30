% Code to create Figure 4 of NCLIM-20051128B, "A Mechanism for Regional
% Variations in Snowpack Melt Under Rising Temperature" by Evan and
% Eisenman

clearvars
coast = load('coast.mat');


coast.long(coast.long>180) = coast.long(coast.long>180)-360;
d = diff(coast.long);
coast.long(abs(d)>180) = NaN;
    
% Load in global estimates of dzdt made using Eq 2 from the paper and
% reanalysis data

omega = 2*pi/365;
R = 0.34;
Tm = 0.18;
dzdT0 = -1/omega * (1+R) * 1./sqrt(T1.^2-(T0-Tm).^2);
dzdT0(abs(T0-Tm)>T1) = NaN;

dzdT0(land<.1) = NaN;

figure
subplot(2,1,1)
p = pcolor(lon,lat,dzdT0');
hold on
plot(coast.long,coast.lat,'k','linewidth',1)
hold off

grid on
p.LineStyle = 'none';

set(gca,'ticklabelinterpreter','latex','fontsize',12)
set(gca,'ylim',[25 85]+1,'xlim',[-180 180],'ytick',0:15:90, ...
    'xtick',-180:30:180)
xlabel('Longitude ($^\circ$E)','interpreter','latex')
ylabel('Latitude ($^\circ$N)','interpreter','latex')

cb = colorbar;
cb.YTick = -30:3:0;
caxis([-21 -3])
ylabel(cb,'$\partial\zeta/\partial T_0$ (d C$^{-1}$)', ...
    'interpreter','latex','fontsize',12)
cb.TickLabelInterpreter = 'Latex';

pos = get(gca,'Position');
cb.Position = cb.Position.*[1 1 .67 1];
set(gca,'Position',pos.*[1 1 1.02 1])

cm = colormap('Parula');
ind = cumsum(.01:.01:.64);
ind = 64*ind/max(ind);
ind = ceil(ind);
cm = cm(ind,:);
colormap((cm))

ta = text(-185-35,86,'\textbf{a.}','fontsize',12,'Interpreter','Latex');


%---- Argentina
subplot(2,2,4)
p = pcolor(lon,lat,dzdT0');
hold on
plot(coast.long,coast.lat,'k','linewidth',1)
hold off

grid on
p.LineStyle = 'none';

set(gca,'ticklabelinterpreter','latex','fontsize',12)
set(gca,'ylim',[-60 -20],'xlim',[-80 -60],'ytick',-60:10:-20, ...
    'xtick',-80:5:-60)
xlabel('Longitude ($^\circ$E)','interpreter','latex')
ylabel('Latitude ($^\circ$N)','interpreter','latex')

cb4 = colorbar;
caxis([-40 -10])
ylabel(cb4,'$\partial\zeta/\partial T_0$ (d C$^{-1}$)','interpreter','latex','fontsize',12)

pos = get(gca,'Position');
cb4.Position = cb4.Position.*[1 1 .67 1];
set(gca,'Position',pos.*[1 1 1.02 1])

cb4.TickLabelInterpreter = 'Latex';

tc = text(-80-6.5,-20,'\textbf{c.}','fontsize',12,'Interpreter','Latex');



%----- Latitude transect
nh = find(lat>0.1);
sh = find(lat<-0.1);
ynh = nanmean(dzdT0(:,nh));
ysh = nanmean(dzdT0(:,sh));

subplot(2,2,3)
plot(ynh,lat(nh),'linewidth',2);grid
set(gca,'ylim',[20 90],'xlim',[-50 0],'fontsize',12, ...
    'TickLabelInterpreter','Latex','ytick',30:15:90)
ylabel('Latitude ($^\circ$N)','interpreter','latex')
xlabel('$\partial\zeta/\partial T_0$ (d C$^{-1}$)','interpreter','latex')

tb = text(-50-12,90,'\textbf{b.}','fontsize',12,'Interpreter','Latex');


