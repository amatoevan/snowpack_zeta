% Code to create Figure 3 of NCLIM-20051128B, "A Mechanism for Regional
% Variations in Snowpack Melt Under Rising Temperature" by Evan and
% Eisenman
% 
% This code requires the output from Fig2.m and Fig2_vic.m

fs = 10; % fontsize
w = 2*pi/365;


lsb = [173 216 255]/255;
db = [0 0 255]/255;
lim = -40;


%--- Plot the observations 
figure(1)
subplot(2,2,1)

x = fig2.Eq2;
xs = fig2.Eq2_e;
y = fig2.dzdT;
ys = fig2.sigma;

errorbar(x,y,ys,ys,xs,xs,'.','Color',lsb)
grid on; box on; hold on
scatter(x,y,10,db,'fill')
plot([lim 0],[lim 0],'--k')
hold off

set(gca,'xlim',[lim 0],'ylim',[lim 0],'xtick',lim:10:0, ...
    'ytick',lim:10:0,'TickLabelInterpreter','Latex','FontSize',8)
xlabel('$\partial\zeta/\partial T_0$ from Eq. 2 (d C$^{-1}$)', ...
    'Interpreter','Latex','FontSize',fs)
ylabel('$\partial\zeta/\partial T_0$ (d C$^{-1}$)', ...
    'Interpreter','Latex','FontSize',fs)
t1 = text(lim-8,0,'\bf{a.}','FontSize',fs,'Interpreter','Latex');
title('SNOTEL/NARR','Interpreter','Latex')

g = ~isnan(x.*y);
r = corr(x(g),y(g));
bias = mean(y(g)-x(g));
rmse = sqrt(mean( (y(g)-x(g)).^2 ));
disp([r bias rmse])



%--- Plot the VIC output 
subplot(2,2,2)

x = vic.Eq2;
xs = vic.Eq2_e;
y = vic.dzdT;
ys = vic.dzdT_e;

errorbar(x,y,ys,ys,xs,xs,'.','Color',lsb)
hold on
scatter(x,y,10,db,'fill')

grid on; box on; hold on
plot([lim 0],[lim 0],'--k')
hold off

set(gca,'xlim',[lim 0],'ylim',[lim 0],'xtick',lim:10:0, ...
    'ytick',lim:10:0,'TickLabelInterpreter','Latex','FontSize',8)
xlabel('$\partial\zeta/\partial T_0$ from Eq. 2 (d C$^{-1}$)', ...
    'Interpreter','Latex','FontSize',fs)
ylabel('$\partial\zeta/\partial T_0$ (d C$^{-1}$)', ...
    'Interpreter','Latex','FontSize',fs)
t2 = text(lim-8,0,'\bf{b.}','FontSize',fs,'Interpreter','Latex');
ttl = title('VIC','Interpreter','Latex');




