% Code to create Figure 2 of NCLIM-20051128B, "A Mechanism for Regional
% Variations in Snowpack Melt Under Rising Temperature" by Evan and
% Eisenman
%
% Also requires reanalysis T and P interpolated to SNOTEL stations

coast = load('coast.mat');
w = 2*pi/365;

% Load SNOTEL data. This can be acquired from:
% https://doi.pangaea.de/10.1594/PANGAEA.896396

% Do a little more QC on snotel data
t = (1:365)';
for i = 1:size(sno.a,1)
    for j = 1:size(sno.a,2)
        
        s = sno.swe(i,t>sno.f(i,j)&t<sno.m(i,j),j)';
        if numel(find(isnan(s))) >= 0.9*numel(s) || ...
                numel(find(s==0)) >= 0.9*numel(s)
            sno.swe(i,:,j) = NaN;
            sno.p(i,j) = NaN;
            sno.z(i,j) = NaN;
        end
        
        s = sno.swe(i,t>sno.m(i,j)&t<sno.z(i,j),j)';
        if numel(find(isnan(s))) >= 0.9*numel(s) || ...
                numel(find(s==0)) >= 0.9*numel(s)
            sno.swe(i,:,j) = NaN;
            sno.p(i,j) = NaN;
            sno.z(i,j) = NaN;
        end
        
    end
end
clearvars i j s t


% Calculate Ra/Rm
b = sno.p>500 | sno.p<10;
sno.p(b) = NaN;
sno.m(b) = NaN;
sno.f(b) = NaN;
sno.Ra = sno.p./(sno.m-sno.f);
sno.Rm = sno.p./(sno.z-sno.m);
sno.Ram = nanmean(sno.Ra./sno.Rm,2);
clearvars b


% Calculate LTM T0, T1, and PHI from reanalysis
load('../../data/T.daily.mat');
T0 = NaN(size(T,1),1);
T1 = NaN(size(T,1),1);
PHI = NaN(size(T,1),1);
Trmse = NaN(size(T,1),1);
for i = 1:size(T,1)
    [res,~] = fitsine((1:365)', T(i,:)'-273.15);
    T0(i) = res.c;
    T1(i) = res.a;
    PHI(i) = res.b;
    Trmse(i) = sqrt(mean( (res((1:365)') - (T(i,:)' - 273.15)).^2 ));
end
narr.T0 = T0;
narr.T1 = T1;
narr.PHI = PHI;
narr.RMSE = Trmse;


% Load annual mean T & P
load('../../data/T.annual.mat');
load('../../data/P.annual.mat'); % Same as using Oct-Mar/Apr
T = T - 273.15;

% Correct biases in T via SNOTEL
snotelT = load('SNOTEL.T0T1PHI.mat');

% Via SNOTEL
T0 = mean(T(:,20:37),2); % only 2001-2018
dT = snotelT.T0 - T0; % via snotel

% apply correction
T = T + dT;
T0 = mean(T,2);

% Also correct T1
T1 = narr.T1;
dT1 = snotelT.T1 - T1;
T1 = T1 + dT1;


% Observational dzdTm
dzdT_p = NaN(size(T,1),1);
dzdT_e = NaN(size(T,1),1);
dSdT_p = NaN(size(T,1),1);  % trend in 4/1 snowpack
dTdt = NaN(size(T,1),1);    % actual trend (C/yr)
dmxdT = NaN(size(T,1),1);   % partial(peak snowpack day)/
for i = 1:size(T,1)
    
    if ~isnan(sum(T(i,:)))
    y = sno.z(i,:)';
    
    [b,bi] = regress(y,[ones(size(y)) T(i,:)' P(i,:)']);
    dzdT_e(i) = b(2)-bi(2,1);
    dzdT_p(i) = b(2);
    
    y = squeeze(sno.swe(i,184,:));
    b = regress(y,[ones(size(y)) T(i,:)' P(i,:)']);
    dSdT_p(i) = b(2);
    
    y = sno.m(i,:)';
    b = regress(y,[ones(size(y)) T(i,:)' P(i,:)']);
    dmxdT(i) = b(2);
    
    y = T(i,:)';
    x = (1:numel(T(i,:)))';
    b = regress(y,[ones(size(y)) x P(i,:)']);
    dTdt(i) = b(2);
    
    end
    
end
T_ = T;
P_ = P;
clearvars T P i y b



%--- CALCULATE Tm 
% Use the LT means to get Tm
Tm = NaN(size(sno.Ram));
T = nanmean(T,3);
wd = 360/365;
for i = 1:length(sno.Ram)
    s = nanmean(sno.swe(i,:,:),3)';
    d = diff(s);
        
    % Remove the fitted sine wave (via narr coefficients), then add back in
    % the fitted sine wave (via corrected coefficients).
    t = T(i,:)' - ...
        (-narr.T1(i)*sind(wd*(1:365)'+narr.PHI(i))+narr.T0(i));
    t = t + (-T1(i)*sind(wd*(1:365)'+narr.PHI(i))+T0(i));    
    t = t(1:end-1) + diff(t)/2;
    t = t(d<0);
    Tm(i) = min(t); % coolest temps when snow melts
    
end
clearvars s d i



% Calculate Eq2
Ra = nanmean(sno.Ra,2);
Rm = nanmean(sno.Rm,2);
Ram = sno.Ram;
Eq2 = NaN(size(dzdT_p));


tm = nanmean(Tm);
for i = 1:398
    ram = Ram(i); %Ra(i)./Rm(i);
    eq2 = -1/w.*(1+ram)*1./sqrt(T1(i)^2-(T0(i)-tm).^2);
    eq2(T1(i)^2 < (T0(i)-tm).^2) = NaN;
    Eq2(i) = nanmean(eq2);
end





%------------------------------------------------------------------
% Map the data
%------------------------------------------------------------------
figure
x = sno.slon;
y = sno.slat;

b = isnan(dzdT_p.*Eq2);
dzdT_p(b) = NaN;
Eq2(b) = NaN;

%--- Observations
subplot(2,2,1)

% setup the Map environment
ax1 = usamap({'CA','MT'});
stxt = {'California','Oregon','Washington','Arizona','Nevada', ...
    'Utah','Idaho','Montana','Wyoming','Colorado','New Mexico', ...
    'North Dakota','South Dakota','Nebraska'};
states = shaperead('usastatelo', 'UseGeoCoords', true,'Selector',...
  {@(name) any(strcmp(name,stxt)), 'Name'});
faceColors = makesymbolspec('Polygon',{'INDEX', [1 numel(states)],...
   'FaceColor',repmat([1 1 1],numel(states),1)});
gs = geoshow(ax1,states,'DisplayType','polygon','SymbolSpec',faceColors,...
    'EdgeColor',[0 0 0]);
setm(ax1,'fontsize',8)
ax1.TickLabelInterpreter = 'Latex';

t1 = title('Observations','Interpreter','Latex');
t1t = text(-1.5485e6,6.2737e6,'\bf{a.}','Interpreter','Latex', ...
    'fontsize',8);

% Plot the data
[~,ix] = sort(dzdT_p);
ix = flipud(ix);
s1 = scatterm(y(ix),x(ix),8,dzdT_p(ix),'filled');
caxis([-30 0])



%--- Theory
subplot(2,2,2)

% setup the Map environment
ax2 = usamap({'CA','MT'});
stxt = {'California','Oregon','Washington','Arizona','Nevada', ...
    'Utah','Idaho','Montana','Wyoming','Colorado','New Mexico', ...
    'North Dakota','South Dakota','Nebraska'};
states = shaperead('usastatelo', 'UseGeoCoords', true,'Selector',...
  {@(name) any(strcmp(name,stxt)), 'Name'});
faceColors = makesymbolspec('Polygon',{'INDEX', [1 numel(states)],...
   'FaceColor',repmat([1 1 1],numel(states),1)});
geoshow(ax2,states,'DisplayType','polygon','SymbolSpec',faceColors, ...
    'EdgeColor',[0 0 0]);
setm(ax2,'fontsize',8)
ax2.TickLabelInterpreter = 'Latex';
t2 = title('Idealized Model','Interpreter','Latex');
t2t = text(-1.5485e6,6.2737e6,'\bf{b.}','Interpreter','Latex', ...
    'fontsize',8);


% Plot the data
[~,ix] = sort(Eq2);
ix = flipud(ix);
s2b = scatterm(y(ix),x(ix),8,Eq2(ix),'filled');



% Create the colorbar
cb2 = colorbar;
caxis([-30 0])
ylabel(cb2,'$\partial\zeta/\partial T_0$ (d C$^{-1}$)', ...
    'Interpreter','Latex','FontSize',8);
cb2.TickLabelInterpreter = 'Latex';
cb2.FontSize = 8;
cb2pos = cb2.Position;
cb2.Position = cb2.Position.*[1 1 0.67 1];
cb2.Position = cb2.Position+[0.04 0 0 0];



% Move over the first plot
ax1.Position = ax1.Position + [0.1 0 0 0];


