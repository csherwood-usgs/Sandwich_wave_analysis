% capecodbay_anal - Analyze stationary SWAN runs for Sandwich

clear
fnam='../swan_stationary/grids/CCBay_FG2.nc';
%% 2017 runs - peak of Feb storm
datt = '20170213_000000'

wlev = 2.21
rundir = '../swan_stationary/H2017_0d_2p21_stat/'
rname = '2017 3p44m, 0\circ, +2.21 m'
pname = '2017_3p44m_0d_2p21.png'
% wlev = 1.43
% rundir = './H2017_0d_1p43_stat/'
% rname = '2017 3p44m, 0\circ, +1.43 m'
% pname = '2017_3p44m_0d_1p433.png'

%% 2017 runs - event before of Feb storm
% datt = '20170213_000000'
% 
% wlev = 1.73
% rundir = './H2017_0d_1p73_stat/'
% rname = '2017 2p87m, 0\circ, +1.73 m'
% pname = '2017_2p87m_0d_1p73.png'

%% 2016 runs
% set appropriate water level to match SWAN simulation
% datt = '20000101_000000'
% wlev = 1.43;
% wlev = 0;
% if (wlev-1.43)<0.001
% rundir = './H3g0d_1p43_stat/'
% rname = '4.45m, 0\circ, +1.43 m'
% pname = 'Sta4p45_0d_1p433.png'
% 
% % rundir = './H3g40d_1p43_stat/'
% % rname = '4.45m, 40\circ, +1.43 m'
% % pname = 'Sta4p45_40d_1p433.png'
% else
% % rundir = './H3g40d_0_stat/'
% % rname = '4.45m, 40\circ, 0 m'
% % pname = 'Sta4p45_40d_0.png'
% 
% % rundir = './H3g0d_0_stat/'
% % rname = '4.45m, 0\circ, 0 m'
% % pname = 'Sta4p45_0d_0.png'
% end

rhow = 1030.;
g = 9.81;
dx = 5;
dy = 5;
%% Read the grid
h=ncread(fnam,'h')';
[nrows, ncols] = size(h)
lon=ncread(fnam,'lon_rho')';
lat=ncread(fnam,'lat_rho')';
xg = 0:dx:5*ncols-1;
yg = 0:dy:5*nrows-1;
[xgg,ygg]=meshgrid(xg,yg);

%%

hsa=load([rundir,'hsigF.mat']);
eval(['hs(:,:,:)=hsa.Hsig_',datt,';']);
hs=double(squeeze(hs));
clear hsa

wda=load([rundir,'wdirF.mat']);
eval(['wd(:,:,:)=wda.Dir_',datt,';']);
wd=squeeze((double(wd)));
clear wda
v=sind(-90-wd+40);
u=cosd(-90-wd+40);

dia=load([rundir,'dissipF.mat']);
eval(['diss(:,:,:)=dia.Dissip_',datt,';']);
diss=squeeze(double(diss));
clear dia

qbfa=load([rundir,'qbF.mat']);
eval(['qb(:,:,:)=qbfa.Qb_',datt,';']);
qb=squeeze(double(qb));
clear qbfa

wpera=load([rundir,'wperF.mat']);
eval(['Tp(:,:,:)=wpera.Period_',datt,';']);
Tp=squeeze(double(Tp));
clear wpera

%% calculate wavelength, n, celerity, and wave power
% Komar, eqn. 5.19-5.21

hv = h(:)+wlev;
w = (2.*pi)./Tp;
wv = w(:);
khv = zeros(size(hv));

for i=1:length(hv)
   khv(i) = qkhfs( wv(i), hv(i) );
   %fprintf('%f %f %f\n',wv(i),hv(i),khv(i))
end
kh = reshape(real(khv),126,288);
clear hv wv khv
k = kh./h;
L = 2.*pi ./k; % n=[m]
n = 0.5*( 1. + 2.*kh./sinh(2.*kh) ); []
E = (1./8.)*rhow*g*hs.^2; % [N/m2]
C = (g*Tp)./(2.*pi) .* tanh(kh); % [m/s]
Cg = C.*n;
P = E.*Cg; % [ N m-1 s-2] I think 
%% find breakpoint
ibr = zeros(ncols,1);      % breakpoint (50% broken)
idiss = zeros(ncols,1);    % peak dissipation
idry = zeros(ncols,1);     % edge of water
dissmx = zeros(ncols,1);   % value of peak dissipation
slopedry = zeros(ncols,1); % slope at edge of water
for i=1:ncols
   idry(i) = find(h(:,i) <= -wlev,1,'last')+1;
   slopedry(i) = max(0.02, -(h(idry(i)-1,i)-h(idry(i),i))./dy );
   idx = find( qb(:,i) >= .5,1,'last');
   if(any(idx))
   ibr(i) = idx;
   else
      ibr(i)=0;
   end
   [dissmx(i), idiss(i)] = max(diss(:,i));
end
%% Find profile locations for various quantities

Tpb = nan*ones(ncols,1);
Tpdiss = nan*ones(ncols,1);
slopeb = nan*ones(ncols,1);
slopediss = nan*ones(ncols,1);
hsb = nan*ones(ncols,1);
hsdiss = nan*ones(ncols,1);
khb = nan*ones(ncols,1);
khdiss = nan*ones(ncols,1);
hb = nan*ones(ncols,1);
hdiss = nan*ones(ncols,1);
Pb = nan*ones(ncols,1);
Pdiss = nan*ones(ncols,1);
dirb = nan*ones(ncols,1);
dirdiss = nan*ones(ncols,1);
for i=1:ncols
   if( ibr(i)>0 )
      Tpb(i) = Tp(ibr(i),i);
      slopeb(i) = max(0.02, -(h(ibr(i)-1,i)-h(ibr(i),i))./dy );
      hsb(i) = hs(ibr(i),i);
      khb(i) = kh(ibr(i),i);
      hb(i) = h(ibr(i),i)+wlev;
      Pb(i) = P(ibr(i),i);
      dirb(i) = wd(ibr(i),i);
   end
   Tpdiss(i) = Tp(idiss(i),i);
   hsdiss(i) = hs(idiss(i),i);
   slopediss(i) = max(0.02, -(h(idiss(i)-1,i)-h(idiss(i),i))./dy );
   khdiss(i) = kh(idiss(i),i);
   hdiss(i) = h(idiss(i),i)+wlev;
   Pdiss(i) = P(idiss(i),i);
   dirdiss(i) = wd(idiss(i),i);
end
dalongb = cosd(-90-dirb+40);
dalongdiss = cosd(-90-dirdiss+40);
dacrossb = sind(-90-dirb+40);
dacrossdiss = sind(-90-dirdiss+40);
%% Calculate runup
% reverse shoal to get Ho (Nielsen, 2009, eq. 1.7.5)
Ksb = 1../ sqrt( tanh(khb) .* (1.+2.*khb./sinh(2.*khb)) );
Hob = hsb./Ksb;
Ksdiss = 1../ sqrt( tanh(khdiss) .* (1.+2.*khdiss./sinh(2.*khdiss)) );
Hodiss = hsdiss./Ksdiss;
R2b = nan*ones(ncols,1);
R2diss = nan*ones(ncols,1);

% Stockdon formula
for i=1:ncols
   [R2b(i),S,setup, Sinc, SIG, ir, R16] = calcR2(Hob(i),Tpb(i),atan(slopedry(i)),0);
   [R2diss(i),S,setup, Sinc, SIG, ir, R16] = calcR2(Hodiss(i),Tpdiss(i),atan(slopedry(i)),0);

end

%% calculate Iribarren number at breakpoint according to Battjes (1974); Komar (1988, p. 210) Holthuijsen (2007) and Wikipedia
for i=1:ncols
Lo(i) = (g/(2.*pi))*Tpb(i).^2;
eb(i) = slopeb(i)./( sqrt( hsb(i)/Lo(i) ) );
eo(i) = slopeb(i)./( sqrt( Hob(i)/Lo(i) ) );
% surf-scaling paramter (Guza and Inman, 1975)
ep(i) = 0.5*hsb(i).*(2.*pi/Tpb(i)).^2/(g*slopeb(i).^2);
end
%% Bathy figure
minush = -h;
figure(1); clf
subplot(211)
pcolorjw(xg,yg,minush)
hold on
shading flat
xlim([0,1400])
ylim([50,500])
caxis([-8, 10])
ylabel('Cross-shore distance [m]','fontsize',14)
xlabel('Alongshore distance [m]','fontsize',14)
title('Topography/bathymetry [m NAVD88]','fontweight','normal','fontsize',14);
colorbar('northoutside')
%colormap(cmocean('topo'))
%colormap(cmocean('-deep'))
colormap(cmocean('-tarn'))

% colormap('jet')
[c,hndle] = contour(xgg,ygg,minush,[-8:2:6],'linecolor',[.4 .4 .4]);
clabel(c,hndle)
h1=quiver(150,350,-100*sind(40),100*cosd(40),'color',[.2 .2 .2],'linewidth',2);
text(125,320,'North')
% %axis equal
% for i=1:ncols
%    plot(xg(i),yg(idiss(i)),'.','color',[.7 .7 .7],'markersize',12)
%    %plot(xg(i),yg(idry(i)) ,'.','color',[.7 .7 .7],'markersize',12)
%    if (ibr(i)>0)
%       plot(xg(i),yg(ibr(i)),'.','color',[.8 .2 .2],'markersize',12)
%    end
% end

%% Wave power plot
figure(2);clf
subplot(211)
pcolorjw(xg,yg,P./1000)
hold on
shading flat
xlim([0,1400])
ylim([50,500])
%colormap(cmocean('haline'))
colormap('parula')
caxis([0 40])
ylabel('Cross-shore distance [m]','fontsize',14)
[c,hndle] = contour(xgg,ygg,minush,[-8:2:6],'linecolor',[.8 .8 1],'linewidth',2);
clabel(c,hndle,'Fontsize',12,'color',[.8 .8 1])

%axis equal
for i=1:ncols
   plot(xg(i),yg(idiss(i)),'.','color',[.7 .7 .7],'markersize',12)
   %plot(xg(i),yg(idry(i)) ,'.','color',[.7 .7 .7],'markersize',12)
   if (ibr(i)>0)
      plot(xg(i),yg(ibr(i)),'.','color',[.8 .2 .2],'markersize',12)
   end
end

h1=quiver(xgg(1:8:end,1:8:end),ygg(1:8:end,1:8:end),u(1:8:end,1:8:end),v(1:8:end,1:8:end));
set(h1,'color',[.3 .3 .4],'linewidth',1)
title('Wave Power {\itP} [kW m^{-1}]','fontweight','normal','fontsize',14);
set(gca,'xticklabels',[])
colorbar('northoutside')
h1=quiver(150,350,-100*sind(40),100*cosd(40),'color',[.2 .2 .2],'linewidth',2);
text(125,320,'North')
set(gca, 'fontsize', 12)

subplot(413)
plot(xg,smoothdata(Pdiss.*dalongdiss./1000,'movmean',3),'-','linewidth',2,'color',[.7 .7 .7]);
hold on
plot(xg,smoothdata(Pb.*dalongb./1000,'movmean',3),'-','linewidth',2,'color',[.8 .2 .2]);
xlim([0,1400])
ylabel('{\itP_x}  [kW m^{-1}]','fontsize',14)
xlim([0,1400])
set(gca,'xticklabels',[])
grid on
pos = get(gca, 'Position');
pos(4) = pos(4)*1.3;
set(gca, 'Position', pos)
set(gca, 'fontsize', 12)

subplot(414)
h2=plot(xg(1:end-1)+dx/2,-diff(smoothdata(Pdiss.*dalongdiss./1000,'movmean',5)./dx),'.k');
set(h2,'color',[.7 .7 .7],'markersize',14)
hold on
h1=plot(xg(1:end-1)+dx/2,-diff(smoothdata(Pb.*dalongb./1000,'movmean',11)./dx),'.r');
set(h1,'color',[.8 .2 .2],'markersize',14)
ylim([-0.0500 0.05001])
xlim([0,1400])

grid on
ylabel('-{\Delta} {\itP_x}/{\Delta}{\itx}  [kW m^{-2}]','fontsize',14)
xlabel('Alongshore distance [m]','fontsize',14)
pos = get(gca, 'Position');
pos(4) = pos(4)*1.3;
set(gca, 'Position', pos)
set(gca, 'fontsize', 12)
print('wave_power_plot.png','-dpng','-r300') 
%% Plot runup and slope
figure(3);clf
ax1=subplot(211);
h2=plot(xg,smoothdata(R2b,'movmean',3),'-','linewidth',2,'color',[.7 .7 .7]);
hold on
h1=plot(xg,smoothdata(R2b+wlev,'movmean',3),'-','linewidth',2,'color',[.8 .2 .2]);

set(gca, 'fontsize', 12)
ylabel('Runup {\itR}_2, Elevation [m]','fontsize',14)
legend([h1,h2],'Runup {\itR}_2','Runup + water level','location','northeast')
set(gca,'xticklabels',[])
xlim([0,1400])

grid on

ax2=subplot(212);
h2=plot(xg,smoothdata(slopediss,'movmean',3),'-','linewidth',2,'color',[.7 .7 .7]);
hold on
h1=plot(xg,smoothdata(slopeb,'movmean',3),'-','linewidth',2,'color',[.8 .2 .2]);
set(gca, 'fontsize', 12)
ylabel('Slope tan{\alpha} [ ]','fontsize',14)
legend([h1,h2],'Slope at breakpoint','Slope at max. diss.', 'location','northeast')
xlabel('Alongshore distance [m]','fontsize',14)
grid on
xlim([0,1400])
print('R2_slopeplot.png','-dpng','-r300') 

%% Plot hs and h
figure(4);clf
ax1=subplot(211);
h2=plot(xg,smoothdata(hsdiss,'movmean',3),'-','linewidth',2,'color',[.7 .7 .7]);
hold on
h1=plot(xg,smoothdata(hsb,'movmean',3),'-','linewidth',2,'color',[.8 .2 .2]);
set(gca, 'fontsize', 12)
ylabel('{\itH_s} [m]','fontsize',14)
legend([h1,h2],'{\itH_s} at breakpoint','{\itH_s} at max. diss.','location','northwest')
set(gca,'xticklabels',[])
xlim([0,1400])

grid on

ax2=subplot(212);

h2=plot(xg,smoothdata(hdiss,'movmean',3),'-','linewidth',2,'color',[.7 .7 .7]);
hold on
h1=plot(xg,smoothdata(hb,'movmean',3),'-','linewidth',2,'color',[.8 .2 .2]);
set(gca, 'fontsize', 12)
ylabel('Depth [m]','fontsize',14)
legend([h1,h2],'{\ith} at breakpoint','{\ith} at max. diss.', 'location','northwest')
xlabel('Alongshore distance [m]','fontsize',14)
grid on
xlim([0,1400])
print('hs_hplot.png','-dpng','-r300') 

%% plot hs and runup
figure(5);clf
ph = .41
gap = (1-(2*ph))./4
px = [2*gap, 3*gap+ph]

ax1=subplot(211);
h2=plot(xg,smoothdata(hsdiss,'movmean',3),'-','linewidth',2,'color',[.7 .7 .7]);
hold on
h1=plot(xg,smoothdata(hsb,'movmean',3),'-','linewidth',2,'color',[.8 .2 .2]);
set(gca, 'fontsize', 12)
ylabel('{\itH_s} [m]','fontsize',14)
plot([390 500],[2.4 2.4],'-r','linewidth',6)
plot([230 390],[2.4 2.4],'-r','linewidth',6,'color',[.9 .7 .7])
plot([100 230],[2.4 2.4],'-r','linewidth',6)
%text( 245, 2.3,'geotubes','fontsize',12)

legend([h1,h2],'{\itH_s} at breakpoint','{\itH_s} at max. diss.','location','southwest')
set(gca,'xticklabels',[])
ylim([0 2.5])
xlim([0,1400])
grid on
pos = get(gca, 'Position')
pos(2)=px(2) 
pos(4)=ph
set(gca, 'Position', pos)

ax1=subplot(212);
h2=plot(xg,smoothdata(R2b,'movmean',3),'-','linewidth',2,'color',[.7 .7 .7]);
hold on
h1=plot(xg,smoothdata(R2b+wlev,'movmean',3),'-','linewidth',2,'color',[.8 .2 .2]);
set(gca, 'fontsize', 12)
ylabel('Runup {\itR}_2, Elevation [m]','fontsize',14)
h3=plot([390 500],[4.8 4.8],'-r','linewidth',6);
h4=plot([230 390],[4.8 4.8],'-r','linewidth',6,'color',[.9 .7 .7]);
h5=plot([100 230],[4.8 4.8],'-r','linewidth',6);
h6=plot([1070 1100],[2.7 2.7],'-r','linewidth',4,'color',[.4 .4 .4]);
legend([h1,h2,h3,h4,h6],'TWL = {\itR}_2 + tides + surge','Runup {\itR}_2','Max. Erosion','Geotubes','Overwash channel','location','northeast');

% text( 245, 4.8,'geotubes','fontsize',12)
ylim([0, 5])
xlim([0,1400])
xlabel('Alongshore distance [m]','fontsize',14)
grid on
pos = get(gca, 'Position')
pos(2)=px(1) 
pos(4)=ph
set(gca, 'Position', pos)
print('hs_runup.png','-dpng','-r300')
%% beach slope calcs and plot
dhx = diff(h);
dhy = diff(h')';
dhxm = [NaN*ones(1,288);dhx];
dhym = [NaN*ones(126,1) dhy];
slopemap = sqrt(dhxm.^2 + dhym.^2 );

dhx = diff(h);
dhy = diff(h')';
dhxm = [NaN*ones(1,288);dhx];
dhym = [NaN*ones(126,1) dhy];
slopemap = sqrt(dhxm.^2 + dhym.^2 );
nanmaskh=ones(size(h));
nanmaskh( -h >=3 )=NaN;
figure(8); clf
subplot(211)
pcolor(xg,yg,slopemap.*nanmaskh)
hold on
shading flat
xlim([0,1400])
ylim([50,500])
%colormap(cmocean('haline'))
colormap('parula')
caxis([0 .6])
ylabel('Cross-shore distance [m]','fontsize',14)
[c,hndle] = contour(xgg,ygg,minush,[-8:2:6],'linecolor',[.8 .8 1],'linewidth',2);
clabel(c,hndle,'Fontsize',12,'color',[.8 .8 1])
%axis equal
for i=1:ncols
   plot(xg(i),yg(idiss(i)),'.','color',[.7 .7 .7],'markersize',12)
   %plot(xg(i),yg(idry(i)) ,'.','color',[.7 .7 .7],'markersize',12)
   if (ibr(i)>0)
      plot(xg(i),yg(ibr(i)),'.','color',[.8 .2 .2],'markersize',12)
   end
end
title('Bottom slope tan{\alpha} [ ]','fontweight','normal','fontsize',14);
set(gca,'xticklabels',[])
colorbar('northoutside')

subplot(413)
plot(xg,slopeb,'linewidth',2)
ylabel('Slope at breakpoint tan(\alpha)_b')
xlim([0,1400])
set(gca,'xticklabels',[])
grid on
pos = get(gca, 'Position');
pos(4) = pos(4)*1.3;
set(gca, 'Position', pos)
set(gca, 'fontsize', 12)

subplot(414)
plot(xg,2.5*ones(size(xg)),'--k')
hold on
plot(xg,20*ones(size(xg)),'--k')
plot(xg,ep,'linewidth',2)
ylabel('Surf-scaling parameter \epsilon')
xlim([0,1400])
set(gca,'yscale','log')
pos = get(gca, 'Position');
pos(4) = pos(4)*1.3;
text(10, 30, 'Dissipative','fontsize',12)
text(10, 1, 'Reflective','Fontsize',12)
set(gca, 'Position', pos)
set(gca, 'fontsize', 12)
xlabel('Alongshore distance [m]','fontsize',14)
print('slope_map_surf_scaling.png','-dpng','-r300')
%% slope stats
slopemapNaN = slopemap.*nanmaskh;
hmapNaN = h.*nanmaskh;
ok = ~isnan(slopemapNaN(:));
slopelist = slopemapNaN(ok);
hlist = hmapNaN(ok);
figure(9); clf
plot(-hlist,slopelist,'.')
ylim([0 3])