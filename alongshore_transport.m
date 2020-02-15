% Breakpoint formula Ashton, Muarray, Arnoult (2001)
clear
echo off
rhow = 1030.;
g = 9.81;
K = 6.4e4; % m^3/day

dx = 5;
dy = 5;

mrun(1).datt = '20170213_000000';
mrun(1).wlev = 2.21;
mrun(1).rundir = '../swan_stationary/H2017_0d_2p21_stat/';
mrun(1).rname = '2017 3p44m, 0\circ, +2.21 m, FG2';
mrun(1).fnam='../swan_stationary/grids/CCBay_FG2.nc';
mrun(1).pname = '2017_3p44m_0d_2p21_FG2';
mrun(1).descrip = '2013 USACE lidar + 2017-01-25 SfM' ;

mrun(2).datt = '20170213_000000';
mrun(2).wlev = 2.21;
mrun(2).rundir = '../swan_stationary/H2017_0d_2p21_stat5/';
mrun(2).rname = '2017 3p44m, 0\circ, +2.21 m, FG5';
mrun(2).fnam='../swan_stationary/grids/CCBay_FG5.nc';
mrun(2).pname = '2017_3p44m_0d_2p21_FG5';
mrun(2).descrip = '2016-06-06 USGS jetyak + 2016-09-21 SfM'; 

mrun(3).datt = '20170213_000000';
mrun(3).wlev = 2.21;
mrun(3).rundir = '../swan_stationary/H2017_0d_2p21_stat3/';
mrun(3).rname = '2017 3p44m, 0\circ, +2.21 m, FG3';
mrun(3).fnam='../swan_stationary/grids/CCBay_FG3.nc';;
mrun(3).pname = '2017_3p44m_0d_2p21_FG3';
mrun(3).descrip = '2017-04-27 WHOI jetyak + 2017-04-28 SfM' ;

%% loop over all model results

for i=1:length(mrun)
   
   rundir = mrun(i).rundir;
   datt = mrun(i).datt;
   fnam = mrun(i).fnam;
   % Read the grid
   h=ncread(fnam,'h')';
   ss = size(h); if(ss(1)~=126); h=h'; end
   mrun(i).h = h;
   [nrows, ncols] = size(h);
   lon=ncread(fnam,'lon_rho')';
   lat=ncread(fnam,'lat_rho')';
   ss = size(lat); if(ss(1)~=126); lat=lat'; lon=lon'; end
   xg = 0:dx:5*ncols-1;
   yg = 0:dy:5*nrows-1;
   [xgg,ygg]=meshgrid(xg,yg);
   
   hsa=load([rundir,'hsigF.mat']);
   eval(['hs=double(squeeze(hsa.Hsig_',datt,'));']);
   ss = size(hs); if(ss(1)~=126); hs=hs'; end
   mrun(i).hs=hs
   clear hsa
   
   wda=load([rundir,'wdirF.mat']);
   eval(['wd=double(squeeze(wda.Dir_',datt,'));']);
   ss = size(wd); if(ss(1)~=126); wd=wd'; end;
   mrun(i).v=sind(-90-wd+40);
   mrun(i).u=cosd(-90-wd+40);
   clear wda
   
   dia=load([rundir,'dissipF.mat']);
   eval(['diss=double(squeeze(dia.Dissip_',datt,'));']);
   ss = size(diss); if(ss(1)~=126); diss=diss'; end;
   mrun(i).diss=diss;
   clear dia
   
   qbfa=load([rundir,'qbF.mat']);
   eval(['qb=double(squeeze(qbfa.Qb_',datt,'));']);
   ss = size(qb); if(ss(1)~=126); qb=qb'; end;
   mrun(i).qb=qb;
   clear qbfa
   
   wpera=load([rundir,'wperF.mat']);
   eval(['Tp=double(squeeze(wpera.Period_',datt,'));']);
   ss = size(Tp); if(ss(1)~=126); Tp=Tp'; end;
   mrun(i).Tp=Tp;
   clear wpera
   
   
   % calculate wavelength, n, celerity, and wave power
   % Komar, eqn. 5.19-5.21
   
   hv = h(:)+mrun(i).wlev;
   w = (2.*pi)./Tp;
   wv = w(:);
   khv = zeros(size(hv));
   
   for ii=1:length(hv)
      khv(ii) = qkhfs( wv(ii), hv(ii) );
      %fprintf('%f %f %f\n',wv(i),hv(i),khv(i))
   end
   kh = reshape(real(khv),size(h));
   clear hv wv khv
   k = kh./h;
   L = 2.*pi ./k; % n=[m]
   n = 0.5*( 1. + 2.*kh./sinh(2.*kh) ); []
   E = (1./8.)*rhow*g*hs.^2; % [N/m2]
   C = (g*Tp)./(2.*pi) .* tanh(kh); % [m/s]
   Cg = C.*n;
   mrun(i).P = E.*Cg; % [ N m-1 s-2] I think
   
   % find breakpoint
   ibr = zeros(ncols,1);      % breakpoint (50% broken)
   idiss = zeros(ncols,1);    % peak dissipation
   idry = zeros(ncols,1);     % edge of water
   dissmx = zeros(ncols,1);   % value of peak dissipation
   slopedry = zeros(ncols,1); % slope at edge of water
   for ii=1:ncols
      idry(ii) = find(h(:,ii) <= -mrun(i).wlev,1,'last')+1;
      slopedry(ii) = max(0.02, -(h(idry(ii)-1,ii)-h(idry(ii),ii))./dy );
      idx = find( qb(:,ii) >= .5,1,'last');
      if(any(idx))
         ibr(ii) = idx;
      else
         ibr(ii)=0;
      end
      [dissmx(ii), idiss(ii)] = max(diss(:,ii));
   end
   
   % Find profile locations for various quantities
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
   for ii=1:ncols
      if( ibr(ii)>0 )
         Tpb(ii) = Tp(ibr(ii),ii);
         slopeb(ii) = max(0.02, -(h(ibr(ii)-1,ii)-h(ibr(ii),ii))./dy );
         hsb(ii) = hs(ibr(ii),ii);
         khb(ii) = kh(ibr(ii),ii);
         hb(ii) = h(ibr(ii),ii)+mrun(i).wlev;
         Pb(ii) = mrun(i).P(ibr(ii),ii);
         dirb(ii) = wd(ibr(ii),ii);
      end
      Tpdiss(ii) = Tp(idiss(ii),ii);
      hsdiss(ii) = hs(idiss(ii),ii);
      slopediss(ii) = max(0.02, -(h(idiss(ii)-1,i)-h(idiss(ii),ii))./dy );
      khdiss(ii) = kh(idiss(ii),ii);
      hdiss(ii) = h(idiss(ii),ii)+mrun(i).wlev;
      Pdiss(ii) = mrun(i).P(idiss(ii),ii);
      dirdiss(ii) = wd(idiss(ii),ii);
   end
%    dalongb = cosd(-90-dirb+40);
%    dalongdiss = cosd(-90-dirdiss+40);
%    dacrossb = sind(-90-dirb+40);
%    dacrossdiss = sind(-90-dirdiss+40);
   % Sine terms normalizes P for distance alongshore, second term determines
   % alongshore or onshore component per Komar p. 390-391. First version
   % submitted did not include first sine term.
   dalongb = sind(-90-dirb+40).*cosd(-90-dirb+40);
   dalongdiss = sind(-90-dirdiss+40).*cosd(-90-dirdiss+40);
   dacrossb = sind(-90-dirb+40).*sind(-90-dirb+40);
   dacrossdiss = sind(-90-dirdiss+40).*sind(-90-dirdiss+40);
   
   % Calculate runup
   % reverse shoal to get Ho (Nielsen, 2009, eq. 1.7.5)
   Ksb = 1../ sqrt( tanh(khb) .* (1.+2.*khb./sinh(2.*khb)) );
   Hob = hsb./Ksb;
   Ksdiss = 1../ sqrt( tanh(khdiss) .* (1.+2.*khdiss./sinh(2.*khdiss)) );
   Hodiss = hsdiss./Ksdiss;
   R2b = nan*ones(ncols,1);
   R2diss = nan*ones(ncols,1);
   
   % Stockdon formula
   for ii=1:ncols
      [R2b(ii),S,setup, Sinc, SIG, ir, R16] = calcR2(Hob(ii),Tpb(ii),atan(slopedry(ii)),0);
      [R2diss(ii),S,setup, Sinc, SIG, ir, R16] = calcR2(Hodiss(ii),Tpdiss(ii),atan(slopedry(ii)),0);
      
   end
   % calculate Iribarren number at breakpoint according to Battjes (1974); Komar (1988, p. 210) Holthuijsen (2007) and Wikipedia
   for ii=1:ncols
      Lo(ii) = (g/(2.*pi))*Tpb(i).^2;
      eb(ii) = slopeb(ii)./( sqrt( hsb(ii)/Lo(ii) ) );
      eo(ii) = slopeb(ii)./( sqrt( Hob(ii)/Lo(ii) ) );
      % surf-scaling paramter (Guza and Inman, 1975)
      ep(ii) = 0.5*hsb(ii).*(2.*pi/Tpb(ii)).^2/(g*slopeb(ii).^2);
   end
   
   % assign these to the model structure
   varlist = {'ibr','idiss','idry','dissmx','slopedry','Tpb','Tpdiss','slopeb','slopediss',...
      'hsb','hsdiss','khb','khdiss','hb','hdiss','Pb','Pdiss','dirb','dirdiss','R2b','R2diss','eb','eo','ep',...
      'dalongb','dalongdiss','dacrossb', 'dacrossdiss'}
   for iv = 1:length(varlist)
      eval( char( strcat('mrun(i).',varlist(iv),'=',varlist(iv),';')))
   end
%    % ...then clear them
%    for iv = 1:length(varlist)
%       disp( char( strcat(''clear '','',varlist(iv),'')))
%    end
end
%% other numbers we need
nrun = length(mrun)

rhow = 1030.;
g = 9.81;

[nrows, ncols] = size(mrun(1).h)
xg = 0:dx:5*ncols-1;
yg = 0:dy:5*nrows-1;
[xgg,ygg]=meshgrid(xg,yg);

% List of groin locations
gbuf = 20;
gloc = [138 399 582 767 946 1293];

% make a nan list to blank out groins in grid coords
gnang = ones(size(xg));
hold on
for i=1:length(gloc)
   gnang(xg>=(gloc(i)-gbuf) & xg<=(gloc(i)+gbuf)) = NaN;
end
gnang = (gnang');
% Wave power plot
gray = [.8 .8 .8];
dkgray = [.7 .7 .7];
rust = [.8 .2 .2];
pink = [.9 .7 .7];
blue = [.4 .4 .8];
ylims = [-50 20];

load sandwich_vars
% gng = ones(size(xg)); % don't hide near-groin data
% gn = ones(size(xf));
gng = gnang;          % do hide near groin data
gn = gnan;


%% calc maps statistics
hstack = zeros(3,nrows,ncols);
pstack = zeros(3,nrows,ncols);
dstack = zeros(3,nrows,ncols);

for i=1:nrun
   hstack(i,:,:)=-mrun(i).h;
   pstack(i,:,:)=mrun(i).P;
   dstack(i,:,:)=mrun(i).diss;
end
hmin = min(hstack,[],1);
hmax = max(hstack,[],1);
hrange = squeeze( hmax-hmin );
hdiff = diff(hstack,1,1);

pmin = nanmin(pstack,[],1);
pmax = nanmax(pstack,[],1);
pmean = squeeze(nanmean(pstack,1));
prange = squeeze( pmax-pmin )./pmean;
prange(isnan(prange))=0.;

dmin = nanmin(dstack,[],1);
dmax = nanmax(dstack,[],1);
dmean = squeeze(nanmean(dstack,1));
drange = squeeze( dmax-dmin )./dmean;
%% map plots
figure(1)
subplot(2,2,1)
pcolorjw(-mrun(1).h)
caxis([-8 8])
colorbar
title(strcat(mrun(1).descrip,' Bathymetry [m]'))

subplot(2,2,2)
pcolorjw(-mrun(2).h)
caxis([-8 8])
colorbar
title(strcat(mrun(2).descrip,' Bathymetry [m]'))

subplot(2,2,3)
pcolorjw(-mrun(3).h)
caxis([-8 8])
colorbar
title(strcat(mrun(3).descrip,' Bathymetry [m]'))

subplot(2,2,4)
pcolorjw(hrange)
caxis([0,2])
colorbar
title(strcat('Range in Bathymetry [m]'))

%% difference plots
figure(10)
subplot(311)
pcolorjw(squeeze(hdiff(1,:,:)))
caxis([-2 2])
colorbar

subplot(312)
pcolorjw(squeeze(hdiff(2,:,:)))
caxis([-2 2])
colorbar

subplot(313)
pcolorjw(squeeze(hstack(3,:,:)-hstack(1,:,:)))
caxis([-2 2])
colorbar
%%
figure(11); clf
subplot(211)
h=pcolorjw(xg,yg,hrange)
hold on
shading flat
set(h,'facealpha',.5)
xlim([0,1400])
ylim([90,500])
caxis([0, 2])
ylabel('Cross-shore distance [m]','fontsize',14)
xlabel('Alongshore distance [m]','fontsize',14)
%title('Topography/bathymetry [m NAVD88]','fontweight','normal','fontsize',14);
h=colorbar('eastoutside')
ylabel(h,'Range in bathymetry [m]','fontsize',12)
%colormap(cmocean('topo'))
%colormap(cmocean('-deep'))
%colormap(cmocean('-tarn'))
%colormap(cmocean('curl'))
colormap(cmocean('balance'))

% colormap('jet')
almost_white = [.9,.9,.9];
yellow = [.9, .9, 0]
almost_black = [.1, .1, .4]
[c,hndle] = contour(xgg,ygg,-mrun(1).h,[-8:2:6],'linecolor',almost_black,'linewidth',2);
clabel(c,hndle,'fontsize',12,'color',almost_black)
%h1=quiver(200,350,-100*sind(40),100*cosd(40),'color',almost_black,'linewidth',4);
h1=plot_arrow(200,350,200-100*sind(40),350+100*cosd(40),'color','k','linewidth',3,'headwidth',.8);
text(175,320,'North','fontsize',14,'color','k')

subplot(223)
pcolorjw(xg,yg,100*prange)
hold on
caxis([0 200])
xlim([0,1400])
ylim([90,500])
h=colorbar;
ylabel(h,'Range in Wave Power (%)','fontsize',12)
ylabel('Cross-shore distance [m]','fontsize',14)
xlabel('Alongshore distance [m]','fontsize',14)

subplot(224)
pcolorjw(xg,yg,100*drange)
xlim([0,1400])
ylim([90,500])
caxis([0 200])
h=colorbar;
ylabel(h,'Range in Dissipation Rate (%)','fontsize',12)
xlabel('Alongshore distance [m]','fontsize',14)
% print('range_maps.png','-dpng','-r300')
%% plot P grids
figure(2)
subplot(2,2,1)
pcolorjw(mrun(1).P)
title(strcat(mrun(1).descrip,' Power'))
%caxis([-8 8])
colorbar

subplot(2,2,2)
pcolorjw(mrun(2).P)
%caxis([-8 8])
colorbar
title(strcat(mrun(2).descrip,' Power'))

subplot(2,2,3)
pcolorjw(mrun(3).P)
%caxis([-8 8])
colorbar
title(strcat(mrun(3).descrip,' Power'))

subplot(2,2,4)
pcolorjw(prange*100)
hold on
contour(prange*100,[0,100.,200],'-w')
caxis([0,200])
colorbar
title('Range in Power (%)')
%% plot dissipation grids
figure(3)
subplot(2,2,1)
pcolorjw(mrun(1).diss)
title(strcat(mrun(1).descrip,' Dissipation'))
caxis([0 500])
colorbar

subplot(2,2,2)
pcolorjw(mrun(2).diss)
caxis([0 500])
colorbar
title(strcat(mrun(2).descrip,' Dissipation'))

subplot(2,2,3)
pcolorjw(mrun(3).diss)
caxis([0 500])
colorbar
title(strcat(mrun(3).descrip,' Dissipation'))

subplot(2,2,4)
pcolorjw(drange*100)
hold on
contour(drange*100,[0,100.,200],'-w')
caxis([0,200])
colorbar
title('Range in Dissipation (%)')
%%
theta = 0;
phib = [0:5:90.]'
Hb = 2.
Qs = K*Hb.^(5/2)*sind(phib-theta).*cosd(phib-theta)

plot(phib,Qs)
%% plot alongshore components
figure(4); clf;
subplot(513)
hf=fill([65; 65; 125; 125],[ylims'; flipud(ylims')],gray);
set(hf,'edgecolor','none','facealpha',.6)
hold on
hf=fill([265; 265; 340; 340],[ylims'; flipud(ylims')],gray);
set(hf,'edgecolor','none','facealpha',.6)
hf=fill([430; 430; 555; 555],[ylims'; flipud(ylims')],pink);
set(hf,'edgecolor','none','facealpha',.6)
for i=1:nrun
plot(xg,smoothdata(mrun(i).Pdiss.*mrun(i).dalongdiss./1000,'movmean',3).*gng,'-','linewidth',2,'color',[.7 .7 .7]);
plot(xg,smoothdata(mrun(i).Pb.*mrun(i).dalongb./1000,'movmean',3).*gng,'-','linewidth',2,'color',[.8 .2 .2]);
end

% calculate the mean and std dev
Pdiss_sum = zeros(size(mrun(1).Pdiss));
Pdiss_ss = Pdiss_sum;
Pb_sum = Pdiss_sum;
Pb_ss = Pdiss_sum;
for i=1:nrun
   Pdiss_sum = Pdiss_sum + mrun(i).Pdiss.*mrun(i).dalongdiss./1000;
   Pb_sum = Pb_sum + mrun(i).Pb.*mrun(i).dalongdiss./1000;
end
Pdiss_mean = Pdiss_sum ./nrun;
Pb_mean = Pb_sum ./nrun;

plot(xg,Pb_mean.*gng,'-','linewidth',2,'color',[.2 .2 .2]);
plot(xg,Pdiss_mean.*gng,'-','linewidth',2,'color',[1 .2 .2]);

xlim([0,1400])
text(.02,.9,'b','Fontsize',14,'Units','normalized')
ylabel('{\itP_x}  [kW m^{-1}]','fontsize',14)
ylim([-5 15])
set(gca,'xticklabels',[])
grid on
set(gca, 'fontsize', 12)

posb = get(gca, 'Position');
posb(4) = .13;
posb(2) = .27+.16
set(gca, 'Position', posb)

subplot(514)
hf=fill([265; 265; 340; 340],[ylims'; flipud(ylims')],gray);
set(hf,'edgecolor','none','facealpha',.6)
hold on
hf=fill([65; 65; 125; 125],[ylims'; flipud(ylims')],gray);
set(hf,'edgecolor','none','facealpha',.6)
hf=fill([430; 430; 555; 555],[ylims'; flipud(ylims')],pink);
set(hf,'edgecolor','none','facealpha',.6)
for i=1:nrun
   h2=plot(xg(1:end-1)+dx/2,-diff(smoothdata(mrun(i).Pdiss.*mrun(i).dalongdiss./1000,'movmean',5)./dx).*gng(1:end-1),'.k');
   set(h2,'color',[.7 .7 .7],'markersize',14)
   hold on
   h1=plot(xg(1:end-1)+dx/2,-diff(smoothdata(mrun(i).Pb.*mrun(i).dalongb./1000,'movmean',11)./dx).*gng(1:end-1),'.r');
   set(h1,'color',[.8 .2 .2],'markersize',14)
end
% calculate the mean and std dev
dPdiss_sum = zeros(size(diff(mrun(1).Pdiss)));
dPb_sum = dPdiss_sum;
for i=1:nrun
   dPdiss_sum = dPdiss_sum + -diff(gnang.*smoothdata(((mrun(i).Pdiss.*mrun(i).dalongdiss./1000)./dx ),'movmean',13));
   dPb_sum = dPb_sum +       -diff(gnang.*smoothdata(((mrun(i).Pb.*mrun(i).dalongdiss./1000)./dx),'movmean',13));
end
dPdiss_mean = dPdiss_sum ./nrun;
dPb_mean = dPb_sum ./nrun

% plot(xg(1:end-1)+dx/2,smoothdata(dPb_mean,'movmean',13),'-','linewidth',2,'color',[.2 .2 .2]);
% plot(xg(1:end-1)+dx/2,smoothdata(dPdiss_mean,'movmean',13),'-','linewidth',2,'color',[1 .2 .2]);

% plot(xg(1:end-1)+dx/2,smoothdata(-diff(Pb_mean),'movmean',5),'-','linewidth',2,'color',[.2 .2 .2]);
% plot(xg(1:end-1)+dx/2,smoothdata(-diff(Pdiss_mean),'movmean',5),'-','linewidth',2,'color',[1 .2 .2]);

ylim([-0.150 0.1501])
xlim([0,1400])
grid on
set(gca,'xticklabels',[])
text(.02,.9,'c','Fontsize',14,'Units','normalized')
ylabel('-{\Delta} {\itP_x}/{\Delta}{\itx}  [kW m^{-2}]','fontsize',14)
posc = get(gca, 'Position');
posc(4) = .13;
posc(2) = .27

set(gca, 'Position', posc)
set(gca, 'fontsize', 12)

subplot('position',[0.1300    0.1100    0.7750    0.1300])
hf=fill([65; 65; 125; 125],[ylims'; flipud(ylims')],gray);
set(hf,'edgecolor','none','facealpha',.6)
hold on
hf=fill([265; 265; 340; 340],[ylims'; flipud(ylims')],gray);
set(hf,'edgecolor','none','facealpha',.6)
hf=fill([430; 430; 555; 555],[ylims'; flipud(ylims')],pink);
set(hf,'edgecolor','none','facealpha',.6)
h1=plot(xf,medfilt1(dall_vols(:,7),7).*gn,'linewidth',3,'color',rust);
h2=plot(xf,medfilt1(dall_vols(:,7),7).*gn+lb_err(:,7),'--','color',rust);
h3=plot(xf,medfilt1(dall_vols(:,7),7).*gn-lb_err(:,7),'--','color',rust);
xlim([0,1400])
ylim([-50 20])
set(gca, 'fontsize', 12)
grid on
text(.02,.9,'d','Fontsize',14,'Units','normalized')
hy=ylabel('Volume Change [m^3/m]','fontsize',14);
hyp=get(hy,'position')
hyp=hyp
xlabel('Alongshore distance [m]','fontsize',14)
posd = get(gca, 'Position');
posd(4) = .13;
posd(2) = 0.11;
set(gca, 'Position', posd)
print('alongshore_gradients.png','-dpng','-r300')
%%
%% plot hs and runup

% first, load data from 10 Feb survey of dune toe
plot_survey_data

figure(15);clf
ph = .41
gap = (1-(2*ph))./4
px = [2*gap, 3*gap+ph]

ax1=subplot(211);
hf=fill([65; 65; 125; 125],[ylims'; flipud(ylims')],gray);
set(hf,'edgecolor','none','facealpha',.6)
hold on
hf=fill([265; 265; 340; 340],[ylims'; flipud(ylims')],gray);
set(hf,'edgecolor','none','facealpha',.6)
hf=fill([430; 430; 555; 555],[ylims'; flipud(ylims')],pink);
set(hf,'edgecolor','none','facealpha',.6)

for i=1:nrun
h2=plot(xg,smoothdata(mrun(i).hsdiss.*gnang,'movmean',3),'-','linewidth',2,'color',[.7 .7 .7]); hold on
h1=plot(xg,smoothdata(mrun(i).hsb.*gnang,'movmean',3),'-','linewidth',2,'color',[.8 .2 .2]);
end
set(gca, 'fontsize', 12)
ylabel('{\itH_s} [m]','fontsize',14)
plot([430 555],[2.45 2.45],'-r','linewidth',6)
plot([265 340],[2.45 2.45],'-k','linewidth',6,'color',[.6 .6 .6])
text(.02,.95,'a','Fontsize',14,'Units','normalized')

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
hf=fill([65; 65; 125; 125],[ylims'; flipud(ylims')],gray);
set(hf,'edgecolor','none','facealpha',.6)
hold on
hf=fill([265; 265; 340; 340],[ylims'; flipud(ylims')],gray);
set(hf,'edgecolor','none','facealpha',.6)
hf=fill([430; 430; 555; 555],[ylims'; flipud(ylims')],pink);
set(hf,'edgecolor','none','facealpha',.6)
h5=plot(xg,1.43*ones(size(xg)),':k','linewidth',2); %MHHW
h6=plot(xsurvey,medfilt(elev,5),'--k','linewidth',2,'color',[.3 .3 .3]);
for i=1:nrun
h2=plot(xg,smoothdata(mrun(i).R2b.*gnang,'movmean',3),'-','linewidth',2,'color',[.7 .7 .7]); hold on
h1=plot(xg,smoothdata(mrun(i).R2b.*gnang+mrun(i).wlev,'movmean',3),'-','linewidth',2,'color',[.8 .2 .2]);
end
% % h6=plot(xxf,dtoey,'--k'); % from my attempts to find inflection point
% h2=plot(xg,smoothdata(R2b,'movmean',3),'-','linewidth',2,'color',[.7 .7 .7]);
% 
% h1=plot(xg,smoothdata(R2b+wlev,'movmean',3),'-','linewidth',2,'color',[.8 .2 .2]);
set(gca, 'fontsize', 12)
ylabel('Runup {\itR}_2, Elevation [m]','fontsize',14)
h3=plot([430 555],[5.9 5.9],'-r','linewidth',6)
h4=plot([265 340],[5.9 5.9],'-k','linewidth',6,'color',[.6 .6 .6])
%h6=plot([1070 1100],[2.7 2.7],'-r','linewidth',4,'color',[.4 .4 .4]);
legend([h1,h2,h4,h3,h6,h5],'TWL = {\itR}_2 + tides + surge','Runup {\itR}_2','Geotubes','Max. Erosion','Dune Toe','MHHW','location','northeast');
text(.02,.95,'b','Fontsize',14,'Units','normalized')

% text( 245, 4.8,'geotubes','fontsize',12)
ylim([0, 6])
xlim([0,1400])
xlabel('Alongshore distance [m]','fontsize',14)
grid on
pos = get(gca, 'Position')
pos(2)=px(1) 
pos(4)=ph
set(gca, 'Position', pos)
print('hs_runup_revised.png','-dpng','-r300')
