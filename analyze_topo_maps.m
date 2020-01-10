% sandwich_map_analysis
% run this before capecodbay_anal.m

clear
close all

%fnam = 'one_meter_test.nc'
%fnam = 'D:\crs\src\Raster_calcs\20cm_test2.nc'
fnam = 'D:\crs\src\Raster_calcs\one_meter_test2.nc'
%fnam = 'one_meter_test2.nc'

x=ncread(fnam,'Alongshore')';
y=ncread(fnam,'Cross-shore')';
z=ncread(fnam,'__xarray_dataarray_variable__');
z(z<-20)=nan;
[ny,nx,nmap]=size(z);
zd = diff(z,1,3);
zdd = squeeze(zd(:,:,7));
xf = fliplr(x)+110;
zs = nan*z;
% smoothed version of maps
for is=1:nmap
   zs(:,:,is) = conv2(squeeze(z(:,:,is)),ones(3)./9,'same');
end
% work with this range of transects
lgd=1274
fgd=82
%% dates of the 14 maps
titles = ['22-Jan-2016';...
   '25-Jan-2016';...
   '11-Feb-2016';...
   '30-Mar-2016';...
   '21-Sep-2016';...
   '09-Jan-2017';...
   '25-Jan-2017';...
   '14-Feb-2017';...
   '16-Mar-2017';...
   '28-Apr-2017';...
   '04-May-2017';...
   '18-Sep-2017';...
   %         '03-Jan-2017',...
   '10-Jan-2018';...
   '09-Mar-2018'];
titles = cellstr(titles);
dn_map = datenum(titles);
%% find MHW for average map (not including 2018 maps)
zm = mean(z(:,:,1:12),3);
mhwm = nan*ones(nx,1);
for ir=1:nx
   zt = zm(:,ir);
   zok = flipud(zt(~isnan(zt)));
   yok = flipud(y(~isnan(zt))');
   iy = find(zok>1.28,1,'first');
   if(~isempty(iy))
      if(iy>1)
         y0 = interp1(zok(iy-1:iy),yok(iy-1:iy),1.28);
         %y0 = interp1(zok,yok,1.)
      else
         y0 = yok(1);
      end
      mhwm(ir)=y0;
   end
end
% get rid of some outliers
mhwmf = medfilt1(mhwm,15);
ireplace=find(abs(mhwmf-mhwm)>2);
mhwm(ireplace)=mhwmf(ireplace);
fprintf(1,'Replaced %d points.\n',length(ireplace))
% linear interpolation across groin by beachcam house
mhwm(1005:1017) = linspace(126.5,127.3,13);
imhwm = round(mhwm);
%% find Dune toe (2.5 m), MHHW (1.43 m NAVD88) MHW (1.28 m NAVD88) and  MWL (0. m NAVD88) on each map
dt = nan*ones(nx,nmap);
mhhw = nan*ones(nx,nmap);
mhw = nan*ones(nx,nmap);
mwl = nan*ones(nx,nmap);
slope2=nan*ones(nx,nmap);
serr=nan*ones(nx,nmap);
for is=1:nmap 
   for ir=fgd:lgd
      zt = z(:,ir,is);
      %zt = zs(:,ir); % work on smoothed version
      zok = flipud(zt(~isnan(zt)));
      yok = flipud(y(~isnan(zt))');
      iy = find(zok>1.43,1,'first');
      % MHHW
      if(~isempty(iy))
         if(iy>1)
            y0 = interp1(zok(iy-1:iy),yok(iy-1:iy),1.43);
         else
            y0 = yok(1);
         end
         mhhw(ir,is)=y0;
      end
      % MHW
      iy = find(zok>1.28,1,'first');
      if(~isempty(iy))
         if(iy>1)
            y0 = interp1(zok(iy-1:iy),yok(iy-1:iy),1.28);
         else
            y0 = yok(1);
         end
         mhw(ir,is)=y0;
      end
      % MWL
      iy = find(zok>0,1,'first');
      if(~isempty(iy))
         if(iy>1)
            y0 = interp1(zok(iy-1:iy),yok(iy-1:iy),0);
         else
            y0 = yok(1);
         end
         mwl(ir,is)=y0;
      end
      slope2(ir,is)=1.28/abs(mhw(ir,is)-mwl(ir,is));
      serr(ir,is)=0.08/slope2(ir,is);
   end
end
%% List of groin locations
gbuf = 20;
gloc = [138 399 582 767 946 1293];

% make a nan list to blank out groins
gnan = ones(size(xf));
hold on
for i=1:length(gloc)
   gnan(xf>=(gloc(i)-gbuf) & xf<=(gloc(i)+gbuf)) = NaN;
end
gnan = (gnan');
% plot locations of MWL, MHW, and MHHW for maps 7 and 8
figure(1); clf
plot(xf,mwl(:,7).*gnan)
hold on
plot(xf,mwl(:,8).*gnan);
plot(xf,mhw(:,7).*gnan);
plot(xf,mhw(:,8).*gnan);
plot(xf,mhhw(:,7).*gnan);
plot(xf,mhhw(:,8).*gnan);
title('Locations of MWL, MHW, and MHHW')

%% fit linear trend to shoreline MWL and MHW position
mwl_rate = nan*ones(nx,1);
mwl_rate_se = nan*ones(nx,1);
mwl_rate_r2 = nan*ones(nx,1);
mhw_rate = nan*ones(nx,1);
mhw_rate_se = nan*ones(nx,1);
mhw_rate_r2 = nan*ones(nx,1);
mhhw_rate = nan*ones(nx,1);
mhhw_rate_se = nan*ones(nx,1);
mhhw_rate_r2 = nan*ones(nx,1);

% dates of maps
X = dn_map/365.25;
for ir = fgd:lgd
   Y = mwl(ir,:)';
   [a,b,r2,sa,sb,hdot]=lsfit(X,Y,0);
   mwl_rate(ir)=b;
   mwl_rate_se(ir)=sb;
   mwl_rate_r2(ir)=r2;
   Y = mhw(ir,:)';
   [a,b,r2,sa,sb,hdot]=lsfit(X,Y,0);
   mhw_rate(ir)=b;
   mhw_rate_se(ir)=sb;
   mhw_rate_r2(ir)=r2;
   Y = mhhw(ir,:)';
   [a,b,r2,sa,sb,hdot]=lsfit(X,Y,0);
   mhhw_rate(ir)=b;
   mhhw_rate_se(ir)=sb;
   mhhw_rate_r2(ir)=r2;
end
%% plot shoreline retreat rates
figure(2); clf
plot(xf,mhhw_rate);
hold on
plot(xf,mhhw_rate+mhhw_rate_se);
plot(xf,mhhw_rate-mhhw_rate_se);
%% set the distance from MWL to the back of the profile
% these are meant to capture the changing bits
xoff = 80*ones(nx,1);

xoff(1:136)=45;
xoff(81:136)=linspace(45,70,length(81:136));
xoff(137:550)=70;
xoff(550:635)=linspace(70,50,length(550:635));
xoff(636:1005)=50;
xoff(1005:1010)=linspace(50,35,length(1005:1010));
xoff(1011:1046)=35;
xoff(1047:1052)=linspace(35,22,length(1047:1052));
xoff(1052:end)=22;

figure(3); clf
pcolorjw(xf,y,zm)
shading interp
hold on
plot(xf,mhwm,'.b')
plot(xf,mhwm-xoff,'.k')
title('Mean elevation and back of transects')
%% find volume between mhwm and back of profile
dzerr = 0.13 % 1.96 * sqrt ( 0.063^2 + 0.013^2 ) maybe a little optimistic
dune_vols = nan*ones(nx,nmap);
dv_err  = nan*ones(nx,nmap);
hb_vols = nan*ones(nx,nmap);
hb_err  = nan*ones(nx,nmap);
lb_vols = nan*ones(nx,nmap);
lb_err  = nan*ones(nx,nmap);
for is=1:nmap
   for ir=fgd:lgd
      % start on land side, end at MHHW, MHW, or MWL
      istrt= int32(imhwm(ir)-xoff(ir));
      % dune volume (above MHHW)
      iend = int32(round(mhhw(ir,is)));
      if(~isnan(istrt+iend))
         dune_vols(ir,is)=sum(z(istrt:iend,ir,is));
         dv_err(ir,is) = dzerr*(iend-istrt);
      end
      % high beach volume (above MHW)
      iend = int32(round(mhw(ir,is)));
      if(~isnan(istrt+iend))
         hb_vols(ir,is)=sum(z(istrt:iend,ir,is));
         hb_err(ir,is) = dzerr*(iend-istrt);
      end
      % low beach volume (above MHW)
      iend = int32(round(mwl(ir,is)));
      if(~isnan(istrt+iend))
         lb_vols(ir,is)=sum(z(istrt:iend,ir,is));
         lb_err(ir,is) = dzerr*(iend-istrt);
      end
   end
end
all_vols = lb_vols;
hb_vols = hb_vols-dune_vols;
lb_vols = all_vols-hb_vols-dune_vols;
dall_vols = diff(all_vols,1,2);
ddune_vols = diff(dune_vols,1,2);
dhb_vols = diff(hb_vols,1,2);
dlb_vols = diff(lb_vols,1,2);
%% plot of volume change and inferred transport rate
fig=figure(4);clf
fig.PaperUnits = 'inches';
fig.PaperPosition = [ 0 0 9 6 ];

subplot(211)
h1=plot(xf,medfilt(dall_vols(:,7),7).*gnan,'linewidth',3);
hold on
set(h1,'color',[.8 .2 .2])
h2=plot(xf,medfilt1(dall_vols(:,7),7).*gnan+lb_err(:,7),'--');
set(h2,'color',[.8 .2 .2])
h3=plot(xf,medfilt1(dall_vols(:,7),7).*gnan-lb_err(:,7),'--');
set(h3,'color',[.8 .2 .2]);

h1b=plot(xf,medfilt1(ddune_vols(:,7),7).*gnan,'linewidth',3);
hold on
set(h1b,'color',[.2 .2 .8])
h2b=plot(xf,medfilt1(ddune_vols(:,7),7).*gnan+dv_err(:,7),'--');
set(h2b,'color',[.2 .2 .8])
h3b=plot(xf,medfilt1(ddune_vols(:,7),7).*gnan-dv_err(:,7),'--');
set(h3b,'color',[.2 .2 .8]);
xlim([0,1400])
grid on
text(.02,.95,'a','Fontsize',14,'Units','normalized')
ylabel('Volume Change [m^3/m]','fontsize',14)
set(gca,'fontsize',12)

legend([h1b;h1],'Dunes','All','location','southeast')

subplot(212)
dv = flipud((dall_vols(fgd:lgd,7)).*gnan(fgd:lgd));
dv(isnan(dv))=0;
h1c=plot(xf(fgd:lgd),flipud(-cumsum(dv))/(3600*3),'linewidth',2,'color',[.8 .2 .2]);
hold on
dv = flipud((ddune_vols(:,7)));
dv(isnan(dv))=0;
h2c=plot(xf,flipud(-cumsum(dv))/(3600*3),'linewidth',2,'color',[.2 .2 .8]);

xlim([0,1400])
grid on
text(.02,.95,'b','Fontsize',14,'Units','normalized')
ylabel('Inferred Alongshore Flux [m^3/s]','fontsize',14)
xlabel('Alongshore Distance [m]')
set(gca,'fontsize',12)
print('vol_change_cum_flux.png', '-dpng', '-r 200')
%% shoreline location differences
dmwl = diff(mwl,1,2);
dmhw = diff(mhw,1,2);
% edit some wild points
dmwl(115:124,7)=NaN;
dmwl(455:466,7)=NaN;
dmwl(469,7)=NaN;
dmwl(1261:1273,7)=NaN;

dmhw(118,7)=NaN;
dmhw(120,7)=NaN;
dmhw(122:123,7)=NaN;
dmhw(1276:1278,7)=NaN;
dmhw(1284,7)=NaN;
%% plot shoreline change and volume change
fig=figure(3);clf
fig.PaperUnits = 'inches';
fig.PaperPosition = [ 0 0 9 6 ];
subplot(211)
% the median horizontal error of shoreline locations,based on the slope and vertical precision of 8 cm
% is +/-1.9 m...try to make this band about 4-m wide

dmwlf = medfilt1(dmwl(:,7).*gnan);

h1=plot(xf,dmwlf,'linewidth',12,'color',[.9 .8 .8]);
hold on
h1=plot(xf,dmwlf,'linewidth',2,'color',[.9 .2 .2]);
hold on
plot(xf,dmwlf+sqrt(2*serr(:,7)),'--r')
plot(xf,dmwlf-sqrt(2*serr(:,7)),'--r')

dmhwf = medfilt(dmhw(:,7).*gnan);
h2=plot(xf,dmhwf,'linewidth',14,'color',[.8 .8 .9]);
h2=plot(xf,dmhwf,'linewidth',2,'color',[.2 .2 .9]);
plot(xf,dmhwf+sqrt(2*serr(:,7)),'--b')
plot(xf,dmhwf-sqrt(2*serr(:,7)),'--b')
legend([h1;h2],'MWL','MHW','location','southwest')
grid on
ylabel('Shoreline Change [m]','fontsize',14)
text(.02,.95,'a','Fontsize',14,'Units','normalized')
ylim([-40, 10])
xlim([0 1400])
set(gca, 'fontsize', 12)
set(gca,'xticklabels',[])

subplot(212)
h1=plot(xf,medfilt(dall_vols(:,7),7).*gnan,'linewidth',3);
hold on
h2=plot(xf,medfilt(dall_vols(:,7),7).*gnan+lb_err(:,7),'--');
set(h2,'color',[.4 .4 .7])
h3=plot(xf,medfilt(dall_vols(:,7),7).*gnan-lb_err(:,7),'--');
set(h3,'color',[.4 .4 .7]);
xlim([0,1400])
grid on
text(.02,.95,'b','Fontsize',14,'Units','normalized')
ylabel('Volume Change [m^3/m]','fontsize',14)
set(gca, 'fontsize', 12)
xlabel('Alongshore distance [m]','fontsize',14)
print('vol_change_shoreline_change.png', '-dpng', '-r 200')

%% find profile points
% dimension arrays
peaks = nan*ones(nx,2,nmap);
peaklocs = peaks;
inflecs = peaks;
infleclocs = peaks;
vols = nan*ones(nx,nmap);
mhw = nan*ones(nx,nmap);

is = 7;
ir = 970;
for is=1:nmap
   figure(8);clf
   for ir=fgd:lgd
      zt = zs(:,ir);
      zok = flipud(zt(~isnan(zt)));
      yok = flipud(y(~isnan(zt))');
      iy = find(zok>1.28,1,'first');
      if(isempty(iy)), break, end
      if(iy>1)
         y0 = interp1(zok(iy-1:iy),yok(iy-1:iy),1.28);
      else
         y0 = yok(1);
      end
      mhw(ir,is)=y0;
      yr = y0-yok;
      
      iy = find(yr>=0.,1,'first');
      iend = min(length(yr),iy+80);
      yr = yr(iy:iend);
      zok = zok(iy:iend);
      vol(ir,is)=sum(zok); % times dz, which is 1
      
      %plot(yr,zok)
      hold on
      dzy = diff(zok)./diff(yr);
      ydzy = yr(1:end-1)+0.5*diff(yr);
      ddzy = diff(dzy)./diff(yr(2:end));
      yddzy = ydzy(1:end-1)+0.5*diff(ydzy);
      %plot(ydzy,dzy)
      %plot(yddzy,ddzy)
      
      [pks,locs]=findpeaks(zok,yr,'minpeakheight',1.43,'minpeakwidth',5,'npeaks',2);
      if(~isempty(pks))
         %plot(locs,pks,'*r'); % need to add y0 to compensate for offset
         peaks(ir,:,is)=pks;
         peaklocs(ir,:,is)=y0-locs;
         [infl,ilocs]=findpeaks(ddzy,yddzy,'minpeakheight',.05,'npeaks',2,'minpeakwidth',2);
         
         if(~isempty(infl))
            infz = interp1(yr,zok,ilocs);
            %plot(ilocs,infl,'*b');
            inflecs(ir,:,is)=infz;
            infleclocs(ir,:,is)=y0-infl;
         end
      end
   end
end
%% Plot topography and change
is = 7;
ir = 970;
ir2 = 840;
ir3 = 600;
% ir = 5*970
% ir2 = 5*840
% ir3 = 5*600
figure(4); clf
ph = .26
gap = (1-(3*ph))./5
px = [2*gap, 3*gap+ph, 4*gap+2*ph]

ax1=subplot(311);
pcolorjw(xf,y,z(:,:,is))
hold on
if(0)
plot(xf,squeeze(peaklocs(:,:,is)),'.')
plot(xf,squeeze(mhw(:,is)),'.b')
plot(xf,squeeze(infleclocs(:,:,is)),'.')
end
plot([xf(ir) xf(ir)],[50, 150],'--k','linewidth',2)
plot([xf(ir2) xf(ir2)],[50, 150],'--k','linewidth',2)
plot([xf(ir3) xf(ir3)],[50, 150],'--k','linewidth',2)
htp1=text(xf(ir),170,'P1');
set(htp1,'horizontalalignment','center','fontsize',12);
htp2=text(xf(ir2),170,'P2');
set(htp2,'horizontalalignment','center','fontsize',12);
htp3=text(xf(ir3),170,'P3');
set(htp3,'horizontalalignment','center','fontsize',12);
h1=quiver(175,10,-50*sind(52),100*cosd(52),'color',[.2 .2 .2],'linewidth',2);
text(170,40,'North')
text(200,230,'Elevation [m NAVD88] 25-Jan-2017','fontsize',14)
set(gca,'ytick',[0:100:200])
set(gca,'xticklabels',[])
set(gca, 'fontsize', 12)
colormap(ax1,cmocean('-tarn'))
caxis([-8 14])
xlim([100,1400])
ylim([0,250])
colorbar

posa = get(gca, 'Position')
pos(2)=px(3);
posa(3)=0.78
posa(4)=ph;
set(gca, 'Position', posa);


ax2=subplot(312);
pcolorjw(xf,y,z(:,:,is+1))
hold on
plot([xf(ir) xf(ir)],[50, 150],'--r','linewidth',2)
plot([xf(ir2) xf(ir2)],[50, 150],'--r','linewidth',2)
plot([xf(ir3) xf(ir3)],[50, 150],'--r','linewidth',2)
h1=quiver(175,10,-50*sind(52),100*cosd(52),'color',[.2 .2 .2],'linewidth',2);
text(170,40,'North')
ht=text(200,230,'Elevation [m NAVD88] 14-Feb-2017','fontsize',14,'color','r');
set(gca,'ytick',[0:100:200])
set(gca,'xticklabels',[])
ylabel('Cross-shore distance [m]','fontsize',14)
set(gca, 'fontsize', 12)

colormap(ax2,cmocean('-tarn'))
caxis([-8 14])
xlim([100,1400])
ylim([0,250])
cbb = colorbar;

posb = get(gca,'Position');
posb(2)=px(2)
posb(3)=0.78
posb(4)=ph
set(gca, 'Position', posb)

ax3=subplot(313);
v = [-2.5:.25:2.5];
% for plotting purposes, zero out +/- 13 cm
zdd2 = zdd;
zdd2(abs(zdd)<.13)=0;
[C,hc]=contourf(xf,y,zdd2,v,'linestyle','none');
hold on
colormap(ax3,cmocean('-balance'));
caxis([-2.5,2.5])
ts = 'Elevation change [m] 14-Feb minus 25-Jan-2017'
h1=quiver(175,10,-50*sind(52),100*cosd(52),'color',[.2 .2 .2],'linewidth',2);
text(170,40,'North')
ht=text(200,230,ts,'fontsize',14,'color','k');
xlabel('Alongshore distance [m]','fontsize',14)
boxx = [xf(fgd:lgd)';flipud(xf(fgd:lgd)'); xf(fgd)];
boxy = [mhwm(fgd:lgd); flipud(mhwm(fgd:lgd)-xoff(fgd:lgd)); mhwm(fgd)];
plot(boxx,boxy,'-k')
set(gca, 'fontsize', 12);
%axis equal
xlim([100,1400])
ylim([0,250])
cbc=colorbar;
cbc_pos = get(cbc,'position')
cbc_pos=[.9205 px(1) 0.022, .26];
posc =get(gca, 'Position');
posc(2)=px(1);
posc(3)=posb(3);
posc(4)=ph;
set(cbc, 'position',cbc_pos)
set(gca, 'Position', posc);



print('compare_topo.png','-dpng','-r300')
%%
figure(5); clf
subplot(411)
hp1=plot(y,z(:,ir,is),'-k','linewidth',2)
hold on
hp2=plot(y,z(:,ir,is+1),'-r','linewidth',2)
xlim([50, 160])
ylim([-0.5, 6.5])
set(gca,'xticklabels',[])

set(gca, 'fontsize', 12)
text(54,5.7,'P1','fontsize',14)
legend([hp1,hp2],'25-Jan-2017','14-Feb-2017','location','northeast')
grid on

subplot(412)
plot(y,z(:,ir2,is),'-k','linewidth',2)
hold on
plot(y,z(:,ir2,is+1),'-r','linewidth',2)
xlim([50, 160])
ylim([-0.5, 6.5])
set(gca,'xticklabels',[])
set(gca, 'fontsize', 12)
text(54,5.7,'P2','fontsize',14)

ylabel('Elevation [m NAVD88]','fontsize',14)
grid on

subplot(413)
plot(y,z(:,ir3,is),'-k','linewidth',2)
hold on
plot(y,z(:,ir3,is+1),'-r','linewidth',2)
xlim([50, 160])
ylim([-0.5, 6.5])
set(gca, 'fontsize', 12)
text(54,5.7,'P3','fontsize',14)

xlabel('Cross-shore distance [m]','fontsize',14)
grid on
print('compare_profiles.png','-dpng','-r300')
%% save stuff so this dones not have to be re-run prior to stat. analysis
save sandwich_vars xf gnan dall_vols ddune_vols mwl_rate mhw_rate mhhw_rate lb_err