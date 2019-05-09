% sandwich_map_analysis
clear
%fnam = 'one_meter_test.nc'
%fnam = 'D:\crs\src\Raster_calcs\20cm_test2.nc'
fnam = 'D:\crs\src\Raster_calcs\one_meter_test2.nc'

x=ncread(fnam,'Alongshore')';
y=ncread(fnam,'Cross-shore')';
z=ncread(fnam,'__xarray_dataarray_variable__');
z(z<-20)=nan;
[ny,nx,nmap]=size(z);
zd = diff(z,1,3);
zdd = squeeze(zd(:,:,7));
xf = fliplr(x)+110;

% first good profile:
irg = 82;
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
mhwmf = medfilt(mhwm,15);
ireplace=find(abs(mhwmf-mhwm)>2)
mhwm(ireplace)=mhwmf(ireplace);
fprintf(1,'Replaced %d points.\n',length(ireplace))
% linear interpolation across groin by beachcam house
mhwm(1005:1017) = linspace(126.5,127.3,13)
imhwm = round(mhwm);
%% set the distance to the back of the profile
% these are meant to capture the changing bits
xoff = 80*ones(nx,1);
xoff(1046:end)=20;
%xoff(1006:end)=20;
xoff(1006:1045)=35;
xoff(1:136)=45;
xoff(636:1005)=50;
xoff(137:635)=70;
figure(1); clf
pcolorjw(xf,y,zm)
shading interp
hold on
plot(xf,mhwm,'.b')
plot(xf,mhwm-xoff,'.k')
title('Mean elevation')
%% find volume between mhwm and back of profile
vols = nan*ones(nx,nmap);
for is=1:nmap
for ir=irg:nx
   istrt= imhwm(ir)-xoff(ir)
   iend = imhwm(ir)
   if(~isnan(istrt+iend))
   vols(ir,is)=sum(z(imhwm(ir)-xoff(ir):imhwm(ir),ir,is));
   end
end
end
%%
plot(xf,dvols(:,7),'linewidth',3)
hold on
plot(xf,medfilt(dvols(:,7),7),'linewidth',3)

%% find profile points


% dimension
peaks = nan*ones(nx,2,nmap);
peaklocs = peaks;
inflecs = peaks;
infleclocs = peaks;
vols = nan*ones(nx,nmap);
mhw = nan*ones(nx,nmap);

is = 7
ir = 970;
for is=1:nmap
   % smooth the map with 3x3 boxcar
   zs = conv2(squeeze(z(:,:,is)),ones(3)./9,'same');
   figure(8);clf
   for ir=20:nx
      zt = zs(:,ir);
      zok = flipud(zt(~isnan(zt)));
      yok = flipud(y(~isnan(zt))');
      iy = find(zok>1.28,1,'first');
      if(isempty(iy)), break, end
      if(iy>1)
         y0 = interp1(zok(iy-1:iy),yok(iy-1:iy),1.28);
         %y0 = interp1(zok,yok,1.)
      else
         y0 = yok(1)
      end
      mhw(ir,is)=y0;
      yr = y0-yok;
      
      iy = find(yr>=0.,1,'first');
      iend = min(length(yr),iy+80);
      yr = yr(iy:iend);
      zok = zok(iy:iend);
      vol(ir,is)=sum(zok); % times dz, which is 1
      
      plot(yr,zok)
      hold on
      dzy = diff(zok)./diff(yr);
      ydzy = yr(1:end-1)+0.5*diff(yr);
      ddzy = diff(dzy)./diff(yr(2:end));
      yddzy = ydzy(1:end-1)+0.5*diff(ydzy);
      %plot(ydzy,dzy)
      plot(yddzy,ddzy)
      
      [pks,locs]=findpeaks(zok,yr,'minpeakheight',1.43,'minpeakwidth',5,'npeaks',2)
      if(~isempty(pks))
         plot(locs,pks,'*r'); % need to add y0 to compensate for offset
         peaks(ir,:,is)=pks;
         peaklocs(ir,:,is)=y0-locs;
         [infl,ilocs]=findpeaks(ddzy,yddzy,'minpeakheight',.05,'npeaks',2,'minpeakwidth',2);
      end
      if(~isempty(infl))
         infz = interp1(yr,zok,ilocs);
         plot(ilocs,infl,'*b');
         inflecs(ir,:,is)=infz;
         infleclocs(ir,:,is)=y0-infl;
      end
   end
end
%plot(yr(3:end),ddzy)
%%
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
titles = cellstr(titles)

%%
is = 7;
ir = 970;
ir2 = 840;
ir3 = 600;
% ir = 5*970
% ir2 = 5*840
% ir3 = 5*600
figure(6); clf
ph = .27
gap = (1-(3*ph))./5
px = [2*gap, 3*gap+ph, 4*gap+2*ph]

ax1=subplot(311);
pcolorjw(xf,y,z(:,:,is))
hold on
plot(xf,squeeze(peaklocs(:,:,is)),'.')
plot(xf,squeeze(mhw(:,is)),'.b')
plot(xf,squeeze(infleclocs(:,:,is)),'.')
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
set(gca, 'fontsize', 12)
%axis equal
xlim([100,1400])
ylim([0,250])
caxis([-8 14])
%set(gca,'xticklabels',[])
colormap(ax1,cmocean('-tarn'))

pos = get(gca, 'Position')
pos(2)=px(3);
pos(4)=ph;
set(gca, 'Position', pos);
set(gca, 'fontsize', 12);
colorbar
text(200,230,'Elevation [m NAVD88] 25-Jan-2017','fontsize',14)

ax2=subplot(312);
pcolorjw(xf,y,z(:,:,is+1))
hold on
plot([xf(ir) xf(ir)],[50, 150],'--r','linewidth',2)
plot([xf(ir2) xf(ir2)],[50, 150],'--r','linewidth',2)
plot([xf(ir3) xf(ir3)],[50, 150],'--r','linewidth',2)
h1=quiver(175,10,-50*sind(52),100*cosd(52),'color',[.2 .2 .2],'linewidth',2);
text(170,40,'North')
%axis equal
xlim([100,1400])
ylim([0,250])
%set(gca,'xticklabels',[])
pos = get(gca, 'Position')
pos(2)=px(2)
pos(4)=ph
set(gca, 'Position', pos)
set(gca, 'fontsize', 12)
caxis([-8 14])
colormap(ax2,cmocean('-tarn'))
ylabel('Cross-shore distance [m]','fontsize',14)
ht=text(200,230,'Elevation [m NAVD88] 14-Feb-2017','fontsize',14,'color','r');
colorbar

ax3=subplot(313);
v = [-2.5:.25:2.5];
%pcolorjw(xf,y,zd(:,:,is)
[C,hc]=contourf(xf,y,zdd,v,'linestyle','none');
hold on
%ts = [char(titles(is+1)),' minus ',char(titles(is))]
ts = 'Elevation change [m] 14-Feb minus 25-Jan-2017'
h1=quiver(175,10,-50*sind(52),100*cosd(52),'color',[.2 .2 .2],'linewidth',2);
text(170,40,'North')
ht=text(200,230,ts,'fontsize',14,'color','k');
xlabel('Alongshore distance [m]','fontsize',14)

%axis equal
xlim([100,1400])
ylim([0,250])
colormap(ax3,cmocean('-balance'));
pos = get(gca, 'Position')
pos(2)=px(1);
pos(4)=ph;
set(gca, 'Position', pos);
set(gca, 'fontsize', 12);
caxis([-2.5,2.5])
colorbar
print('compare_topo.png','-dpng','-r300')
%%
figure(7); clf
subplot(311)
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

subplot(312)
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

subplot(313)
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
