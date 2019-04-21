clear
%%
fnam='../grids/CCBay_MG1.nc';
nc=ncgeodataset(fnam);
lonm=nc{'lon_rho'}(:);
latm=nc{'lat_rho'}(:);
hm=nc{'h'}(:);
fnam='../grids/CCBay_CG1.nc';
nc=ncgeodataset(fnam);
lonc=nc{'lon_rho'}(:);
latc=nc{'lat_rho'}(:);
hc=nc{'h'}(:);
fnam='../grids/CCBay_FG1.nc';
nc=ncgeodataset(fnam);
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
h=nc{'h'}(:);

fnam='hsigF.nc';
nc=ncgeodataset(fnam);

%lon=nc{'longitude'}(:);
%lat=nc{'latitude'}(:);
%hs=double(nc{'hs'}(:));
time=double(nc{'time'}(:));time=time/24/3600+datenum(1970,01,01);
%%
hsa=load('hsigF.mat');
for i=1:length(time)
    [yy,mm,dd,hh,du,dum]=datevec(time(i));
    datt=datestr(time(i),'yyyymmdd_HHMMSS');
    eval(['hs(i,:,:)=hsa.Hsig_',datt,';']);
end
hs=double(hs);
%%
wda=load('wdirF.mat');
for i=1:length(time)
    [yy,mm,dd,hh,du,dum]=datevec(time(i));
    datt=datestr(time(i),'yyyymmdd_HHMMSS');
    eval(['wd(i,:,:)=wda.Dir_',datt,';']);
end
wd=(double(wd));
v=sind(-90-squeeze(wd));
u=cosd(-90-squeeze(wd));

wdac=load('wdirC.mat');
for i=1:length(time)
    [yy,mm,dd,hh,du,dum]=datevec(time(i));
    datt=datestr(time(i),'yyyymmdd_HHMMSS');
    eval(['wdc(i,:,:)=wdac.Dir_',datt,';']);
end
wdc=(double(wdc));
vc=sind(-90-squeeze(wdc));
uc=cosd(-90-squeeze(wdc));
%%

figure(1);clf;
for i=1:length(time)
    pcolor(lon,lat,squeeze(hs(i,:,:)));shading interp;hold on
    caxis([0,3]);colorbar;dasp(mean(lat(:)));
    quiver(lon(1:8:end,1:8:end),lat(1:8:end,1:8:end),u(1:8:end,1:8:end),v(1:8:end,1:8:end));
    drawnow;%pause(0.01);
    title('Hs (m)')
    print -dpng -painters Hs_stat_fine.png

end
%%
dia=load('dissipF.mat');
for i=1:length(time)
    [yy,mm,dd,hh,du,dum]=datevec(time(i));
    datt=datestr(time(i),'yyyymmdd_HHMMSS');
    eval(['diss(i,:,:)=dia.Dissip_',datt,';']);
end
diss=double(diss);

figure(2);clf;
for i=1:length(time)
    pcolor(lon,lat,squeeze(diss(i,:,:)));shading interp
    caxis([0,400]);colorbar
    drawnow;%pause;
    title('Dissipation')
    print -dpng -painters Diss_stat_fine.png
end
%%
hsam=load('hsigM.mat');
for i=1:length(time)
    datt=datestr(time(i),'yyyymmdd_HHMMSS');
    eval(['hsm(i,:,:)=hsam.Hsig_',datt,';']);
end
hsm=double(hsm);
hsac=load('hsigC.mat');
for i=1:length(time)
    datt=datestr(time(i),'yyyymmdd_HHMMSS');
    eval(['hsc(i,:,:)=hsac.Hsig_',datt,';']);
end
hsc=double(hsc);

figure(4);clf;
for i=1%:length(time)
    pcolor(lonc,latc,squeeze(hsc(i,:,:)));shading interp;hold on
    pcolor(lonm,latm,squeeze(hsm(i,:,:)));shading interp;
    pcolor(lon,lat,squeeze(hs(i,:,:)));shading interp
    caxis([0,4]);colorbar
    quiver(lonc(1:8:end,1:8:end),latc(1:8:end,1:8:end),uc(1:8:end,1:8:end),vc(1:8:end,1:8:end));
    contour(lonc,latc,hc,[0,0]);
    contour(lonm,latm,hm,[0,0]);
    contour(lon,lat,h,[0,0]);
    drawnow;%pause;
    title('Hs (m)')
    print -dpng -painters Hs_stat_all.png
end






