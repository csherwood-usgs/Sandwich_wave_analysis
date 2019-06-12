% scatterplot_anal - Look a cross correlations
% first:
% 1) run analyze_topo_maps
% 2) run capecodbay_anal



% xg
% gnang
Pbx = Pb./1000 .* dalongb.*gnang;
Pdissx = Pdiss.*dalongdiss.*gnang;
dPbx = [NaN; -diff( Pb.*dalongb./1000/dx).*gnang(1:end-1)];
% slopediss
% slopeb
% hsdiss
% hsb
% R2b
%
% xf
% gnan
dvol = dall_vols(:,7);
ddvol = ddune_vols(:,7);
% mwl_rate
% mhw_rate
% mhhw_rate

%% bin variables on 1-m topo grid by bsize-m segments
bsize = 35
nb = fix( 1300/bsize)

gn = ones(length(xf),1);
% if you want to omit data near groins, uncomment next line
gn = gnan;

xf_b = nan*ones(nb,1);
dvol_b = nan*ones(nb,1);
dvol_sb = nan*ones(nb,1);
ddvol_b = nan*ones(nb,1);
ddvol_sb = nan*ones(nb,1);
mwlr_b = nan*ones(nb,1);
mwlr_sb = nan*ones(nb,1);
mhwr_b = nan*ones(nb,1);
mhwr_sb = nan*ones(nb,1);
mhhwr_b = nan*ones(nb,1);
mhhwr_sb = nan*ones(nb,1);
for i=1:nb
   iok = (xf>=(bsize-1)*i & xf<bsize*i);
   xfn = xf';
   xf_b(i) = nanmean(xfn(iok));
   dvoln = dvol.*gn;
   dvol_b(i) = nanmean(dvoln(iok));
   dvol_sb(i) = nanstd(dvoln(iok));
   ddvoln = ddvol.*gn;
   ddvol_b(i) = nanmean(ddvoln(iok));
   ddvol_sb(i) = nanstd(ddvoln(iok));
   mwl_raten = mwl_rate.*gn;
   mwlr_b(i) = nanmean(mwl_raten(iok));
   mwlr_bs(i) = nanstd(mwl_raten(iok));
   mhw_raten = mhw_rate.*gn;
   mhwr_b(i) = nanmean(mhw_raten(iok));
   mhwr_bs = nanstd(mhw_raten(iok));
   mhhw_raten = mhw_rate.*gn;
   mhhwr_b = nanmean(mhhw_raten(iok));
   mhhwr_sb = nanstd(mhhw_raten(iok));
end
%% bin variables on 5-m SWAN grid by bsize-m segments
gn = ones(length(xg),1);
% if you want to omit data near groins, uncomment next line
gn = gnang;

xg_b = nan*ones(nb,1);
Pb_b = nan*ones(nb,1);
Pb_sb = nan*ones(nb,1);
Pbx_b = nan*ones(nb,1);
Pbx_sb = nan*ones(nb,1);
dPbx_b = nan*ones(nb,1);
dPbx_sb = nan*ones(nb,1);
hsb_b = nan*ones(nb,1);
hsb_sb = nan*ones(nb,1);
R2b_b = nan*ones(nb,1);
R2b_sb = nan*ones(nb,1);
slopeb_b = nan*ones(nb,1);
slopeb_sb = nan*ones(nb,1);
for i=1:nb
   iok = (xg>=(bsize-1)*i & xg<bsize*i);
   xgn = xg';
   xg_b(i) = nanmean(xgn(iok));
   Pbn = Pb.*gn;
   Pb_b(i) = nanmean(Pbn(iok));
   Pb_sb(i) = nanstd(Pbn(iok));
   Pbxn = Pbx.*gn;
   Pbx_b(i) = nanmean(Pbxn(iok));
   Pbx_sb(i) = nanstd(Pbxn(iok));
   dPbxn = dPbx.*gn;
   dPbx_b(i) = nanmean(dPbxn(iok));
   dPbx_sb(i) = nanstd(dPbxn(iok));
   R2bxn = R2b.*gn;
   R2b_b(i) = nanmean(R2bxn(iok));
   R2b_sb(i) = nanstd(R2bxn(iok));
   hsbxn = hsb.*gn;
   hsb_b(i) = nanmean(hsbxn(iok));
   hsb_sb(i) = nanstd(hsbxn(iok));
end
isarmor = (xg_b<=400);
isdune = (xg_b>400 & xg_b<=800);
iseast = (xg_b>800);
%%
figure(21);clf
% Pbx v total vol change
subplot(221)
for i=1:length(Pbx_b)
plot([Pbx_b(i)-Pbx_sb(i),Pbx_b(i)+Pbx_sb(i)],[dvol_b(i), dvol_b(i)],'-','color',[.6 .6 .6])
hold on
plot([Pbx_b(i), Pbx_b(i)],[dvol_b(i)-dvol_sb(i), dvol_b(i)+dvol_sb(i)],'-','color',[.6 .6 .6])
end
h1=plot(Pbx_b(isarmor),dvol_b(isarmor),'.k','markersize',25);
h2=plot(Pbx_b(isdune),dvol_b(isdune),'.r','markersize',25);
h3=plot(Pbx_b(iseast),dvol_b(iseast),'.b','markersize',25);
ylabel('Total Volume Change [m^3/m]','fontsize',14)
text(.02,.95,'a','Fontsize',14,'Units','normalized')
LM = fitlm(Pbx_b(isdune),dvol_b(isdune))
ts = sprintf('{r}^2 = %.3f',LM.Rsquared.Ordinary)
text(.8,.1,ts,'Fontsize',12,'Units','normalized')
ts = sprintf('{p} = %.3f',LM.Coefficients.pValue(2))
text(.8,.05,ts,'Fontsize',12,'Units','normalized')
if(LM.Rsquared.Ordinary >= 0.3)
plot(Pbx_b(isdune),LM.Fitted,'--r','linewidth',2)
end


% dPbx v total vol change
subplot(222)
for i=1:length(dPbx_b)
plot([dPbx_b(i)-dPbx_sb(i),dPbx_b(i)+dPbx_sb(i)],[dvol_b(i), dvol_b(i)],'-','color',[.6 .6 .6])
hold on
plot([dPbx_b(i), dPbx_b(i)],[dvol_b(i)-dvol_sb(i), dvol_b(i)+dvol_sb(i)],'-','color',[.6 .6 .6])
end
h1=plot(dPbx_b(isarmor),dvol_b(isarmor),'.k','markersize',25);
h2=plot(dPbx_b(isdune),dvol_b(isdune),'.r','markersize',25);
h3=plot(dPbx_b(iseast),dvol_b(iseast),'.b','markersize',25);
text(.02,.95,'b','Fontsize',14,'Units','normalized')
LM = fitlm(dPbx_b(isdune),dvol_b(isdune))
ts = sprintf('{r}^2 = %.3f',LM.Rsquared.Ordinary)
text(.8,.1,ts,'Fontsize',12,'Units','normalized')
ts = sprintf('{r}^2 = %.3f',LM.Rsquared.Ordinary)
text(.8,.1,ts,'Fontsize',12,'Units','normalized')
ts = sprintf('{p} = %.3f',LM.Coefficients.pValue(2))
text(.8,.05,ts,'Fontsize',12,'Units','normalized')
if(LM.Rsquared.Ordinary >= 0.3)
plot(dPbx_b(isdune),LM.Fitted,'--r','linewidth',2)
end

% Pbx v dune vol change
subplot(223)
for i=1:length(Pbx_b)
plot([Pbx_b(i)-Pbx_sb(i),Pbx_b(i)+Pbx_sb(i)],[ddvol_b(i), ddvol_b(i)],'-','color',[.6 .6 .6])
hold on
plot([Pbx_b(i), Pbx_b(i)],[ddvol_b(i)-ddvol_sb(i), ddvol_b(i)+ddvol_sb(i)],'-','color',[.6 .6 .6])
end
h1=plot(Pbx_b(isarmor),ddvol_b(isarmor),'.k','markersize',25);
h2=plot(Pbx_b(isdune),ddvol_b(isdune),'.r','markersize',25);
h3=plot(Pbx_b(iseast),ddvol_b(iseast),'.b','markersize',25);
ylabel('Dune Volume Change [m^3/m]','fontsize',14)
xlabel('{\itP_x}  [kW m^{-1}]','fontsize',14)
LM = fitlm(Pbx_b(isdune),ddvol_b(isdune))
ts = sprintf('{r}^2 = %.3f',LM.Rsquared.Ordinary)
text(.8,.1,ts,'Fontsize',12,'Units','normalized')
ts = sprintf('{p} = %.3f',LM.Coefficients.pValue(2))
text(.8,.05,ts,'Fontsize',12,'Units','normalized')
if(LM.Rsquared.Ordinary >= 0.3)
plot(Pbx_b(isdune),LM.Fitted,'--r','linewidth',2)
end
text(.02,.95,'c','Fontsize',14,'Units','normalized')

% dPbx v dune vol change
subplot(224)
for i=1:length(dPbx_b)
plot([dPbx_b(i)-dPbx_sb(i),dPbx_b(i)+dPbx_sb(i)],[ddvol_b(i), ddvol_b(i)],'-','color',[.6 .6 .6])
hold on
plot([dPbx_b(i), dPbx_b(i)],[ddvol_b(i)-ddvol_sb(i), ddvol_b(i)+ddvol_sb(i)],'-','color',[.6 .6 .6])
end
h1=plot(dPbx_b(isarmor),ddvol_b(isarmor),'.k','markersize',25);
h2=plot(dPbx_b(isdune),ddvol_b(isdune),'.r','markersize',25);
h3=plot(dPbx_b(iseast),ddvol_b(iseast),'.b','markersize',25);
xlabel('-{\Delta} {\itP_x}/{\Delta}{\itx}  [kW m^{-2}]','fontsize',14)
text(.02,.95,'d','Fontsize',14,'Units','normalized')
LM = fitlm(dPbx_b(isdune),ddvol_b(isdune))
ts = sprintf('{r}^2 = %.3f',LM.Rsquared.Ordinary)
text(.8,.1,ts,'Fontsize',12,'Units','normalized')
ts = sprintf('{p} = %.3f',LM.Coefficients.pValue(2))
text(.8,.05,ts,'Fontsize',12,'Units','normalized')
if(LM.Rsquared.Ordinary >= 0.3)
plot(dPbx_b(isdune),LM.Fitted,'--r','linewidth',2)
end

%%
figure(22);clf
% R2b v total vol change
subplot(221)
for i=1:length(R2b_b)
plot([R2b_b(i)-R2b_sb(i),R2b_b(i)+R2b_sb(i)],[dvol_b(i), dvol_b(i)],'-','color',[.6 .6 .6])
hold on
plot([R2b_b(i), R2b_b(i)],[dvol_b(i)-dvol_sb(i), dvol_b(i)+dvol_sb(i)],'-','color',[.6 .6 .6])
end
h1=plot(R2b_b(isarmor),dvol_b(isarmor),'.k','markersize',25);
h2=plot(R2b_b(isdune),dvol_b(isdune),'.r','markersize',25);
h3=plot(R2b_b(iseast),dvol_b(iseast),'.b','markersize',25);
ylabel('Total Volume Change [m^3/m]','fontsize',14)
text(.02,.95,'a','Fontsize',14,'Units','normalized')
LM = fitlm(R2b_b(isdune),dvol_b(isdune))
ts = sprintf('{r}^2 = %.3f',LM.Rsquared.Ordinary)
text(.8,.1,ts,'Fontsize',12,'Units','normalized')
ts = sprintf('{p} = %.3f',LM.Coefficients.pValue(2))
text(.8,.05,ts,'Fontsize',12,'Units','normalized')
if(LM.Rsquared.Ordinary >= 0.3)
plot(R2b_b(isdune),LM.Fitted,'--r','linewidth',2)
end

% Hsb v total vol change
subplot(222)
for i=1:length(hsb_b)
plot([hsb_b(i)-hsb_sb(i),hsb_b(i)+hsb_sb(i)],[dvol_b(i), dvol_b(i)],'-','color',[.6 .6 .6])
hold on
plot([hsb_b(i), hsb_b(i)],[dvol_b(i)-dvol_sb(i), dvol_b(i)+dvol_sb(i)],'-','color',[.6 .6 .6])
end
h1=plot(hsb_b(isarmor),dvol_b(isarmor),'.k','markersize',25);
h2=plot(hsb_b(isdune),dvol_b(isdune),'.r','markersize',25);
h3=plot(hsb_b(iseast),dvol_b(iseast),'.b','markersize',25);
text(.02,.95,'b','Fontsize',14,'Units','normalized')
LM = fitlm(hsb_b(isdune),dvol_b(isdune))
ts = sprintf('{r}^2 = %.3f',LM.Rsquared.Ordinary)
text(.8,.1,ts,'Fontsize',12,'Units','normalized')
ts = sprintf('{r}^2 = %.3f',LM.Rsquared.Ordinary)
text(.8,.1,ts,'Fontsize',12,'Units','normalized')
ts = sprintf('{p} = %.3f',LM.Coefficients.pValue(2))
text(.8,.05,ts,'Fontsize',12,'Units','normalized')
if(LM.Rsquared.Ordinary >= 0.3)
plot(hsb_b(isdune),LM.Fitted,'--r','linewidth',2)
end

% R2b v dune vol change
subplot(223)
for i=1:length(R2b_b)
plot([R2b_b(i)-R2b_sb(i),R2b_b(i)+R2b_sb(i)],[ddvol_b(i), ddvol_b(i)],'-','color',[.6 .6 .6])
hold on
plot([R2b_b(i), R2b_b(i)],[ddvol_b(i)-ddvol_sb(i), ddvol_b(i)+ddvol_sb(i)],'-','color',[.6 .6 .6])
end
h1=plot(R2b_b(isarmor),ddvol_b(isarmor),'.k','markersize',25);
h2=plot(R2b_b(isdune),ddvol_b(isdune),'.r','markersize',25);
h3=plot(R2b_b(iseast),ddvol_b(iseast),'.b','markersize',25);
ylabel('Dune Volume Change [m^3/m]','fontsize',14)
xlabel('{\itR}_2 [m]','fontsize',14)
LM = fitlm(R2b_b(isdune),ddvol_b(isdune))
ts = sprintf('{r}^2 = %.3f',LM.Rsquared.Ordinary)
text(.8,.1,ts,'Fontsize',12,'Units','normalized')
ts = sprintf('{p} = %.3f',LM.Coefficients.pValue(2))
text(.8,.05,ts,'Fontsize',12,'Units','normalized')
if(LM.Rsquared.Ordinary >= 0.3)
plot(R2b_b(isdune),LM.Fitted,'--r','linewidth',2)
end
text(.02,.95,'c','Fontsize',14,'Units','normalized')

% hsb v dune vol change
subplot(224)
for i=1:length(hsb_b)
plot([hsb_b(i)-hsb_sb(i),hsb_b(i)+hsb_sb(i)],[ddvol_b(i), ddvol_b(i)],'-','color',[.6 .6 .6])
hold on
plot([hsb_b(i), hsb_b(i)],[ddvol_b(i)-ddvol_sb(i), ddvol_b(i)+ddvol_sb(i)],'-','color',[.6 .6 .6])
end
h1=plot(hsb_b(isarmor),ddvol_b(isarmor),'.k','markersize',25);
h2=plot(hsb_b(isdune),ddvol_b(isdune),'.r','markersize',25);
h3=plot(hsb_b(iseast),ddvol_b(iseast),'.b','markersize',25);
xlabel('{\itH_s} [m]','fontsize',14)
text(.02,.95,'d','Fontsize',14,'Units','normalized')
LM = fitlm(hsb_b(isdune),ddvol_b(isdune))
ts = sprintf('{r}^2 = %.3f',LM.Rsquared.Ordinary)
text(.8,.1,ts,'Fontsize',12,'Units','normalized')
ts = sprintf('{p} = %.3f',LM.Coefficients.pValue(2))
text(.8,.05,ts,'Fontsize',12,'Units','normalized')
if(LM.Rsquared.Ordinary >= 0.3)
plot(hsb_b(isdune),LM.Fitted,'--r','linewidth',2)
end
