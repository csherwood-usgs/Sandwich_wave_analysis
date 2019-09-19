% Breakpoint formula Ashton, Muarray, Arnoult (2001)
clear
rhow = 1030.;
g = 9.81;
K = 6.4e4 % m^3/day

dx = 5;
dy = 5;

mrun(1).datt = '20170213_000000'
mrun(1).wlev = 2.21
mrun(1).rundir = '../swan_stationary/H2017_0d_2p21_stat/'
mrun(1).rname = '2017 3p44m, 0\circ, +2.21 m, FG2'
mrun(1).fnam='../swan_stationary/grids/CCBay_FG2.nc';
mrun(1).pname = '2017_3p44m_0d_2p21_FG2'

mrun(2).datt = '20170213_000000'
mrun(2).wlev = 2.21
mrun(2).rundir = '../swan_stationary/H2017_0d_2p21_stat4/'
mrun(2).rname = '2017 3p44m, 0\circ, +2.21 m, FG4'
mrun(2).fnam='../swan_stationary/grids/CCBay_FG4.nc';
mrun(2).pname = '2017_3p44m_0d_2p21_FG4'

for i=1:length(mrun)
   
   rundir = mrun(i).rundir
   datt = mrun(i).datt
   fnam = mrun(i).fnam
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
   for ii=1:ncols
      [R2b(ii),S,setup, Sinc, SIG, ir, R16] = calcR2(Hob(ii),Tpb(ii),atan(slopedry(ii)),0);
      [R2diss(ii),S,setup, Sinc, SIG, ir, R16] = calcR2(Hodiss(ii),Tpdiss(ii),atan(slopedry(ii)),0);
      
   end
   %% calculate Iribarren number at breakpoint according to Battjes (1974); Komar (1988, p. 210) Holthuijsen (2007) and Wikipedia
   for ii=1:ncols
      Lo(ii) = (g/(2.*pi))*Tpb(i).^2;
      eb(ii) = slopeb(ii)./( sqrt( hsb(ii)/Lo(ii) ) );
      eo(ii) = slopeb(ii)./( sqrt( Hob(ii)/Lo(ii) ) );
      % surf-scaling paramter (Guza and Inman, 1975)
      ep(ii) = 0.5*hsb(ii).*(2.*pi/Tpb(ii)).^2/(g*slopeb(ii).^2);
   end
   
   % assign these to the model structure
   varlist = {'ibr','idiss','idry','dissmx','slopedry','Tpb','Tpdiss','slopeb','slopediss',...
      'hsb','hsdiss','khb','khdiss','hb','hdiss','Pb','Pdiss','dirb','dirdiss','R2b','R2diss','eb','eo','ep'}
   for iv = 1:length(varlist)
      eval( char( strcat('mrun(i).',varlist(iv),'=',varlist(iv),';')))
   end
   % ...then clear them
   for iv = 1:length(varlist)
      disp( char( strcat(''clear '','',varlist(iv),'')))
   end
end

theta = 0;
phib = [0:5:90.]'
Hb = 2.
Qs = K*Hb.^(5/2)*sind(phib-theta).*cosd(phib-theta)

plot(phib,Qs)
