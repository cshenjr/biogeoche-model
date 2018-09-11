% axis along points for vertical profile 
close all;
clear;
clc;

load coastline 
xo = [-75.8161 -75.8079 -75.7831 -75.6095 -75.4236 -75.3657 -75.4401 -75.8161];
yo = [39.8968 39.4543 39.2419 38.8938 38.8997 39.5723 39.9027 39.8968];
in = inpolygon(coastline(:,1),coastline(:,2),xo,yo);
coastline(in,1) = nan; coastline(in,2) = nan; 

load pco2_model
load pco2
lat = ncread('grid_new.nc','lat_rho');
lon = ncread('grid_new.nc','lon_rho');
x1 = [-75.5486 -75.7148 -75.8358 -75.9718 -76.0021 -75.9517 ...
    -75.7803 -75.3621 -75.0648 -75.0598 -75.5486];
y1 = [37.8207 37.5331 37.2839 37.0858 36.9260 36.3509 ...
    36.1527 36.1336 36.4978 37.9997 37.8207];
in = inpolygon(lon,lat,x1,y1);
xx = zeros(82,122); xx(in) = nan;
% column 3 lon 4 lat, 12 temperature 14 pco2 presure
bry_time = ncread('Y2016_eutr_bry.nc','bry_time');
bt = bry_time(2:end)-15;
may1 = find(tt>=bt(5) & tt<bt(6));
jun1 = find(tt>=bt(6) & tt<bt(7));
aug1 = find(tt>=bt(8) & tt<bt(9));
oct1 = find(tt>=bt(9) & tt<bt(10));
%oct1 = find(tt>=bt(10) & tt<bt(11));
%%  may
figure(1);
set(gcf,'position',[787    75   635   919]);
h11 = subplot(2,2,1)
pco21 = mean(pco2(2:end-1,2:end-1,may1),3);
pco22 = pco21 + xx;  
pcolor(lon,lat,pco22);
shading interp
% title('\it pCO2 May','fontweight','bold')
hold on
plot(coastline(:,1),coastline(:,2),'k-')
colormap jet;caxis([200 800]);h1 = colorbar;
hold off
set(gca,'color',[0.8 0.8 0.8]);
set(gca,'xlim',[-77.3608 -75.6],'ylim',[36.8 39.6456],'fontweight','bold',...
    'fontsize',12,'ytick',[37 38 39])
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);

set(h11,'position',[0.1300 0.56 0.38 0.37])

set(h1,'location','southoutside');
set(h1,'position',[0.1307    0.5133    0.78    0.015])

text(-77.1844,39.5489,'\itSusquehanna River','fontweight','bold');
text(-77.3256,39.3066,'\itpCO2 May','fontweight','bold');
%--------------------



%% june
h22 = subplot(2,2,2)
pco21 = mean(pco2(2:end-1,2:end-1,jun1),3);
pco22 = pco21 + xx;  
pcolor(lon,lat,pco22);
shading interp
% title('\it pCO2 June','fontweight','bold')
hold on
plot(coastline(:,1),coastline(:,2),'k-')
colormap jet;caxis([200 800]);
hold off
set(gca,'color',[0.8 0.8 0.8]);
set(gca,'xlim',[-77.3608 -75.6],'ylim',[36.8 39.6456],'fontweight','bold',...
    'fontsize',12,'ytick',[37 38 39])
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(h22,'position',[0.53 0.56 0.38 0.37])

text(-77.1844,39.5489,'\itSusquehanna River','fontweight','bold');
text(-77.3256,39.3066,'\itpCO2 June','fontweight','bold');
%--------------------


%% august
h33 = subplot(2,2,3)
pco21 = mean(pco2(2:end-1,2:end-1,aug1),3);
pco22 = pco21 + xx;  
pcolor(lon,lat,pco22);
shading interp
% title('\it pCO2 August','fontweight','bold')
hold on
plot(coastline(:,1),coastline(:,2),'k-')
colormap jet;caxis([200 800]);
hold off
set(gca,'color',[0.8 0.8 0.8]);
set(gca,'xlim',[-77.3608 -75.6],'ylim',[36.8 39.6456],'fontweight','bold',...
    'fontsize',12,'ytick',[37 38 39])
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(h33,'position',[0.1300 0.11 0.38 0.37])
text(-77.1844,39.5489,'\itSusquehanna River','fontweight','bold');
text(-77.3256,39.3066,'\itpCO2 August','fontweight','bold');
%--------------------


%% october
h44 = subplot(2,2,4)
pco21 = mean(pco2(2:end-1,2:end-1,oct1),3);
pco22 = pco21 + xx;  
pcolor(lon,lat,pco22);
shading interp
% title('\it pCO2 October','fontweight','bold')
hold on
plot(coastline(:,1),coastline(:,2),'k-')
colormap jet;caxis([200 800]);
hold off
set(gca,'color',[0.8 0.8 0.8]);
set(gca,'xlim',[-77.3608 -75.6],'ylim',[36.8 39.6456],'fontweight','bold',...
    'fontsize',12,'ytick',[37 38 39])
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(h44,'position',[0.5300 0.11 0.38 0.37])
text(-77.1844,39.5489,'\itSusquehanna River','fontweight','bold');
text(-77.3256,39.3066,'\itpCO2 October','fontweight','bold');
%--------------------



