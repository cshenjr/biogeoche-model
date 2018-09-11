% axis along points for vertical profile 
close all;
clear;
clc;

load coastline 
xo = [-75.8161 -75.8079 -75.7831 -75.6095 -75.4236 -75.3657 -75.4401 -75.8161];
yo = [39.8968 39.4543 39.2419 38.8938 38.8997 39.5723 39.9027 39.8968];
in = inpolygon(coastline(:,1),coastline(:,2),xo,yo);
coastline(in,1) = nan; coastline(in,2) = nan; 
xo = [-75.5445 -75.7026 -75.8081 -75.9102 -75.9432 -75.9629 -75.9695 -75.9036 ...
    -75.7916 -75.6697 -75.5972 -75.5379 -75.5049 -75.4984 -75.5445];
yo = [37.9355 37.7064 37.5491 37.3829 37.3065 37.2031 37.0998 37.0504 ...
    37.1088 37.2436 37.4592 37.6300 37.7558 37.8726 37.9085];
in = inpolygon(coastline(:,1),coastline(:,2),xo,yo);
coastline(in,1) = nan; coastline(in,2) = nan; 

load bathy
rt = 1;
lat = lat(1:rt:end);
lon = lon(1:rt:end);
dep = H(1:rt:end,1:rt:end);
% [a,b] = meshgrid(lon,lat);

id = (dep>=-0.5);
dep(id) = nan;
id = (dep<=-45);
dep(id) = nan;

load abc
%%  may
h1 = figure(1);
set(gcf,'position',[615    41   783   956]);

h = pcolor(lon,lat,dep);
set(h,'edgecolor','none')
shading interp;
c1 = colorbar;caxis([-40 0])
set(c1,'location','southoutside','direction','reverse');
set(c1,'position',[0.1403    0.8662    0.3052    0.0200]);
set(c1,'ytick',[-40:10:0],'fontsize',10,'yticklabel',{'40  (m)' '30' '20' '10' '0'})

hold on
plot(coastline(1:end,1),coastline(1:end,2),'-','color',[0.1 0.1 0.1])
colormap(h1,ab);
plot([-77.0316,-76.6457] ,[38.9242,39.2675],'kp','markersize',14)

% xs = [-76.02599 -76.17579 -76.24050 -76.30631 -76.38000 -76.39445 -76.42127 -76.42794 -76.29215 -76.22787 ...
%     -76.17137 -76.15633 -76.20799 -76.17466 -76.16216 -76.15966 -76.02048 -75.85692];
% ys = [39.44149 39.34873 39.24950 39.16369 38.96210 38.82593 38.64618 38.55505 38.31870 38.13705 37.91011 ...
%     37.48680 37.23653 37.80013 37.58847 37.41153 36.99570 37.06183];
xs = [-76.02599 -76.17579 -76.30631 -76.38000 -76.39445 -76.42794 -76.29215 -76.22787 ...
    -76.17137 -76.15633 -76.20799 -75.85692];
ys = [39.44149 39.34873  39.16369 38.96210 38.82593 38.55505 38.31870 38.13705 37.91011 ...
    37.48680 37.23653 37.06183];
hh1 = plot(xs,ys,'ko','linewidth',1.5);

plot(-77.2462,39.3406,'ko','linewidth',1.5);
text(-77.1962,39.3406,'Field Stations','fontweight','bold','fontsize',12);
% plot([-77.2810 -77.2404],[39.2529 39.2529],'r-','linewidth',3)
% text(-77.1962,39.2529,'Freshwater Boundary','fontweight','bold','fontsize',12);

plot([-76.6922,-76.0424],[38.9936,38.7890],'k--','linewidth',3)
plot([-76.5442,-75.7117],[37.900,37.9],'k--','linewidth',3)

% plot([-76.1381 -76.0801],[39.5671 39.6182],'k-','linewidth',3)
% plot([-76.5993 -76.6022],[39.2200 39.2894],'k-','linewidth',3)
% plot([-77.0635 -76.9938],[38.8547 38.8547],'k-','linewidth',3)
% plot([-76.8923 -76.9271],[38.0219 37.9634],'k-','linewidth',3)
% plot([-76.7676 -76.7995],[37.5434 37.4886],'k-','linewidth',3)
% plot([-76.7096 -76.6486],[38.5990 38.5990],'r-','linewidth',3)
% plot([-77.1766 -77.1766],[37.2657 37.3498],'r-','linewidth',3)
% plot([-75.9757 -75.9292],[38.6685 38.6685],'r-','linewidth',3)

hold off
set(gca,'color',[0.8 0.8 0.8]);
set(gca,'xlim',[-77.3608 -75.6],'ylim',[36.8 39.6456],'fontweight','bold',...
    'fontsize',14,'ytick',[37 38 39],'xtick',[-77 -76])
set(gca,'ticklength',[0 0])

text(-77.1157,39.0082,'Washington D.C','fontweight','bold','fontsize',12)
text(-76.7357,39.3296,'Baltimore','fontweight','bold','fontsize',12)
xlabel('Longitude (W)');ylabel('Latitude (N)');

text(-77.2120,37.2439,'\itJames R.','fontweight','bold','rotation',-15);
text(-76.8720,37.4667,'\itYork R.','fontweight','bold','rotation',-30);
text(-77.1813,38.1245,'\itRappahannock R.','fontweight','bold','rotation',-30);
text(-77.0527,38.6434,'\itPotomac R.','fontweight','bold','rotation',-50);
text(-76.7492,38.7086,'\itPauxent R.','fontweight','bold','rotation',-75);
text(-76.4444,39.5889,'\itSusquehanna R.','fontweight','bold');
text(-75.9670,39.1726,'\itChester R.','fontweight','bold','rotation',20);
text(-75.8857,38.7452,'\itChoptank R.','fontweight','bold','rotation',45);
text(-75.8973,38.4164,'\itNanticoke R.','fontweight','bold','rotation',50);

aa = 12;
text(-76.1311,39.4063,'CB2.1','fontweight','bold','fontsize',aa);
text(-76.3284 ,  39.3114,'CB2.2','fontweight','bold','fontsize',aa);
% text(-76.2193  , 39.2091,'CB3.1','fontweight','bold','fontsize',aa);
text(-76.3782 ,  39.2108,'CB3.2','fontweight','bold','fontsize',aa);
text(-76.390 ,  39.0168,'858','fontweight','bold','fontsize',aa);
text(-76.5014 ,  38.7853,'CB4.1C','fontweight','bold','fontsize',aa);
% text(-76.4202 ,  38.6246,'CB4.2C','fontweight','bold','fontsize',aa);
text(-76.4985 ,  38.5077,'CB4.3C','fontweight','bold','fontsize',aa);
text(-76.3622 ,  38.2776,'CB5.1','fontweight','bold','fontsize',aa);
text(-76.2396 ,  38.0803,'CB5.2','fontweight','bold','fontsize',aa);
text(-76.1623 ,  37.8685,'CB5.3','fontweight','bold','fontsize',aa);
% text(-76.1584,   37.7479,'CB5.4','fontweight','bold','fontsize',aa);
% text(-76.1265  , 37.5945,'CB6.1','fontweight','bold','fontsize',aa);
text(-76.1178  , 37.4776,'CB6.2','fontweight','bold','fontsize',aa);
% text(-76.1265  , 37.3899,'CB6.3','fontweight','bold','fontsize',aa);
text(-76.1584  , 37.2146,'CB6.4','fontweight','bold','fontsize',aa);
% text(-76.0685  , 37.0393,'CB7.4','fontweight','bold','fontsize',aa);
text(-75.8509  , 37.0064,'AO1','fontweight','bold','fontsize',aa);

text(-76.5529 ,39.4684,'\itUpper Bay','fontweight','bold','fontsize',12);
text(-76.6825  , 38.3039,'\itMid Bay','fontweight','bold','fontsize',12);
text( -76.6861  , 37.5020,'\itLower Bay','fontweight','bold','fontsize',12);














