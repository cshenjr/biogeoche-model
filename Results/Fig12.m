 % calculate ph and plot
  clear; clc; close all
  
  load coastline 
  load hori_ph_data_base
  xo = [-75.8161 -75.8079 -75.7831 -75.6095 -75.4236 -75.3657 -75.4401 -75.8161];
  yo = [39.8968 39.4543 39.2419 38.8938 38.8997 39.5723 39.9027 39.8968];
  in = inpolygon(coastline(:,1),coastline(:,2),xo,yo);
  coastline(in,1) = nan; coastline(in,2) = nan; 
  x1 = [-75.5486 -75.7148 -75.8358 -75.9718 -76.0021 -75.9517 ...
    -75.7803 -75.3621 -75.0648 -75.0598 -75.5486];
   y1 = [37.8207 37.5331 37.2839 37.0858 36.9260 36.3509 ...
    36.1527 36.1336 36.4978 37.9997 37.8207];
  in = inpolygon(x,y,x1,y1);
  xx = zeros(82,122); xx(in) = nan;
  fsm = ncread('grid_new.nc','mask_rho');
  aa = find(fsm==0);
  %% base case
  load cmap
  tt = tt-tt(1);
  ph = phx1(:,:,:);
  [a,b,c] = size(ph);
  x1 = zeros(a,b);
  for i = 1:c
      temp = squeeze(ph(:,:,i));
      x2 = temp <7.5;
      x1 = x1 + x2;
  end    
  r2 = x1./(c*ones(a,b));  
  r2 = r2(2:end-1,2:end-1);
  r2(aa) = nan;
  load cp
  figure(1);
  set(gcf,'position',[151         188        1617         794]);  
  hh1 = subplot(1,3,1);
  h1 = pcolor(x,y,r2+xx);  
  shading interp
%   set(h1,'edgecolor','none')
  hh = colorbar;
  set(hh,'ytick',(0:0.25:5))
  set(hh,'yticklabel',{'0%','25%','50%'})
  hold on
  plot(coastline(:,1),coastline(:,2),'k-')
  ylabel('Latitude (N)','fontweight','bold');
  xlabel('Longitude (W)','fontweight','bold')
  set(gca,'color',[0.8 0.8 0.8],'fontweight','bold')   
  set(gca,'xlim',[-77.3608 -75.6],'ylim',[36.8 39.6456],'fontweight','bold',...
    'fontsize',14)     
  colmap1 = cmap1('blue',75,15,0);
  colmap3 = cmap1('red',75,15,0);
  colormap([colmap1;flipud(colmap3)]);
  
  hold off 
  r2_base = r2+xx;
  text(-77.2357,39.5451,'\it1.0X Nutrient','fontweight','bold','fontsize',12)
  
%% 50 reduction
  load hori_ph_data_less
  tt = tt-tt(1);
  ph = phx1(:,:,:);
  [a,b,c] = size(ph);
  x1 = zeros(a,b);
  for i = 1:c
      temp = squeeze(ph(:,:,i));
      x2 = temp <7.5;
      x1 = x1 + x2;
  end    
  r2 = x1./(c*ones(a,b));  
  r2 = r2(2:end-1,2:end-1);
  r2(aa) = nan;
  hh2 = subplot(1,3,2);      
  h2 = pcolor(x,y,r2_base-(r2+xx));  
%   shading interp
  hh = colorbar;
  set(h2,'edgecolor','none')
  set(hh,'ylim',[-0.2 0.2],'ytick',[-0.2 -0.1 0.0 0.1 0.2])
  set(hh,'yticklabel',{'-20%','-10%','0%','10%','20%'})
  hold on
  plot(coastline(:,1),coastline(:,2),'k-')
  xlabel('Longitude (W)','fontweight','bold')
  title('% of annual duration (bottom pH <7.5) ','fontweight','bold','fontsize',20);
  set(gca,'color',[0.8 0.8 0.8],'fontweight','bold')   
  set(gca,'xlim',[-77.3608 -75.6],'ylim',[36.8 39.6456],'fontweight','bold',...
    'fontsize',14,'ytick',[])    
  colmap1 = cmap1('blue',75,15,0);
  colmap3 = cmap1('red',75,15,0);
  colormap([colmap1;flipud(colmap3)]);
  colorbar off
  hold off
  text(-77.2357,39.5451,'\it1.0X-0.5X','fontweight','bold','fontsize',12)
%% 50 increase
  load hori_ph_data_more
  tt = tt-tt(1);
  ph = phx1(:,:,:);
  [a,b,c] = size(ph);
  x1 = zeros(a,b);
  for i = 1:c
      temp = squeeze(ph(:,:,i));
      x2 = temp <7.5;
      x1 = x1 + x2;
  end    
  r2 = x1./(c*ones(a,b));  
  r2 = r2(2:end-1,2:end-1);
  r2(aa) = nan;
  hh3 = subplot(1,3,3);         
  h3 = pcolor(x,y,r2_base-(r2+xx));
  colormap(gcf,cmap3)    
%   shading interp
  hh = colorbar;
  set(h3,'edgecolor','none')
  set(hh,'ylim',[-0.2 0.2],'ytick',(-0.2:0.1:0.2))
  set(hh,'yticklabel',{'-20%','-10%','0%','10%','20%'})
  hold on
  plot(coastline(:,1),coastline(:,2),'k-')
  xlabel('Longitude (W)','fontweight','bold')
  set(gca,'color',[0.8 0.8 0.8],'fontweight','bold')   
  set(gca,'xlim',[-77.3608 -75.6],'ylim',[36.8 39.6456],'fontweight','bold',...
    'fontsize',14,'ytick',[])
  colmap1 = cmap1('blue',75,15,0);
  colmap3 = cmap1('red',75,15,0);
  colormap([colmap1;flipud(colmap3)]);
  hold off
  text(-77.2357,39.5451,'\it1.0X-1.5X','fontweight','bold','fontsize',12)
      
  set(hh1,'position',[0.2300    0.1100    0.2    0.8150]);
  set(hh2,'position',[0.4858    0.1100    0.2    0.8150]);
  set(hh3,'position',[0.6916    0.1100    0.2    0.8150]);
  
  