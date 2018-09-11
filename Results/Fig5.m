 % calculate ph and plot
  close all; clear; clc
  load vertical_model_result_p
  Sz(1,:) = 0.0;
  path = '../../calcium1';
  
  load month_data  
  dep104 = dep2(1,:);   % model depth for 104 vertical profiles
  lon = reshape(lon,104,20);lon = lon';  % latitude matrix
  lat104 = lon(1,:); 
  
  month = {'mar','apr','may','jun','jul','aug','sep','oct','nov','dec'};
  month1 = {'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
  n = length(month);
  load cp
  x1 = 0.8;
  figure(1);
  set(gcf,'position',[410         271        1063         726]);
  k = 0;
  for i = 1:2:7
     k = k+1;
     eval(['mon = ', month{i},';']);
     dep_int = interp1(lat104,dep104,mon(:,2),'linear','extrap');
     rr = dep_int./mon(:,9);  % depth ratio between model and field 
     %%  image comparison  
     gap = 12+i:14+i;
     h_dic = squeeze(mean(dic(:,:,gap),3));
     h_ta = squeeze(mean(ta(:,:,gap),3));
 
     h{2*k-1} = subplot(4,2,2*k-1);
     pcolor(lon,(-1*Sz.*dep2),h_dic);
     shading interp
     colorbar;
     set(gca,'ylim',([-30 0]),'fontweight','bold');
     set(gca,'color',[0.2 0.2 0.2])
     hold on
     yc = -1*mon(:,5).*rr;
     m = length(yc);
     sze = repmat(40,m,1);
     hh = scatter(mon(:,2),yc,sze,mon(:,7),'filled'); 
     hold off;
     colormap(h{2*k-1},cmap3)
     set(gca,'xdir','reverse') 
     caxis([800 2000])
     if 2*k-1<6
     set(gca,'xtick',[])
     end
     
     set(gca,'fontsize',12)
     
     if i == 1
         text(39.4351,-27.4315,'\itRiver','fontweight','bold','fontsize',11,'color','w');
         text(37.4003,-27.4315,'\itOcean','fontweight','bold','fontsize',11,'color','w');
     end
     
     xx = get(h{2*k-1},'position');
     xx(3) = xx(3)+0.03;
     xx(4) = xx(4)+0.03;
     set(h{2*k-1},'position',xx)
     colorbar off
     text(37.4886,-18.9286,['\it',month1{i},' DIC'],'fontweight','bold','fontsize',11,...
     'color','w') 
 
      if k ==4
         ylabel('Depth (m)')
         xlabel('Latitude (N)')
     end
 
     h{2*k} = subplot(4,2,2*k);
     pcolor((lon),(-1*Sz.*dep2),h_ta);
     shading interp
     colorbar;
     set(gca,'ylim',([-30 0]),'fontweight','bold');
     set(gca,'color',[0.2 0.2 0.2])
     hold on
     yc = -1*mon(:,5).*rr;
     m = length(yc);
     sze = repmat(40,m,1);
     scatter(mon(:,2),yc,sze,mon(:,8),'filled'); 
     hold off;
     colormap(h{2*k},cmap3)
     set(gca,'xdir','reverse') 
     caxis([800 2000]);
     if 2*k-1<6
     set(gca,'xtick',[])
     end
  
     set(gca,'ytick',[])
     set(gca,'fontsize',12)
     xx = get(h{2*k},'position');
     xx(1) = xx(1)-0.06;
     xx(3) = xx(3)+0.03;
     xx(4) = xx(4)+0.03;
     set(h{2*k},'position',xx)
     
     if k ==4
         xlabel('Latitude (N)')
     end
     
     text(37.4886,-18.9286,['\it',month1{i},' TA'],'fontweight','bold','fontsize',11,...
     'color','w')
     colorbar off
     
     if k == 4
         c1 = colorbar;
     end
     
  end
  
set(c1,'position',[0.90    0.13    0.0251    0.800])
  
text(36.8307,105.1103,'\itDIC/TA','fontweight','bold','fontsize',10)
text(36.6304,97.8309,'\it(µ mol/kg)','fontweight','bold','fontsize',8)
  
  
  
  
  
  
  
  
  
  
  
  
  