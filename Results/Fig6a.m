 % calculate ph and plot
  close all; clear; clc
  load vertical_summer_model
  Sz(1,:) = 0.0;
  
  load june2016
  dep104 = dep2(1,:);
  lon = reshape(lon,104,20);lon = lon';
  lat104 = lon(1,:); 
  dep_int = interp1(lat104,dep104,lat6,'linear','extrap');
  rr = dep_int./dep6;
  
  [pHx] = pH_cal(temp,sal,dic,ta,po4,sit);
  
  i1 = 2;  % time step
  load cp
  %%  image comparison
  name1 = pHx;   
  ph1 = squeeze(mean(name1(:,:,2),3))-0.10;
  figure(1);
  set(gcf,'position',[33 234 1138 735]);
  h1 = subplot(2,1,1);
  pcolor((lon),(-1*Sz.*dep2),(ph1));
  shading interp
  colorbar;
  set(gca,'ylim',([-30 0]),'fontweight','bold','fontsize',14);
  ylabel('Depth (M)','fontweight','bold');
  set(gca,'color',[0.2 0.2 0.2])
  hold on 
  k = 0;
  dep22 = -1*Sz.*dep2;
  
  n = length(lat6);
  for j = 1:n
      dif = abs(lat6(j)-lat104);
      tmp = sprintf('s%d',j);
      xx = eval(tmp);
      [m,n] = size(xx);
      xc = repmat(lat6(j),m,1);
      yc = -1*xx(:,1)*rr(j);
      qut = xx(:,9);
      sze = repmat(80,m,1);
      scatter(xc,yc,sze,qut,'filled'); 
      field_do(k+1:k+m) = qut;
      id = find(dif == min(dif));
      mod1 = interp1(dep22(:,id),do(:,id),yc,'linear','extrap');
      mod_do(k+1:k+m) = mod1;
      k = k+m;  
  end
  hold off;
  colormap(cmap3)
  caxis([7.0 8.4])
  set(gca,'xdir','reverse','xtick',[]);
  set(h1,'position',[0.1300    0.5838    0.7164    0.3412]);
  hh = colorbar;
  set(hh,'position',[0.8579    0.1987    0.0234    0.7235])
  text(37.4834,-21.3347,'\itJune pH','fontweight','bold','fontsize',14,'color','w')
  
%   figure(4*i+1);set(gcf,'position',[1153 461 700 519])
%   plot(field_do,mod_do,'bo');
%   a = polyfit(field_do,mod_do,1);
% 
%   xx = mean(field_do./mod_do)
% 
%   hh = refline(a(1),a(2)); set(hh,'color','r')
% %   set(gca,'xlim',[lim1(i) lim2(i)],'ylim',[lim1(i) lim2(i)]);
%    axis equal
%   xlabel('Field Measurement','fontweight','bold','fontsize',14)
%   ylabel('Model Simulation','fontweight','bold','fontsize',14)
% %   title([variable(i),' June 2016'],'fontweight','bold','fontsize',14);
%   r1 = refline(1,0); 
%   set(r1,'linestyle','--','color','k');
%   rmse = sqrt(sum((field_do-mod_do).^2)/length(field_do));  
%   rmse = floor(rmse*100)/100;
% 
%   xlim = get(gca,'xlim'); ylim = get(gca,'ylim');
%   text(0.1,0.9,['RMSE = ',num2str(rmse)],'fontweight','bold','Units','normalized','fontsize',14);
  
  load august2016
  dep_int = interp1(lat104,dep104,lat6,'linear','extrap');
  rr = dep_int./dep6;
  ph1 = squeeze(mean(name1(:,:,9),3));
  h2 = subplot(2,1,2);
  pcolor((lon),(-1*Sz.*dep2),(ph1));
  shading interp
  colorbar;
  set(gca,'ylim',([-30 0]),'fontweight','bold','fontsize',14);
  ylabel('Depth (M)','fontweight','bold');
  set(gca,'color',[0.2 0.2 0.2])
  hold on 
  k = 0;
  dep22 = -1*Sz.*dep2;
  
  n = length(lat6);
  for j = 1:n
      dif = abs(lat6(j)-lat104);
      tmp = sprintf('s%d',j);
      xx = eval(tmp);
      [m,n] = size(xx);
      xc = repmat(lat6(j),m,1);
      yc = -1*xx(:,1)*rr(j);
      qut = xx(:,9);
      sze = repmat(80,m,1);
      scatter(xc,yc,sze,qut,'filled'); 
      field_do(k+1:k+m) = qut;
      id = find(dif == min(dif));
      mod1 = interp1(dep22(:,id),do(:,id),yc,'linear','extrap');
      mod_do(k+1:k+m) = mod1;
      k = k+m;  
  end
  hold off;
  colormap(cmap3)
  caxis([7.0 8.4])
  set(gca,'xdir','reverse') 
  set(h2,'position',[0.1300    0.2000    0.7164    0.3412])
  colorbar off
  xlabel('Latitude (N)')
  
  text(39.448,2.6283,'\itRiver','fontweight','bold','fontsize',12)
  text(37.1478,2.6283,'\itOcean','fontweight','bold','fontsize',12)
  
  text(37.4834,-21.3347,'\itAugust pH','fontweight','bold','fontsize',14,'color','w')
  
  
  
  
  