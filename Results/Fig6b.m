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
  %%  image comparison
  name1 = pHx;   
  ph1 = squeeze(mean(name1(:,:,2),3))-0.1;
  dep22 = -1*Sz.*dep2;
  k = 0;
  for j = 1:3
      dif = abs(lat6(j)-lat104);
      tmp = sprintf('s%d',j);
      xx = eval(tmp);
      [m,n] = size(xx);
      xc = repmat(lat6(j),m,1);
      yc = -1*xx(:,1)*rr(j);
      qut = xx(:,9);
      field_do1(k+1:k+m) = qut;
      id = find(dif == min(dif));
      mod1 = interp1(dep22(:,id),ph1(:,id),yc,'linear','extrap');
      mod_do1(k+1:k+m) = mod1;
      k = k+m;  
  end
   k = 0;
  for j = 4:11
      dif = abs(lat6(j)-lat104);
      tmp = sprintf('s%d',j);
      xx = eval(tmp);
      [m,n] = size(xx);
      xc = repmat(lat6(j),m,1);
      yc = -1*xx(:,1)*rr(j);
      qut = xx(:,9);
      field_do2(k+1:k+m) = qut;
      id = find(dif == min(dif));
      mod1 = interp1(dep22(:,id),ph1(:,id),yc,'linear','extrap');
      mod_do2(k+1:k+m) = mod1;
      k = k+m;  
  end 
  k = 0;
  for j = 11:14
      dif = abs(lat6(j)-lat104);
      tmp = sprintf('s%d',j);
      xx = eval(tmp);
      [m,n] = size(xx);
      xc = repmat(lat6(j),m,1);
      yc = -1*xx(:,1)*rr(j);
      qut = xx(:,9);
      field_do3(k+1:k+m) = qut;
      id = find(dif == min(dif));
      mod1 = interp1(dep22(:,id),ph1(:,id),yc,'linear','extrap');
      mod_do3(k+1:k+m) = mod1;
      k = k+m;  
  end  
  field_do = [field_do1,field_do2,field_do3];
  mod_do = [mod_do1,mod_do2,mod_do3];
%   figure(1);set(gcf,'position',[1153 461 700 519])
%   plot(field_do,mod_do,'bo');
%   a = polyfit(field_do,mod_do,1);
% %   xx = mean(field_do./mod_do);
%   hh = refline(a(1),a(2)); set(hh,'color','r')
%   set(gca,'xlim',[6.5 8.5],'ylim',[6.5 8.5]);
%   xlabel('Field Measurement','fontweight','bold','fontsize',14)
%   ylabel('Model Simulation','fontweight','bold','fontsize',14)
%   title(['pH June 2016'],'fontweight','bold','fontsize',14);
%   r1 = refline(1,0); 
%   set(r1,'linestyle','--','color','k');
%   rmse = sqrt(sum((field_do-mod_do).^2)/length(field_do));  
%   rmse = floor(rmse*100)/100;
%   text(0.1,0.9,['RMSE = ',num2str(rmse)],'fontweight','bold','Units','normalized','fontsize',14);
%   
  %%
  load august2016
  dep_int = interp1(lat104,dep104,lat6,'linear','extrap');
  rr = dep_int./dep6;
  ph1 = squeeze(mean(name1(:,:,9),3));
  dep22 = -1*Sz.*dep2;
  k = 0;
  for j = 1:3
      dif = abs(lat6(j)-lat104);
      tmp = sprintf('s%d',j);
      xx = eval(tmp);
      [m,n] = size(xx);
      xc = repmat(lat6(j),m,1);
      yc = -1*xx(:,1)*rr(j);
      qut = xx(:,9);
      field_do11(k+1:k+m) = qut;
      id = find(dif == min(dif));
      mod1 = interp1(dep22(:,id),ph1(:,id),yc,'linear','extrap');
      mod_do11(k+1:k+m) = mod1;
      k = k+m;  
  end
    k = 0;
    for j = 4:11
      dif = abs(lat6(j)-lat104);
      tmp = sprintf('s%d',j);
      xx = eval(tmp);
      [m,n] = size(xx);
      xc = repmat(lat6(j),m,1);
      yc = -1*xx(:,1)*rr(j);
      qut = xx(:,9);
      field_do22(k+1:k+m) = qut;
      id = find(dif == min(dif));
      mod1 = interp1(dep22(:,id),ph1(:,id),yc,'linear','extrap');
      mod_do22(k+1:k+m) = mod1;
      k = k+m;  
    end
      k = 0;
    for j = 12:16
      dif = abs(lat6(j)-lat104);
      tmp = sprintf('s%d',j);
      xx = eval(tmp);
      [m,n] = size(xx);
      xc = repmat(lat6(j),m,1);
      yc = -1*xx(:,1)*rr(j);
      qut = xx(:,9);
      field_do33(k+1:k+m) = qut;
      id = find(dif == min(dif));
      mod1 = interp1(dep22(:,id),ph1(:,id),yc,'linear','extrap');
      mod_do33(k+1:k+m) = mod1;
      k = k+m;  
    end
    
  field_do8 = [field_do,field_do11,field_do22,field_do33];
  mod_do8 = [mod_do,mod_do11,mod_do22,mod_do33];
  figure(2);set(gcf,'position',[478   210   961   712])
  h1 = plot(field_do1,mod_do1,'ks','markerfacecolor','r');
  hold on
  h2 = plot(field_do2,mod_do2,'ks','markerfacecolor','b');
  h3 = plot(field_do3,mod_do3,'ks','markerfacecolor','g');
  h4 = plot(field_do11,mod_do11,'ko','markerfacecolor','r');
  h5 = plot(field_do22,mod_do22,'ko','markerfacecolor','b');
  h6 = plot(field_do33,mod_do33,'ko','markerfacecolor','g');
  a = polyfit(field_do8,mod_do8,1);
  hh = refline(a(1),a(2)); set(hh,'color','r','linewidth',2,'linestyle','--')
  set(gca,'xlim',[6.5 8.5],'ylim',[6.5 8.5],'fontweight','bold','fontsize',14);
  xlabel('Field Measurement','fontweight','bold','fontsize',14)
  ylabel('Model Simulation','fontweight','bold','fontsize',14)
  title(['pH Summer 2016'],'fontweight','bold','fontsize',14);
  r1 = refline(1,0); 
  set(r1,'linestyle','--','color','k','linewidth',3);
  rmse = sqrt(sum((field_do-mod_do).^2)/length(field_do));  
  rmse = floor(rmse*100)/100;
  text(0.1,0.9,['RMSE = ',num2str(rmse)],'fontweight','bold','Units','normalized','fontsize',14);
  l1 = legend([h1,h2,h3],'Upper Bay','Mid bay','Lower bay','location','southeast');
  plot(7.6262,8.3660,'ko','markersize',7)
  plot(7.9074,8.366,'ks','markersize',7)
  text(7.6662,8.3760,'June','fontweight','bold','fontsize',14)
  text(7.9462,8.3760,'August','fontweight','bold','fontsize',14)
  
  