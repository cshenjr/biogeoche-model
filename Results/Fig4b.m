 % calculate ph and plot
  close all; clear; clc
  load vertical_model_result_p
  Sz(1,:) = 0.0;
  
  load month_data  
  dep104 = dep2(1,:);   % model depth for 104 vertical profiles
  lon = reshape(lon,104,20);lon = lon';  % latitude matrix
  lat104 = lon(1,:); 
  
  month = {'mar','apr','may','jun','jul','aug','sep','oct','nov','dec'};
  n = length(month);
  
  load cp
  
  k = 0;
  for i = 1:n
     eval(['mon = ', month{i},';']);
      
     dep_int = interp1(lat104,dep104,mon(:,2),'linear','extrap');
     rr = dep_int./mon(:,9);  % depth ratio between model and field 

     %%  image comparison  
     gap = 12+i:14+i;
     h_dic = squeeze(mean(ta(:,:,gap),3));

     yc = -1*mon(:,5).*rr;
     m = length(yc);
     sze = repmat(100,m,1);

     dep22 = -1*Sz.*dep2;
     mod_int = zeros(m,1);
     for j = 1:m
       dif = abs(mon(j,2)-lat104);  
       id = find(dif == min(dif));  
       mod_int(j) = interp1(dep22(:,id),h_dic(:,id),yc(j),'linear','extrap');
     end
     
     tot(k+1:k+m,1) = mon(:,8);     %field
     tot(k+1:k+m,2) = mod_int; % model
     tot(k+1:k+m,3) = mon(:,2); % latitude
     tot(k+1:k+m,4) = mon(:,6); % salinit    
     k = k+m;
%      pause;
  end

  [a b] = sort(tot(:,3),1,'descend');
  tot(:,1) = tot(b,1);
  tot(:,2) = tot(b,2);
  tot(:,3) = tot(b,3);
  tot(:,4) = tot(b,4); 
  
  station_lat = [39.4415,39.3487, 39.1637, 38.9621, 38.8259,...
       38.5551, 38.3187, 38.1371, 37.9101, 37.4868, 37.2365];
  sid = {'CB21','CB22','CB32','CB858','CB41','CB43',...
      'CB51','CB52','CB53','CB62','CB64'};
  
  sid1 = {'mCB21','mCB22','mCB32','mCB858','mCB41','mCB43',...
      'mCB51','mCB52','mCB53','mCB62','mCB64'};
  


  for i = 1:length(sid)
      lat1 = station_lat(i);
      lat1 = floor(lat1*1000);
      lat2 = floor(tot(:,3)*1000);
      id = find(lat2 == lat1);
      eval([sid{i},'= tot(id,1);']);
      eval([sid1{i},'= tot(id,2);']);
  end
  figure(1);set(gcf,'position',[402         392        1257         549])
  tt = [CB21;mCB21;CB22;mCB22;CB32;mCB32;CB858;mCB858;CB41;mCB41;CB43;mCB43;...
      CB51;mCB51;CB52;mCB52;CB53;mCB53;CB62;mCB62;CB64;mCB64];
  
  tcluster = nan(22,77);
  for i = 1:11
      eval(['temp = ',sid{i},';']);
      eval(['temp1 = ',sid1{i},';']);
      nn = length(temp);
      tcluster(2*i-1,1:nn) = temp;
      tcluster(2*i,1:nn) = temp1;
  end
  
  x1 = linspace(1,18,11);
  x2 = x1+0.4;
  gp1 = zeros(1,22);
  gp1(1:2:end) = x1;
  gp1(2:2:end) = x2;
  

  h1 = boxplot(tcluster','positions',gp1,'symbol','','width',0.15,'color','k');
  set(h1,'linewidth',1)
  patch(get(h1(50),'XData'),get(h1(50),'YData'),'k','FaceAlpha',.3);

  set(gca,'xtick',x1);
  set(gca,'xticklabel',  {'CB21','CB22','CB32','CB858','CB41','CB43',...
      'CB51','CB52','CB53','CB62','CB64'},'fontweight','bold','fontsize',12)
  xlabel('Stations','fontweight','bold','fontsize',14)
  ylabel('TA (µ mol/kg)','fontweight','bold','fontsize',14)
 
  for i = 12:14:165
  patch(get(h1(i),'XData'),get(h1(i),'YData'),'k','FaceAlpha',.3);
  end
  
set(gca,'ylim',[600 2200])
set(gca,'xlim',[0 19])
hold on;
hh = plot(15.2,1019.6,'ks','markersize',12,'linewidth',1.0);
hh1 = plot(15.2,910.6,'ks','markersize',12,'linewidth',1.0,'markerfacecolor',[0.7 0.7 0.7]);

text(15.5,1019.6,'Field Observation','fontweight','bold','fontsize',12)
text(15.5,910.6,'Model Simulation','fontweight','bold','fontsize',12)



  
  
