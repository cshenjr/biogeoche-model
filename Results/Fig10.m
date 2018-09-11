close all; clc; clear;

load diag_model_ta
load diag_model

% mid = find(lon>38.7 & lon<39.3);
% mid = (10:20);
mid = find(lon>38.7 & lon<39.3);
bot = 10:20; 
n = length(bot);

tt = tt-12051;
aa = 1;
% figure(1); set(gcf,'position',[461 561 1110  361]);
%% part bay
resp_2 = squeeze(mean(mean(resp(bot,:,mid),1),3));
oxid_2 = squeeze(mean(mean(oxid(bot,:,mid),1),3));
cal_2 = squeeze(mean(mean(cal(bot,:,mid),1),3));
sod_2 = squeeze(mean(sod(:,mid),2))/n;
jhs_2 = squeeze(mean(jhs(:,mid),2))/n;
sed_2 = squeeze(mean(sed(:,mid),2))/n;
suf_2 = squeeze(mean(mean(sulf(bot,:,mid),1),3));

% time = (150:240);
ab = 20;
for i = 1:12
  time = ((i-1)*30+1):(i*30);  
  ar1(i) = mean(resp_2(time)+ oxid_2(time))*ab;
  cd1(i) = mean(cal_2(time))*ab;
  sr1(i) = mean(suf_2(time))*ab;
  sod1(i) = mean(jhs_2(time)+sod_2(time))*ab;
  
  ar2(i) = 0;
  cd2(i) = 2*cd1(i);
  sr2(i) = 1.14*sr1(i);
  sod2(i) = mean(jhs_2(time)+sed_2(time))*ab;
end

load diag_model_carb 
dic_2 = squeeze(mean(mean(dic(bot,:,mid),1),3));
ta_2 = squeeze(mean(mean(ta(bot,:,mid),1),3));
po4_2 = squeeze(mean(mean(po4(bot,:,mid),1),3));
sal_2 = squeeze(mean(mean(sal(bot,:,mid),1),3));
temp_2 = squeeze(mean(mean(temp(bot,:,mid),1),3));

for i = 1:12
  time = 30*i;
  dic1(i) = dic_2(time);
  ta1(i) = ta_2(time);
  po41(i) = po4_2(time);
  sal1(i) = sal_2(time);
  temp1(i) = temp_2(time);
end

ph = pH_cal(temp1,sal1,dic1,ta1,po41,zeros(size(sal1)));

dic1 = dic1-ar1; ta1 = ta1-ar2;
ph_ar = pH_cal(temp1,sal1,dic1,ta1,po41,zeros(size(sal1)));

dic1 = dic1-sr1; ta1 = ta1-sr2;
ph_sr = pH_cal(temp1,sal1,dic1,ta1,po41,zeros(size(sal1)));

dic1 = dic1-cd1; ta1 = ta1-cd2;
ph_cd = pH_cal(temp1,sal1,dic1,ta1,po41,zeros(size(sal1)));

dic1 = dic1-sod1; ta1 = ta1-sod2;
ph_sod = pH_cal(temp1,sal1,dic1,ta1,po41,zeros(size(sal1)));

ef_ar = ph-ph_ar;
ef_sr = ph_ar-ph_sr;
ef_cd = ph_sr-ph_cd;
ef_sod = ph_cd-ph_sod;


figure(1);set(gcf,'position',[437   574   975   350])
plot(1:12,ef_ar,'r-o','linewidth',1.5);
hold on;
plot(1:12,ef_sod,'b-^','linewidth',1.5);
plot(1:12,ef_sr,'m-s','linewidth',1.5);
plot(1:12,ef_cd,'k-d','linewidth',1.5);

legend('AR','SD','SR','CD','location','southeast');
set(gca,'xlim',[0 13],'ylim',[-1.5 0.5],'fontweight','bold','fontsize',14);
set(gca,'xtick',0.5:11.5);set(gca,'xticklabel',{'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';...
    'Sep';'Oct';'Nov';'Dec'})
xlabel('Month of 2016');ylabel('\DeltapH')












