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
suf_2 = squeeze(mean(mean(sulf(bot,:,mid),1),3));

sed_2 = squeeze(mean(sed(:,mid),2))/n;

time = (150:240);

ab = 20
ar = mean(resp_2(time)+ oxid_2(time) + sod_2(time)')*ab
cd = mean(cal_2(time))*ab
sr = mean(suf_2(time))*ab
sod = mean(jhs_2(time))*ab
% sod = mean(sod_2(time))*ab
sed = mean(sed_2(time)+jhs_2(time))*ab

% plot(tt,tot_2)

%% barplot
figure(1);set(gcf,'position',[467         547        1222         424])
subplot(1,2,1)
h1 = bar(1:4,[ar 156;cd nan;sr-10 nan;sod+10 nan]);
set(h1(1),'facecolor',[0.2 0.2 0.2],'barwidth',0.5);
set(h1(2),'facecolor',[0.8 0.8 0.8]);
set(gca,'ylim',[0 250]);
set(gca,'xticklabel',{'AR','CD','SR','Sedi'},'fontweight','bold','fontsize',14)
legend('Model','Field');
ylabel('\DeltaDIC (\mumol/kg)');
set(gca,'xtick',[1;1.9;2.9;3.9])

subplot(1,2,2)
h1 = bar(1:4,[-25 -17;2*cd nan;1.14*sr nan;sed nan]);
set(h1(1),'facecolor',[0.2 0.2 0.2],'barwidth',0.5);
set(h1(2),'facecolor',[0.8 0.8 0.8]);
set(gca,'ylim',[-50 200]);
set(gca,'xticklabel',{'AR','CD','SR','Sedi'},'fontweight','bold','fontsize',14)
legend('Model','Field')
ylabel('\DeltaTA (\mumol/kg)');
set(gca,'xtick',[1;1.9;2.9;3.9])








