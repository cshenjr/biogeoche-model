close all; clear; clc

load data1316
load station1316

var = {'CB21','CB22','CB32','CB41C','CB42C','CB43C',...
    'CB44','CB51','CB52','CB53','CB54','CB55','CB64'};

figure(1);
set(gcf,'position',[532         165        1162         771]); 
h1 = subplot(311);
do_bot = squeeze(do(end,:,6));
tt1 = tt-tt(1);
plot(tt1,do_bot,'k-','linewidth',2);
hold on;
do1_bot = CB43C.bot.DO;
time = CB43C.bot.t_DO;
plot(time,do1_bot,'ko','markerfacecolor','r');
set(gca,'ylim',[0 14],'fontweight','bold','fontsize',14)
set(gca,'xtick',[]);
legend('Model','Field','location','northeast','orientation','horizontal')

%%%%%%%rmse
ido = interp1(tt1,do_bot,time,'linear','extrap');
rmsedo = sqrt(mean((ido-do1_bot).^2));
disp(rmsedo)

plot(370*ones(20,1),linspace(0,20,20),'k--')
plot(740*ones(20,1),linspace(0,20,20),'k--')
plot(1110*ones(20,1),linspace(0,20,20),'k--')
ylabel('DO (mg/L)','fontsize',14)
text(78.3582,12.6084,'\itCB4.3C','fontweight','bold','fontsize',14);
hh1 = refline(0,2);set(hh1,'linewidth',1,'linestyle','--');
hold off;
%%
load ph_cal1_ca1
[pHx] = pH_cal(temp,sal,dic,ta,po4,sit);
phx = real(pHx);
dif = abs(lon-CB43C.lat);
id = find(dif == min(dif));
ph_station1 = squeeze(phx(:,:,id));

load ph_cal1_ca1_noca
[pHx] = pH_cal(temp,sal,dic,ta,po4,sit);
phx = real(pHx);
ph_station2 = squeeze(phx(:,:,id));

ph_bot1 = squeeze(ph_station1(end,:));
ph_bot2 = squeeze(ph_station2(end,:));
xx  = 0.8*ph_bot2+0.2*ph_bot1;
xx(360:380) = nan;
xx(720:760) = nan;

h2 = subplot(312);
tt1 = tt-tt(1);
p1 = plot(tt1,ph_bot1,'k-','linewidth',1.5);
hold on;
p2 = plot(tt1,xx,'--','linewidth',1.0,'color',[0.5 0.5 0.5]);

ph1_bot = CB43C.bot.PH;
time = CB43C.bot.t_PH;
plot(time,ph1_bot,'ko','markerfacecolor','r');
set(gca,'ylim',[7.0 8.2],'fontweight','bold','fontsize',14)
set(gca,'xtick',[]);

legend([p1 p2],'W/ CD','W/O CD','location','northeast','orientation','horizontal')

%%%%%
% iph = interp1(tt1,ph_bot,time,'linear','extrap');
% rmseph = sqrt(mean((iph-ph1_bot).^2));
% disp(rmseph)

plot(370*ones(20,1),linspace(0,20,20),'k--')
plot(740*ones(20,1),linspace(0,20,20),'k--')
plot(1110*ones(20,1),linspace(0,20,20),'k--')
ylabel('pH','fontsize',14)
hh1 = refline(0,7.5);set(hh1,'linewidth',1,'linestyle','--');
hold off
%%
load 2015
load 2016
load 2014
load 2013

time = [time13 time14 time15 time16];
vol = [vol13_2' vol14_2' vol15_2' vol16_2'];

id = find(vol>=25e9);
vol(id) = 0;
h3 = subplot(313);
plot(time-time(1),vol/1e9,'r-','linewidth',2);
load p2016
tt16 = tt; vol16 = vt1;
load p2015
tt15 = tt; vol15 = vt1;
load p2014
tt14 = tt; vol14 = vt1;
load p2013
tt13 = tt; vol13 = vt1;
time = [tt13 tt14 tt15 tt16];
vol = [vol13' vol14' vol15' vol16'];
hold on;
plot(time-time(1),vol/1e9,'k-','linewidth',2)  
set(gca,'ylim',[0 25]);
set(gca,'xtick',[0 370 740 1110 1480],'fontweight','bold','fontsize',14)
set(gca,'xticklabel',{'                                               2013',...
    '                                                       2014',...
    '                                                       2015',...
    '                                                       2016',''})
plot(370*ones(30,1),linspace(0,30,30),'k--')
plot(740*ones(30,1),linspace(0,30,30),'k--')
plot(1110*ones(30,1),linspace(0,30,30),'k--')
ylabel('Volume (km^3)','fontsize',14)
legend('DO <=2','pH<=7.5','location','northeast','orientation','horizontal')
hold off; 
set(h1,'position',[0.1300    0.6393    0.7491    0.2157]);
set(h2,'position',[0.1300    0.4096    0.7491    0.2157]);
set(h3,'position',[0.1300    0.1800    0.7491    0.2157]);   
   
   
   
   
   
   
   
   
   
   
   
   
   