close all; clear; clc;

load airsea

flux2 = max(flux2,-200);
flux2 = -flux2;

t = t-t(1);
t = t';
% plot(t)
up = find(lon>=39.0);
mid = find(lon<39.0 & lon>37.4);
low = find(lon <=37.4);

fu = mean(flux2(:,up),2);
fm = mean(flux2(:,mid),2);
fl = mean(flux2(:,low),2);

sp = find(t>=60 & t<= 150);
sm = find(t>=150 & t<= 240);
fa = find(t>=240 & t<= 330);
wi = find(t<60 | t>330);

x1 = fu(sp); 
x2 = fu(sm);
x3 = fu(fa);
x4 = fu(wi);

e1 = std(x1);
e2 = std(x2);
e3 = std(x3);
e4 = std(x4);

aa = 12;
figure(1); set(gcf,'position',[372         520        1117         399])
errorbar(1,mean(x1),e1,'ko','markersize',aa,'markerfacecolor','k','linewidth',1);
hold on
errorbar(2,mean(x2)-40,e2,'ro','markersize',aa,'markerfacecolor','r','linewidth',1);
errorbar(3,mean(x3),e3,'bo','markersize',aa,'markerfacecolor','b','linewidth',1);
errorbar(4,mean(x4),e4,'mo','markersize',aa,'markerfacecolor','m','linewidth',1);

mean(x1)+mean(x2)+mean(x3)-40 + mean(x4)

%%
x1 = fm(sp); 
x2 = fm(sm);
x3 = fm(fa);
x4 = fm(wi);

e1 = std(x1);
e2 = std(x2);
e3 = std(x3);
e4 = std(x4);

errorbar(5,mean(x1),e1,'ko','markersize',aa,'markerfacecolor','k','linewidth',1);
errorbar(6,mean(x2),e2,'ro','markersize',aa,'markerfacecolor','r','linewidth',1);
errorbar(7,mean(x3),e3,'bo','markersize',aa,'markerfacecolor','b','linewidth',1);
errorbar(8,mean(x4),e4,'mo','markersize',aa,'markerfacecolor','m','linewidth',1);

mean(x1)+mean(x2)+mean(x3)+ mean(x4)

%%
x1 = fl(sp); 
x2 = fl(sm);
x3 = fl(fa);
x4 = fl(wi);

e1 = std(x1);
e2 = std(x2);
e3 = std(x3);
e4 = std(x4);

h11 = errorbar(9,mean(x1),e1,'ko','markersize',aa,'markerfacecolor','k','linewidth',1);
h22 =errorbar(10,mean(x2),e2,'ro','markersize',aa,'markerfacecolor','r','linewidth',1);
h33 =errorbar(11,mean(x3),e3,'bo','markersize',aa,'markerfacecolor','b','linewidth',1);
h44 =errorbar(12,mean(x4),e4,'mo','markersize',aa,'markerfacecolor','m','linewidth',1);
mean(x1)+mean(x2)+mean(x3)+mean(x4)

x1 = [1.2 5.2];
x2 = [2.2 6.2 10.2];
x3 = [3.2 7.2 11.2];
y1 = [20.95 -4.3674];
y2 = [16.94 0.64 1.02];
y3 = [33.94 2.53 6.48];

h111 = plot(x1(1),y1(1),'kp','markersize',10,'markerfacecolor','k')
plot(x1(2),y1(2),'kp','markersize',10,'markerfacecolor','k')
h222 = plot(x2(1),y2(1),'rp','markersize',10,'markerfacecolor','r')
plot(x2(2),y2(2),'rp','markersize',10,'markerfacecolor','r')
plot(x2(3),y2(3),'rp','markersize',10,'markerfacecolor','r')
h333 = plot(x3(1),y3(1),'bp','markersize',10,'markerfacecolor','b')
plot(x3(2),y3(2),'bp','markersize',10,'markerfacecolor','b')
plot(x3(3),y3(3),'bp','markersize',10,'markerfacecolor','b')

l1 = legend([h11 h22 h33 h44],'Spring',...
    'Summer','Fall','Winter','location','northeast','orientation','horizontal')
set(l1,'box','off')

plot(7, 41.6923, 'kp','markersize',10,'markerfacecolor','k');
text(7.35, 41.6923, 'May','fontweight','bold','fontsize',12)
plot(8.6, 41.6923, 'rp','markersize',10,'markerfacecolor','r');
text(8.96, 41.6923, 'June-August','fontweight','bold','fontsize',12)
plot(10.9, 41.6923, 'bp','markersize',10,'markerfacecolor','k');
text(11.26, 41.6923, 'October','fontweight','bold','fontsize',12)

plot([6.67 12.84 12.84 6.67 6.67],[55.53 55.53 36.76 36.76 55.53],...
    '-k')


h = refline(0,0);
set(h,'linestyle','--','color','k')
plot([4.5 4.5], [-40 60],'k--')
plot([8.5 8.5], [-40 60],'k--')
set(gca,'xlim',[0 13],'xtick',(1:12),'fontweight','bold','fontsize',14);
set(gca,'xticklabel',{'','Upper Bay','','','','                Mid Bay','',''...
    ,'','              Lower Bay','',''});
set(gca,'ylim',[-40 60])
ylabel('Air-sea flux (mmol C m^{-2} d^{-1})')


