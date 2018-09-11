close all; clc; clear;

load ph_base

do = ph;

s1 = 7.5;
hh = h1(id);
pm = pm1(id);
pn = pn1(id);

xx = squeeze(do(100,1,:));
id = find(xx~=0);
n = length(id);
vt1 = zeros(n,20);
vt4 = zeros(n,20);
for i = 1:n-1
    do1 = squeeze(do(:,:,i));
    for j = 1:20
        do2 = squeeze(do1(:,j));
        a1 = (do2<=s1);
        vol1 = a1.*hh(:)/20./pm(:)./pn(:);
        vt1(i,j) = sum(vol1);
    end
end

vt1 = sum(vt1,2);

figure(1);set(gcf,'position',[680         617        1048         361])
plot(tt-tt(1),2*vt1/1e9,'r-','linewidth',2);
hold on;
xlabel('Year 2016','fontweight','bold','fontsize',14);
ylabel('Acid water volume (km^3)','fontweight','bold','fontsize',14);

%%
load ph_more
do = ph;
hh = h1(id);
pm = pm1(id);
pn = pn1(id);

xx = squeeze(do(100,1,:));
id = find(xx~=0);
n = length(id);
vt1 = zeros(n,20);
for i = 1:n-1
    do1 = squeeze(do(:,:,i));
    for j = 1:20
        do2 = squeeze(do1(:,j));
        a1 = (do2<=s1);
        vol1 = a1.*hh(:)/20./pm(:)./pn(:);
        vt1(i,j) = sum(vol1);
    end
end

vt1 = sum(vt1,2);

plot(tt-tt(1),2*vt1/1e9,'b-','linewidth',2);
hold on;

%%
load ph_less
do = ph;
hh = h1(id);
pm = pm1(id);
pn = pn1(id);

xx = squeeze(do(100,1,:));
id = find(xx~=0);
n = length(id);
vt1 = zeros(n,20);
for i = 1:n-1
    do1 = squeeze(do(:,:,i));
    for j = 1:20
        do2 = squeeze(do1(:,j));
        a1 = (do2<=s1); 
        vol1 = a1.*hh(:)/20./pm(:)./pn(:);
        vt1(i,j) = sum(vol1);
    end
end

vt1 = sum(vt1,2);

plot(tt-tt(1),2*vt1/1e9,'k','linewidth',2);

set(gca,'xlim',[0 370])
set(gca,'xtick',[(15:30:365)],'fontsize',14);
set(gca,'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul',...
       'Aug','Sep','Oct','Nov','Dec'},'fontweight','bold')

legend('1.0X','1.5X','0.5X','fontweight','bold','location','northeast',...
    'orientation','horizontal')

%   text( 28.3251 ,  20.3179, '0.25X','fontweight','bold')
%   text( 28.3251  , 17.3176, 'Base','fontweight','bold')
%   text( 28.3251  , 14.3179, '2.0X','fontweight','bold')

%   plot([50.3251, 80.3251], [20.3179, 20.3179],'k-.','linewidth',2)
%   plot([50.3251, 80.3251], [17.3176, 17.3176],'k-','linewidth',2)
%   plot([50.3251, 80.3251], [14.3179, 14.3179],'k--','linewidth',2)
