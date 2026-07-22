%% 上升流指数
clear;clc;close all;

load CCMP_wind1998_2024.mat;    %加载wind的mat数据
%% 物理参数
[lo,la]=meshgrid(Lon,Lat);  %组成网格
rhoa_sea=1.025*10^3;   %海水密度
rhoa_air=1.293;        %用标准空气密度   
Ommiga=7.292*10^(-5);   %地球自转角速度
f=2*Ommiga*sind(la');   %科氏参数
Re = 6371000;    %地球平均半径6371km
%% 计算公式和结果
U = hypot(u,v);   %u v分量风合成为总风速
Cd=(0.8+0.065*U).*10^-3;  %drag coefficient

taou=rhoa_air.*Cd.*U.*u;   %纬向（东向）风应力 
taov=rhoa_air.*Cd.*U.*v;   %经向（北向）风应力 
tao = hypot(taou, taov);   %合成总的风应力 

VN_angle=10;  %夏季风越南风向岸线角为10° 最强上升流
LZ_angle=45;  %冬季风吕宋风向岸线角为45° 中强度上升流
UI_VN = (tao./(f*rhoa_sea)).*cosd(VN_angle);  %上升流指数 越南 夏季6 7 8月  m^2/s
UI_LZ = (tao./(f*rhoa_sea)).*cosd(LZ_angle);  %上升流指数 吕宋岛 冬季12 1 2月  m^2/s

%注意：UI_LZ这个都是用cos45°计算的 UI_VN使用cos10°计算的 没有区分时间，后续需要划定区域来画图
%% 风应力旋度
[dx_deg,~] = gradient(lo);   %gradient第一个输出 dx： 列方向（东西 → 经度方向）的差分 第二个输出 dy：行方向（南北 → 纬度方向）的差分
[~,dy_deg] = gradient(la);  %dx_deg：每个网格点东西方向的经度间隔（角度） dy_deg：每个网格点南北方向的纬度间隔（角度）
dx_dis = dx_deg.*(cosd(la)*Re)*2*pi/360;  
dy_dis = dy_deg.*Re*2*pi/360;  %经纬度转距离的核心步骤

curl=nan(size(tao,1),size(tao,2),size(tao,3));  %初始化curl矩阵
for n=1:size(taou,3)
tu=taou(:,:,n);  %取出第n层的二维风应力场（空间分布）
tv=taov(:,:,n);
[~,dtaux_dy] = gradient(tu');  %东向风应力沿南北方向的偏导数
[dtauy_dx,~] = gradient(tv');  %北向风应力沿东西方向的偏导数
curlz_latlon = dtauy_dx./dx_dis-dtaux_dy./dy_dis;   %风应力旋度核心公式
curl(:,:,n) = curlz_latlon';    %风应力旋度 N/m^3 （10^-6数量级） 正值：气旋式风应力旋度（海洋中通常引发上升流）负值：反气旋式风应力旋度（海洋中通常引发下沉流）
end

save('UI&curl1998-2024.mat','curl',"UI_LZ","UI_VN","Lon","Lat");



%% 画图看看试验效果  colorbar范围已经确定好了

Ui1=nanmean(UI_VN(:,:,213:243),3); %1998.8 越南 夏季
figure('Position',[200,200,1000,900],'Color','W','Resize','off');
m_proj('equidistant','lon',[105 123],'lat',[5  25]);  
m_contourf(Lon,Lat,Ui1',500,'linestyle','none');
hold on
m_gshhs_i('linewidth',1,'color','k');
m_gshhs_i('patch',[0.6,0.6,0.6]);
m_grid('linestyle','none','box','on','linewidth',1.1,'fontname','times','fontsize',20,'xticklabel', []);
hold on
title('Aug UI(m^{2}s^{-1})');
colormap(nclCM(382,28));
caxis([0   3]); 
cb=colorbar('southoutside');
set(cb,'linewidth',1.2,'fontsize',18,'edgecolor','k');
ylabel(cb,'UI(m^{2} s^{-1})','FontName','times','FontSize',18)


%
Ui2=nanmean(UI_LZ(:,:,335:365),3); %1998.12  吕宋岛 冬季比夏季ui高
figure('Position',[200,200,1000,900],'Color','W','Resize','off');
m_proj('equidistant','lon',[105 123],'lat',[5  25]);   
m_contourf(Lon,Lat,Ui2',500,'linestyle','none');
hold on
m_gshhs_i('linewidth',1,'color','k');
m_gshhs_i('patch',[0.6,0.6,0.6]);
m_grid('linestyle','none','box','on','linewidth',1.1,'fontname','times','fontsize',20,'xticklabel', []);
hold on
title('Dec UI(m^{2}s^{-1})');
colormap(nclCM(382,28));
caxis([1   7]);  
cb=colorbar('southoutside');
set(cb,'linewidth',1.2,'fontsize',18,'edgecolor','k');
ylabel(cb,'UI(m^{2}s^{-1})','FontName','times','FontSize',18)

%%
curl1=nanmean(curl(:,:,213:243),3); %1998.12  吕宋岛 curl

curl_plot1 = curl1 ./ 1e-6;


figure('Position',[200,200,1000,900],'Color','W','Resize','off');
m_proj('equidistant','lon',[105 123],'lat',[5  25]);  
m_contourf(Lon,Lat,curl_plot1',500,'linestyle','none');
shading interp
hold on
m_gshhs_i('linewidth',1,'color','k');
m_gshhs_i('patch',[0.6,0.6,0.6]);
m_grid('linestyle','none','box','on','linewidth',1.1,'fontname','times','fontsize',20,'xticklabel', []);
hold on
colormap(nclCM(382,28));
caxis([-1   1]);  
title('Aug wind stress curl (10^{-6} N m^{-3})');
cb = colorbar;
set(cb,'linewidth',1.2,'fontsize',18,'edgecolor','k');
ylabel(cb,'10^{-6} N m^{-3}');


%
curl2=nanmean(curl(:,:,335:365),3); %1998.12  吕宋岛 curl
curl_plot2 = curl2 ./ 1e-6;
figure('Position',[200,200,1000,900],'Color','W','Resize','off');
m_proj('equidistant','lon',[105 123],'lat',[5  25]);  
m_contourf(Lon,Lat,curl_plot2',500,'linestyle','none');
shading interp
hold on
m_gshhs_i('linewidth',1,'color','k');
m_gshhs_i('patch',[0.6,0.6,0.6]);
m_grid('linestyle','none','box','on','linewidth',1.1,'fontname','times','fontsize',20,'xticklabel', []);
hold on
colormap(nclCM(382,28));
caxis([-2.5   2.5]);  
title('Dec wind stress curl (10^{-6} N m^{-3})');
cb = colorbar;
set(cb,'linewidth',1.2,'fontsize',18,'edgecolor','k');
ylabel(cb,'10^{-6} N m^{-3}');