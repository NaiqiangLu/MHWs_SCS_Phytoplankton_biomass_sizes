clc,clear;close all;
load Ibar1998_2024.mat

Ibarclim=mclim(Ibar,datenum(1998,1,1):datenum(2024,12,31),...
    datenum(1998,1,1),datenum(2024,12,31),31);  %数据:起止时间，气候基期，滑动窗口长度; 默认的窗口半宽为5


save('mhw_climway_Ibar.mat',"Ibarclim","Lon","Lat");



%% 查看图
clc; clear; close all;

load mhw_climway_Ibar.mat

%第一天气候态
Ibar = Ibarclim(:,:,1);

figure('Position',[100 100 900 800]);

m_proj('miller','lon',[105 125],'lat',[0 25]);

% ===== 1. 彩色填色 =====
m_contourf(Lon, Lat, Ibar, 100, 'linestyle','none'); 
shading interp;

hold on;
% [C,h] = m_contour(Lon, Lat, Ibar, [5 8 12],'k', 'LineWidth',1);  %三个等值线统一为黑色表示

% ===== 2. 三条等值线分别画 =====
% ---- 5（蓝色）----
[C1,h1] = m_contour(Lon, Lat, Ibar, [5 5], ...
    'LineWidth',1.8, 'LineColor','b');
clabel(C1,h1,'FontSize',15,'Color','k','LabelSpacing',500);

% ---- 8（绿色）----
[C2,h2] = m_contour(Lon, Lat, Ibar, [8 8], ...
    'LineWidth',1.8, 'LineColor','g');
clabel(C2,h2,'FontSize',15,'Color','k','LabelSpacing',500);

% ---- 12（红色）----
[C3,h3] = m_contour(Lon, Lat, Ibar, [12 12], ...
    'LineWidth',1.8, 'LineColor','y');
clabel(C3,h3,'FontSize',15,'Color','k','LabelSpacing',500);

% ===== 3. 地图要素 =====
m_gshhs_i('patch',[0.7 0.7 0.7]);
m_grid('linestyle','none','box','on','tickdir','in');

colormap(nclCM(156,20));
colorbar('location','southoutside');
caxis([0, 15]);

title('Mixed-layer averaged irradiance');
