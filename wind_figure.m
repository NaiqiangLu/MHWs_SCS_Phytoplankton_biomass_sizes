clc;clear;close all;
load CCMP_wind1998_2024.mat

% ================== 1. 取30天平均 ==================
u1  = nanmean(u(:,:,1:30),3);
v1  = nanmean(v(:,:,1:30),3);
ws1 = nanmean(ws(:,:,1:30),3);


% %% ================== 2. topo掩膜 去掉风向箭头（原方法） ==================
%
% % 读取海洋深度
% lo = ncread('topo_20.1.nc','lon');
% la = ncread('topo_20.1.nc','lat');
% depth = ncread('topo_20.1.nc','z');
% 
% % 裁剪区域
% x2 = find(lo>=100 & lo<=125);
% y2 = find(la>=0 & la<=25);
% 
% la = la(y2);
% lo = lo(x2);
% dep = depth(x2,y2);
% 
% % 构建索引网格
% x = 1:length(la);
% y = 1:length(lo);
% [X,Y] = meshgrid(x,y);
% 
% % 插值到风场网格大小
% xi = 1:(length(la)/length(lat)):length(la);
% yi = 1:(length(lo)/length(lon)):length(lo);
% [Xi,Yi] = meshgrid(xi,yi);
% 
% Z = interp2(X,Y,dep,Xi,Yi,'spline');
% 
% % ===== 构建掩膜 =====
% % 海洋=1，陆地=NaN
% mask = ones(size(Z));
% mask(Z > 0) = NaN;   % 陆地
% mask(Z <= 0) = 1;    % 海洋
% 
% % ===== 应用掩膜（关键）=====
% u1 = u1 .* mask;
% v1 = v1 .* mask;





%% ================== 2. 绘图 ==================
% ================== 3. 风矢量 ==================
[Lo,La] = meshgrid(Lon,Lat);
figure;
set(gcf,'position',[200,200,900,800],'Color','W');

m_proj('equidistant','lon',[100 125],'lat',[0 25]);

% --- 风速填色 ---
m_contourf(Lon,Lat,ws1',100,'linestyle','none'); hold on;


hold on
step = 5;
% 关键：统一缩放系数（自己控制）
scale_factor = 0.25;   % 可调（非常关键）

m_quiver(Lo(1:step:end,1:step:end), ...
         La(1:step:end,1:step:end), ...
         u1(1:step:end,1:step:end)'*scale_factor, ...
         v1(1:step:end,1:step:end)'*scale_factor, ...
         0, ...   % 关闭自动缩放
         'color','k', ...
         'linewidth',1);

hold on
% --- 海岸线 ---
m_gshhs_i('linewidth',1,'color','k');
m_gshhs_i('patch',[0.5 0.5 0.5]);
% --- 网格 ---
m_grid('linestyle','none','box','on','linewidth',1.2,'fontsize',20);

% ================== 4. 色标 ==================
colormap(m_colmap('diverging',18));
cb = colorbar;

caxis([0 12]);

set(cb,'tickdir','in');
set(cb,'linewidth',1,'fontsize',20,'edgecolor','k');
set(cb,'YTick',0:2:12);

ylabel(cb,'Wind Speed (m/s)','FontSize',22,'FontWeight','bold');

%% ================== 5. 真实比例参考箭头（居中版） ==================

Uref = 6;   % 参考风速

% ===== 自动右下角 =====
lon_ref = Lon(end) - 3;
lat_ref = Lat(1) + 2;

% ===== 图例框尺寸 =====
box_w = 4;    % 经度宽度
box_h = 2.4;  % 纬度高度

% ===== 白色背景框（中心 = lon_ref, lat_ref）=====
lon_box = [lon_ref-box_w/2, lon_ref+box_w/2, lon_ref+box_w/2, lon_ref-box_w/2];
lat_box = [lat_ref-box_h/2, lat_ref-box_h/2, lat_ref+box_h/2, lat_ref+box_h/2];

m_patch(lon_box, lat_box, 'w', ...
    'edgecolor','none','linewidth',1);   % 加边框更清晰

% ===== 箭头（居中 & 水平）=====
arrow_len = Uref * scale_factor;   % 长度必须一致

x_start = lon_ref - arrow_len/2;
y_arrow = lat_ref + 0.4;   % 略偏上

m_quiver(x_start, y_arrow, arrow_len, 0, ...
         0, ...
         'color','k','linewidth',1.1,'MaxHeadSize',1.1);

% ===== 文字（居中）=====
m_text(lon_ref, lat_ref - 0.5, ...
    sprintf('%d m/s',Uref), ...
    'fontsize',16, ...
    'fontweight','bold', ...
    'horizontalalignment','center');


% ================== 5. 保存 ==================
print(gcf,'CCMP_day1-30_wind','-dtiff','-r100');
