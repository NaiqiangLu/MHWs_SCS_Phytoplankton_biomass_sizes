%% 插值降低1998-2024年9862每日叶绿素数据矩阵
% 基于经纬度的插值（推荐方法）
clear; clc;close all;

% ===== 目标网格（0.25°）=====
load sst_Lon_Lat.mat
x2 = find(Lon>=100 & Lon<=125); 
y2 = find(Lat>=0 & Lat<=25);
Lon = Lon(x2); 
Lat = Lat(y2);

% ===== 原始叶绿素网格 =====
load chla_lat_lon
x1 = find(lon>=100 & lon<=125); 
y1 = find(lat>=0 & lat<=25);
lon = lon(x1); 
lat = lat(y1);

% ===== 构建网格（关键）=====
[Lon_004, Lat_004] = meshgrid(lon, lat);   % 原始网格
[Lon_025, Lat_025] = meshgrid(Lon, Lat);   % 目标网格

chla_full = [];

for i = 1998:2024
    load(['chla_' num2str(i)]);
    
    chl = Chla(x1, y1, :);
    chl(chl < 0 ) = nan;   % 去异常值
    
    cchl = nan(length(Lat), length(Lon), size(chl,3));
    
    for j = 1:size(chl,3)
        % ===== 核心：经纬度插值 =====
        chll = interp2(Lon_004, Lat_004, chl(:,:,j)', ...
                       Lon_025, Lat_025, 'linear');
        
        % 注意：这里不再需要 flipud + 转置乱操作
        cchl(:,:,j) = chll;
    end
    
    % ===== 时间拼接 =====
    t_index = (datenum(i,1,1):datenum(i,12,31)) ...
              - datenum(1998,1,1) + 1;
    
    chla_full(:,:,t_index) = cchl;
end

save('Chla_19982024_interp025.mat', ...
     'Lon','Lat','chla_full','-v7.3');


%% 插值降低每日气候态366叶绿素分辨率
% 基于经纬度的气候态叶绿素插值（366天）
clear; clc;close all;

% ===== 目标网格（0.25°）=====
load sst_Lon_Lat
x2 = find(Lon>=100 & Lon<=125); 
y2 = find(Lat>=0 & Lat<=25);
Lon = Lon(x2); 
Lat = Lat(y2);

% ===== 原始叶绿素网格 =====
load chla_lat_lon
x1 = find(lon>=100 & lon<=125); 
y1 = find(lat>=0 & lat<=25);
lon = lon(x1); 
lat = lat(y1);

% ===== 气候态数据 =====
load climchl_1998_2024.mat   % climchl (lon × lat × 366 或 lat × lon × 366)

climchl0 = climchl(x1, y1, :);   % 先裁剪（避免重复计算）

% ===== 构建经纬度网格 =====
[Lon_004, Lat_004] = meshgrid(lon, lat);   % 原网格
[Lon_025, Lat_025] = meshgrid(Lon, Lat);   % 目标网格

% ===== 初始化 =====
climchla = nan(length(Lat), length(Lon), size(climchl0,3));

for j = 1:size(climchl0,3)
    
    % ⚠️ 注意维度方向（关键）
    % 如果 climchl 是 (lon,lat)，需要转置
    data = climchl0(:,:,j)';
    
    % ===== 经度-纬度插值 =====
    chll = interp2(Lon_004, Lat_004, data, ...
                   Lon_025, Lat_025, 'linear');
    
    climchla(:,:,j) = chll;
end  

save('Chla_clim.mat','Lon','Lat','climchla','-v7.3');


%%  画图看每日气候态叶绿素值   降低分辨率后的
clc,clear;
load Chla_clim.mat
climchl60=climchla(:,:,60);

figure('Position',[100 100 1000 1000]);  %绘图框
m_proj('equidistant','lon',[100 125],'lat',[0 25]);%投影范围
m_contourf(Lon,Lat,climchl60,500,'linestyle','none');
shading interp;
hold on;
m_gshhs_i('linewidth',1,'color','k');%生成地形
m_gshhs_i('patch',[0.7,0.7,0.7]);%陆地颜色灰度
m_grid('linestyle','none','box','on','linewidth',1,'fontsize',26);%生成经纬度网格

hold on
colormap(m_colmap('jet',18));
cb=colorbar;
caxis([0  0.8]);
set(cb,'tickdir','in'); % 刻度线朝内
set(cb,'linewidth',1,'fontsize',28,'edgecolor','k','FontName','Times New Roman');%设置colorbar字体大小和框线宽度
set(cb,'YTick',0:0.2:0.8); %色标值范围及显示间隔
set(get(cb,'ylabel'),'string','Chla(​mg/m^{3})','FontSize',30,'FontName','Times New Roman','FontWeight','bold');%生成colorbar标题，单位