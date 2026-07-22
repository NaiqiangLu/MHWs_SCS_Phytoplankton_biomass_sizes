%% 计算1998–2024多年日平均Cphyto气候态（2月29日单独计算）
clc; clear; close all;
tic

% 数据路径
data_path = 'G:\第二篇\BBP_1998_2024';
load(fullfile(data_path,'Cphyto_025.mat')); 
% 变量: Cphyto (100×100×9862), Lon, Lat

[ni,nj,ndays] = size(Cphyto);
fprintf('Cphyto size: %d x %d x %d\n',ni,nj,ndays);

% 去除异常值（Cphyto不应为负）
Cphyto(Cphyto< 0) = NaN;

% 构造时间序列
dates_all = datetime(1998,1,1):datetime(2024,12,31);

if length(dates_all) ~= ndays
    error('时间长度与数据天数不一致');
end

% 初始化
Cphyto_sum   = zeros(ni,nj,366);
Cphyto_count = zeros(ni,nj,366);

fprintf('Start accumulating climatology...\n');

% 主循环
for d = 1:ndays

    date_now = dates_all(d);

    y = year(date_now);
    m = month(date_now);

    % 判断闰年
    leap = (mod(y,4)==0 & mod(y,100)~=0) | mod(y,400)==0;

    % 原始DOY
    doy = day(date_now,'dayofyear');

    % 平年3月之后整体后移一位
    if ~leap && m>=3
        doy = doy + 1;
    end

    data = Cphyto(:,:,d);

    valid = ~isnan(data);

    tmp = data;
    tmp(~valid) = 0;

    Cphyto_sum(:,:,doy) = Cphyto_sum(:,:,doy) + tmp;
    Cphyto_count(:,:,doy) = Cphyto_count(:,:,doy) + valid;

    if mod(d,365)==0
        fprintf('Processed %d / %d\n',d,ndays);
    end

end

% 计算气候态
Cphyto_clim = Cphyto_sum ./ Cphyto_count;
Cphyto_clim(Cphyto_count==0) = NaN;

% 保存（按你要求命名）
save(fullfile(data_path,'climCphyto.mat'),'Cphyto_clim','Lon','Lat','-v7.3');

fprintf('Cphyto climatology saved to climbbp.mat\n');

toc


%% 计算1998–2024多年日平均cellular pigmentation气候态（2月29日单独计算）
clc; clear; close all;
tic

% 数据路径
data_path = 'G:\第二篇\BBP_1998_2024';
load(fullfile(data_path,'cellpigment_025.mat')); 
% 变量: cellpig (100×100×9862), Lon, Lat

[ni,nj,ndays] = size(cellpig);
fprintf('cellpig size: %d x %d x %d\n',ni,nj,ndays);

% 去除异常值（一般cellpig不应为负）
cellpig(cellpig< 0) = NaN;

% 构造时间序列
dates_all = datetime(1998,1,1):datetime(2024,12,31);

if length(dates_all) ~= ndays
    error('时间长度与数据天数不一致');
end

% 初始化
cellpig_sum   = zeros(ni,nj,366);
cellpig_count = zeros(ni,nj,366);

fprintf('Start accumulating climatology...\n');

% 主循环
for d = 1:ndays

    date_now = dates_all(d);

    y = year(date_now);
    m = month(date_now);

    % 判断闰年
    leap = (mod(y,4)==0 & mod(y,100)~=0) | mod(y,400)==0;

    % 原始DOY
    doy = day(date_now,'dayofyear');

    % 平年3月以后整体后移一位（保证2月29日独立）
    if ~leap && m>=3
        doy = doy + 1;
    end

    data = cellpig(:,:,d);

    valid = ~isnan(data);

    tmp = data;
    tmp(~valid) = 0;

    cellpig_sum(:,:,doy) = cellpig_sum(:,:,doy) + tmp;
    cellpig_count(:,:,doy) = cellpig_count(:,:,doy) + valid;

    if mod(d,365)==0
        fprintf('Processed %d / %d\n',d,ndays);
    end

end

% 计算气候态
cellpig_clim = cellpig_sum ./ cellpig_count;
cellpig_clim(cellpig_count==0) = NaN;

% 保存
save(fullfile(data_path,'climCellpig.mat'),...
    'cellpig_clim','Lon','Lat','-v7.3');

fprintf('cellpig climatology saved to climcellpig.mat\n');

toc


%% 画图看看
clc,clear;
load climCellpig.mat
load sst_Lon_Lat.mat
cellpig_clim90=nanmean(cellpig_clim(:,:,335:366),3);

figure('pos',[100 100 900 800]);  %绘图框
m_proj('equidistant','lon',[105 123],'lat',[0 25]);%投影范围
m_pcolor(Lon,Lat,cellpig_clim90);
shading interp;
hold on;
m_gshhs_i('linewidth',1,'color','k');%生成地形
m_gshhs_i('patch',[0.7,0.7,0.7]);%陆地颜色灰度
m_grid('linestyle','none','box','on','linewidth',1,'fontsize',18);%生成经纬度网格

hold on
colormap(m_colmap('jet'));
cb=colorbar;
caxis([0  100]);
set(cb,'tickdir','in'); % 刻度线朝内
set(cb,'linewidth',1,'fontsize',22,'edgecolor','k','FontName','Times New Roman');%设置colorbar字体大小和框线宽度
set(get(cb,'ylabel'),'string','cellular pigmentation(θ)','FontSize',20,'FontName','Times New Roman','FontWeight','bold');%生成colorbar标题，单位