%% 查看数据信息，画图
clc,clear;
file='cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D_20030101_20030630.nc';
fileinfo=ncinfo(file);
lon=ncread(file,'longitude');
lat=ncread(file,'latitude');
micro=ncread(file,'MICRO');
nano=ncread(file,'NANO');
pico=ncread(file,'PICO');

% 2003年1月平均

micro=nanmean(micro(:,:,1:31),3);
nano=nanmean(nano(:,:,1:31),3);
pico=nanmean(pico(:,:,1:31),3);


% 画三组分（MICRO / NANO / PICO）2003.1月平均分布
figure('Position',[100 100 900 350]);

m_proj('mercator','lon',[105 125],'lat',[0 25]);

data_all = {micro, nano, pico};
name_used = {'MICRO','NANO','PICO'};

for i = 1:3
    subplot(1,3,i)

    m_pcolor(lon, lat, data_all{i}');
    shading interp

    m_gshhs_i('patch',[0.7 0.7 0.7]);
    m_grid('linestyle','none','box','on','tickdir','in');

    colormap(m_colmap('jet'));        % 叶绿素相关更常用
    colorbar('location','southoutside');
    caxis([0 0.4]);
    title([name_used{i}]);

end

%% PFTs 4km → 0.25° 重网格并拼接 1998–2024
clc; clear;close all;

%================== 1. 路径 ==================
dataDir = 'G:\第二篇\PFTs';
outFile = 'PFTs0.25_1998_2024.mat';

files = dir(fullfile(dataDir,...
'cmems_obs-oc_glo_bgc-plankton_my_l3-multi-4km_P1D_*.nc'));

fnames = sort({files.name});   % 保证时间顺序


% ================== 2. 读取0.25°网格 ==================
load sst_Lon_Lat.mat   % 含 Lon,Lat (100×1)

[Lon025,Lat025] = meshgrid(Lon,Lat);   % 100×100


% ================== 3. 读取4km原始网格 ==================
file0 = fullfile(dataDir,fnames{1});

lon4 = ncread(file0,'longitude');
lat4 = ncread(file0,'latitude');

[Lon4,Lat4] = meshgrid(lon4,lat4);  % 600×480  %meshgrid会转置网格


%================== 4. 计算总天数 ==================
total_days = 0;

for i = 1:numel(fnames)

    info = ncinfo(fullfile(dataDir,fnames{i}),'MICRO');
    total_days = total_days + info.Size(3);

end

fprintf('总天数 = %d\n', total_days);   % 应≈9862


% ================== 5. 预分配0.25°数组 ==================
nx = length(Lon);
ny = length(Lat);

MICRO = nan(nx,ny,total_days,'single');
NANO  = nan(nx,ny,total_days,'single');
PICO  = nan(nx,ny,total_days,'single');


%================== 6. 按文件读取+插值 ==================
t0 = 1;

for i = 1:numel(fnames)

    f = fullfile(dataDir,fnames{i});
    fprintf('Processing %s\n',fnames{i});

    micro = ncread(f,'MICRO');
    nano  = ncread(f,'NANO');
    pico  = ncread(f,'PICO');

    nt = size(micro,3);

    for t = 1:nt

        micro_day = squeeze(micro(:,:,t))';
        nano_day  = squeeze(nano(:,:,t))';
        pico_day  = squeeze(pico(:,:,t))';
        % 去除负值和0
        micro_day(micro_day < 0) = nan;
        nano_day(nano_day < 0)  = nan;
        pico_day(pico_day < 0)  = nan;

        % 0.25°插值
        MICRO(:,:,t0) = interp2(Lon4,Lat4,micro_day,...
            Lon025,Lat025,'linear');

        NANO(:,:,t0)  = interp2(Lon4,Lat4,nano_day,...
            Lon025,Lat025,'linear');

        PICO(:,:,t0)  = interp2(Lon4,Lat4,pico_day,...
            Lon025,Lat025,'linear');

        t0 = t0 + 1;

    end

end


%================== 7. 检查 ==================
if t0-1 ~= total_days
    error('时间维拼接出错');
end

disp('0.25° PFTs数据生成完成');

% ================== 8. 保存 ==================
save(outFile,...
    'MICRO','NANO','PICO','Lon','Lat','-v7.3');

disp('已保存 PFTs0.25_1998_2024.mat');

%% 查看降低分辨率后的图
clc,clear;close all;

load PFTs0.25_1998_2024.mat

% 2003年1月平均
micro=nanmean(MICRO(:,:,1827:1900),3);
nano=nanmean(NANO(:,:,1827:1900),3);
pico=nanmean(PICO(:,:,1827:1900),3);

% 画三组分（MICRO / NANO / PICO）2003.1月平均分布
figure('Position',[100 100 900 400]);

m_proj('mercator','lon',[100  125],'lat',[0 25]);

data_all = {micro, nano, pico};
name_used = {'MICRO','NANO','PICO'};

for i = 1:3
    subplot(1,3,i)

    m_pcolor(Lon, Lat, data_all{i});
    shading interp

    m_gshhs_i('patch',[0.7 0.7 0.7]);
    m_grid('linestyle','none','box','on','tickdir','in');

    colormap(m_colmap('jet'));       %叶绿素相关更常用
    colorbar('location','southoutside');
    caxis([0 0.35]);
    title([name_used{i}]);

end


%% 计算1998–2024多年日平均PFTs气候态（严格366天）
clc; clear; close all;
tic
% 数据路径
data_path = 'G:\第二篇\PFTs';
load(fullfile(data_path,'PFTs0.25_1998_2024.mat'));
% 变量: MICRO, NANO, PICO (100×100×9862), Lon, Lat

[ni,nj,ndays] = size(MICRO);
fprintf('Data size: %d x %d x %d\n',ni,nj,ndays);

% 去除异常值（不应为负）
MICRO(MICRO< 0) = NaN;
NANO(NANO< 0) = NaN;
PICO(PICO< 0) = NaN;

% 构造时间序列
dates_all = datetime(1998,1,1):datetime(2024,12,31);

if length(dates_all) ~= ndays
    error('时间长度与数据天数不一致');
end

% 初始化
MICRO_sum = zeros(ni,nj,366);
NANO_sum  = zeros(ni,nj,366);
PICO_sum  = zeros(ni,nj,366);

MICRO_count = zeros(ni,nj,366);
NANO_count  = zeros(ni,nj,366);
PICO_count  = zeros(ni,nj,366);

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

    % 平年3月以后整体后移（保证2月29日独立）
    if ~leap && m>=3
        doy = doy + 1;
    end

    % ===== MICRO =====
    data = MICRO(:,:,d);
    valid = ~isnan(data);
    tmp = data; tmp(~valid) = 0;
    MICRO_sum(:,:,doy) = MICRO_sum(:,:,doy) + tmp;
    MICRO_count(:,:,doy) = MICRO_count(:,:,doy) + valid;

    % ===== NANO =====
    data = NANO(:,:,d);
    valid = ~isnan(data);
    tmp = data; tmp(~valid) = 0;
    NANO_sum(:,:,doy) = NANO_sum(:,:,doy) + tmp;
    NANO_count(:,:,doy) = NANO_count(:,:,doy) + valid;

    % ===== PICO =====
    data = PICO(:,:,d);
    valid = ~isnan(data);
    tmp = data; tmp(~valid) = 0;
    PICO_sum(:,:,doy) = PICO_sum(:,:,doy) + tmp;
    PICO_count(:,:,doy) = PICO_count(:,:,doy) + valid;

    if mod(d,365)==0
        fprintf('Processed %d / %d\n',d,ndays);
    end

end

% 计算气候态
clim_MICRO = MICRO_sum ./ MICRO_count;
clim_NANO  = NANO_sum  ./ NANO_count;
clim_PICO  = PICO_sum  ./ PICO_count;

clim_MICRO(MICRO_count==0) = NaN;
clim_NANO(NANO_count==0)   = NaN;
clim_PICO(PICO_count==0)   = NaN;

% 保存
save(fullfile(data_path,'climPFTs.mat'),...
    'clim_MICRO','clim_NANO','clim_PICO','Lon','Lat','-v7.3');

fprintf('PFTs climatology saved to climPFTs.mat\n');
toc

%% 查看计算出1998-2024年每日气候态PFTs的图
clc,clear;close all;

load climPFTs.mat;

% 取出某一天气候态的值
micro=clim_MICRO(:,:,61);
nano=clim_NANO(:,:,61);
pico=clim_PICO(:,:,61);

% 画三组分（MICRO / NANO / PICO）2003.1月平均分布
figure('Position',[100 100 900 500]);

m_proj('mercator','lon',[105 125],'lat',[0 25]);

data_all = {micro, nano, pico};
name_used = {'MICRO','NANO','PICO'};

for i = 1:3
    subplot(1,3,i)

    m_pcolor(Lon, Lat, data_all{i});
    shading interp

    m_gshhs_i('patch',[0.7 0.7 0.7]);
    m_grid('linestyle','none','box','on','tickdir','in');

    colormap(m_colmap('jet'));       %叶绿素相关更常用
    colorbar('location','southoutside');
    caxis([0 0.3]);
    title([name_used{i}]);

end

%% 计算PFTs 三种类型植物的比例
clc,clear;close all
load PFTs0.25_1998_2024.mat

Tphyto = sum(cat(4, MICRO, NANO, PICO), 4, 'omitnan');  %把三种组分植物加起来得到总的浓度

% 避免后面除以0
Tphyto(Tphyto==0) = NaN;  

%算三种分粒径的比例

Fmicro = MICRO ./ Tphyto;
Fnano  = NANO  ./ Tphyto;
Fpico  = PICO  ./ Tphyto;

save('Fphyto1998-2024.mat',"Fpico","Fnano","Fmicro","Lon","Lat");

%% 每日气候态的三组分粒径植物比例
clc,clear,close all
load climPFTs.mat

% ===== 总量 =====
Tclim = sum(cat(4, clim_MICRO, clim_NANO, clim_PICO), 4, 'omitnan');

% ===== 避免除以0 =====
Tclim(Tclim==0) = NaN;

% ===== 占比 =====
climFmicro = clim_MICRO ./ Tclim;
climFnano  = clim_NANO  ./ Tclim;
climFpico  = clim_PICO  ./ Tclim;

save('ClimFphyto1998-2024.mat',"climFmicro","climFnano","climFpico");
