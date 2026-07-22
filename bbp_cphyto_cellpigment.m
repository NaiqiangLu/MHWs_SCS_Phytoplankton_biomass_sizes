%% 查看nc文件信息
clc,clear;close all;
file1='L3m_19980125__576853982_4_GSM-SWF_BBP_DAY_00.nc';
fileinfo=ncinfo("L3m_19980101__576853982_4_GSM-SWF_BBP_DAY_00.nc");
BBP=ncread(file1,'BBP_mean');
lat=ncread(file1,'lat');  lon=ncread(file1,'lon');  

%% 读取bbp的nc数据为mat文件
clc; clear; close all;
% 1 设置路径
data_path = 'E:\第二篇\BBP_1998_2024\';

% 2 获取所有nc文件
filelist = dir(fullfile(data_path,'*.nc'));

% 按文件名排序（保证时间顺序）
[~,idx] = sort({filelist.name});
filelist = filelist(idx);

nfile = length(filelist);   % 应为9862个

% 3 读取经纬度（只读一次）
file1 = fullfile(data_path,filelist(1).name);

lat = ncread(file1,'lat');
lon = ncread(file1,'lon');

nlon = length(lon);
nlat = length(lat);

% 4 预分配BBP矩阵
BBP = nan(nlon,nlat,nfile,'double');

% 5 批量读取BBP
for i = 1:nfile
    
    file = fullfile(data_path,filelist(i).name);
    
    try
        bbp_day = ncread(file,'BBP_mean');
        
        BBP(:,:,i) = double(bbp_day);
        
    catch
        warning(['读取失败: ',filelist(i).name])
    end
    
    if mod(i,500)==0
        disp(['已读取 ',num2str(i),' / ',num2str(nfile)])
    end
    
end


% 7 保存mat文件
save('BBP_1998_2024.mat','BBP','lon','lat','-v7.3')

disp('BBP数据读取完成并保存为MAT文件')  %%保存了601*601*9862的高分辨率4km bbp.mat文件


%% 画图 BBP 4km
clc,clear;close all;
load BBP_1998_2024.mat
bbp90=nanmean(BBP(:,:,334:424),3); 

figure('pos',[10 10 1000 1000]);  %绘图框
m_proj('equidistant','lon',[100 125],'lat',[0 25]);  %投影范围
m_contourf(lon,lat,bbp90',300,'linestyle','none');
shading interp;
hold on;
m_gshhs_i('linewidth',1,'color','k');%生成地形  
m_gshhs_i('patch',[0.7,0.7,0.7]);             %陆地颜色灰度
m_grid('linestyle','none','box','on','linewidth',1,'fontsize',26);  %生成经纬度网格

hold on;
colormap(m_colmap('jet',20));
cb=colorbar;
caxis([0  0.01]);
set(cb,'tickdir','in');        %刻度线朝内
set(cb,'linewidth',1,'fontsize',28,'edgecolor','k','FontName','Times New Roman');   %设置colorbar字体大小和框线宽度
% set(cb,'YTick',0:0.2:0.8);   %色标值范围及显示间隔
set(get(cb,'ylabel'),'string','bbp','FontSize',30,'FontName','Times New Roman','FontWeight','bold');  %生成colorbar标题，单位


%% Cphyto 计算 + 4km → 0.25°（经纬度插值标准方法）

clear; clc; close all;

% ================== 1. 目标网格（0.25°） ==================
load sst_Lon_Lat   % Lon, Lat

% 区域裁剪
x2 = find(Lon>=100 & Lon<=125);
y2 = find(Lat>=0 & Lat<=25);

Lon = Lon(x2);
Lat = Lat(y2);

nlon = length(Lon);
nlat = length(Lat);

% 目标网格
[Lon025, Lat025] = meshgrid(Lon, Lat);

% ================== 2. 读取原始BBP ==================
m = matfile('BBP_1998_2024.mat');

lon = m.lon;
lat = m.lat;

% 原始网格裁剪
x1 = find(lon>=100 & lon<=125);
y1 = find(lat>=0 & lat<=25);

lon_cut = lon(x1);
lat_cut = lat(y1);

% 原始4km网格
[Lon04, Lat04] = meshgrid(lon_cut, lat_cut);

ntime = size(m,'BBP',3);

% ================== 3. 预分配 ==================
Cphyto = nan(nlat, nlon, ntime, 'single');  
% ⚠️ 注意顺序：lat × lon × time（与PAR保持一致）

% ================== 4. 主循环 ==================
for j = 1:ntime
    
    % ===== 读取一天 BBP =====
    bbp_day = m.BBP(x1, y1, j);   % (lon,lat)
    
    % ===== 转成 (lat,lon) =====
    bbp_day = bbp_day';  
    
    % ===== 数据质量控制 =====
    bbp_day(bbp_day < 0) = NaN;
    
    % ===== 计算 Cphyto =====
    cphyto_day = 13000 .* (bbp_day - 0.00098);
    cphyto_day(cphyto_day < 0) = NaN;
    
    % ===== 插值（核心改动：经纬度插值）=====
    tmp = interp2(Lon04, Lat04, cphyto_day, ...
                  Lon025, Lat025, 'linear');
    
    % ===== 存储 =====
    Cphyto(:,:,j) = single(tmp);
    
    % ===== 进度 =====
    if mod(j,500)==0
        fprintf('processing %d / %d\n', j, ntime);
    end
    
end

% ================== 5. 保存 ==================
save('Cphyto_025.mat', "Cphyto", "Lat", "Lon", "-v7.3");

disp('✅ Cphyto 0.25°（经纬度插值）计算完成');



%% 画图 Cphyto
clc,clear;close all;
load Cphyto_025.mat
Cphyto90=nanmean(Cphyto(:,:,334:424),3); 

load sst_Lon_Lat.mat

figure('pos',[10 10 1000 1000]);  %绘图框

m_proj('equidistant','lon',[100 125],'lat',[0 25]);  %投影范围
m_contourf(Lon,Lat,Cphyto90,100,'linestyle','none');
shading interp;
hold on;
m_gshhs_i('linewidth',1,'color','k');         %生成地形  
m_gshhs_i('patch',[0.7,0.7,0.7]);             %陆地颜色灰度
m_grid('linestyle','none','box','on','linewidth',1,'fontsize',26);%生成经纬度网格

hold on;
colormap(m_colmap('jet',20));
cb=colorbar;
caxis([0  30]);
set(cb,'tickdir','in');        %刻度线朝内
set(cb,'linewidth',1,'fontsize',28,'edgecolor','k','FontName','Times New Roman');   %设置colorbar字体大小和框线宽度
% set(cb,'YTick',0:0.2:0.8);   %色标值范围及显示间隔
set(get(cb,'ylabel'),'string','Cphyto','FontSize',30,'FontName','Times New Roman','FontWeight','bold');  %生成colorbar标题，单位

%% 计算cellular pigmentation  细胞色素沉着程度 θ
clear,clc,close all;
load Cphyto_025.mat
load Chla_19982024_interp025.mat
chla_full(chla_full<0 )=nan; %去掉异常chla值
cellpig=Cphyto./chla_full;                  %计算细胞色素沉着度

save('cellpigment_025.mat',"Lon",'Lat',"cellpig");  %保存数据

%% 画图cellular pigmentation  细胞色素沉着程度 θ
clc,clear;close all;
load cellpigment_025.mat
Cp90=nanmean(cellpig(:,:,334:424),3); 

load sst_Lon_Lat.mat

figure('pos',[100 100 1000 900]);  %绘图框
m_proj('equidistant','lon',[100 125],'lat',[0 25]);  %投影范围
m_contourf(Lon,Lat,Cp90,100,'linestyle','none');
shading interp;

hold on;
m_gshhs_i('linewidth',1,'color','k');%生成地形  
m_gshhs_i('patch',[0.7,0.7,0.7]);             %陆地颜色灰度
m_grid('linestyle','none','box','on','linewidth',1,'fontsize',26);%生成经纬度网格

hold on;
colormap(m_colmap('jet',20));
cb=colorbar;
caxis([0  140]);
set(cb,'tickdir','in');        %刻度线朝内
set(cb,'linewidth',1,'fontsize',30,'edgecolor','k','FontName','Times New Roman');   %设置colorbar字体大小和框线宽度
set(get(cb,'ylabel'),'string','Cp','FontSize',30,'FontName','Times New Roman','FontWeight','bold');  %生成colorbar标题，单位
