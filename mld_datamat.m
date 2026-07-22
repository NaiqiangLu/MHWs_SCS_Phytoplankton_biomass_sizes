clc,clear;close all;
file='cmems_mod_glo_phy_my_0.083deg_P1D-m_19980101_20031231.nc';
finfo=ncinfo(file);

lon=ncread(file,'longitude'); 
lat=ncread(file,'latitude');

x=find(lon>=100 & lon<=125); 
y=find(lat>=0 & lat<=25);

lon=lon(x); lat=lat(y); 

mld=ncread(file,'mlotst');%读取混合层nc数据

mld=mld(x,y,1); %限制范围




%% MLD 0.083° → 0.25° 重网格并拼接 1998–2024
clc; clear; close all;

%================== 1. 路径 ==================
dataDir = 'E:\第二篇\ML_irradiance\混合层深度daily_1998_2024';
outFile = 'mld1998-2024_0.25.mat';

files = dir(fullfile(dataDir,'cmems_mod_glo_phy_my_0.083deg_P1D-m_*.nc'));
fnames = sort({files.name});   % 时间排序


%================== 2. 读取0.25°目标网格 ==================
load sst_Lon_Lat.mat   % Lon, Lat (100×1)

[Lon025,Lat025] = meshgrid(Lon,Lat);   % 100×100


%================== 3. 读取原始0.083°网格 ==================
file0 = fullfile(dataDir,fnames{1});

lon = ncread(file0,'longitude');
lat = ncread(file0,'latitude');

% 裁剪区域
x = find(lon>=100 & lon<=125);
y = find(lat>=0 & lat<=25);

lon = lon(x);
lat = lat(y);

[Lon083,Lat083] = meshgrid(lon,lat);   % 注意meshgrid


%================== 4. 计算总天数 ==================
total_days = 0;

for i = 1:numel(fnames)

    info = ncinfo(fullfile(dataDir,fnames{i}),'mlotst');
    total_days = total_days + info.Size(3);

end

fprintf('总天数 = %d\n', total_days);   % 应≈9862


%================== 5. 预分配 ==================
nx = length(Lon);
ny = length(Lat);

MLD = nan(nx,ny,total_days,'single');


%================== 6. 逐文件读取 + 插值 ==================
t0 = 1;

for i = 1:numel(fnames)

    f = fullfile(dataDir,fnames{i});
    fprintf('Processing %s\n',fnames{i});

    mld = ncread(f,'mlotst');   % (lon,lat,time)

    % 裁剪区域
    mld = mld(x,y,:);

    nt = size(mld,3);

    for t = 1:nt

        % 转置到 (lat,lon)
        mld_day = squeeze(mld(:,:,t))';

        % 去异常值（MLD一般>=0）
        mld_day(mld_day < 0) = nan; %没有小于0的混合层

        % 插值到0.25°
        MLD(:,:,t0) = interp2(Lon083,Lat083,mld_day,...
                             Lon025,Lat025,'linear');

        t0 = t0 + 1;

    end

end


%================== 7. 检查 ==================
if t0-1 ~= total_days
    error('时间拼接出错');
end

disp('0.25° MLD数据生成完成');


%================== 8. 保存 ==================
save(outFile,'MLD','Lon','Lat','-v7.3');

disp('已保存 mld1998-2024_0.25.mat');





%% 查看降低分辨率后的图
clc,clear;close all;

load mld1998_2024_0.25.mat

% 2003年1月平均
mld=nanmean(MLD(:,:,335:365),3);

figure('Position',[100 100 900 800]);

m_proj('mercator','lon',[105 125],'lat',[0 25]);

m_pcolor(Lon, Lat, mld);  %mld不用转置，原始mld为lon*lat，降低分辨率插值方向是lat*lon，但在插值之前mld矩阵已经提前转置过，因此画图不需要转置
shading interp;

m_gshhs_i('patch',[0.7 0.7 0.7]);
m_grid('linestyle','none','box','on','tickdir','in');

colormap(m_colmap('jet'));       %叶绿素相关更常用
colorbar('location','southoutside');

title('MLD');

