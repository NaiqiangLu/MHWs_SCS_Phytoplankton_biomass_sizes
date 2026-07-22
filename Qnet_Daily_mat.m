clc,clear,close all;
% 潜热通量 （负值）
file_slhf='SLHF2009_daily.nc';
file_slhfin=ncinfo (file_slhf);
SLHF=ncread(file_slhf,'slhf');
SLHF= SLHF(:,:,215); %读取这一年某天的SLHF

SLHF=SLHF./3600; %转换单位 从 J*m^-2 转换为 W*m^-2

% 感热通量  （负值）
file_sshf='SSHF2009_daily.nc'; 
file_sshfin=ncinfo (file_sshf);
SSHF=ncread(file_sshf,'sshf');
SSHF= SSHF(:,:,215); %读取这一年某天的SSHF

SSHF=SSHF./3600; %转换单位 从 J*m^-2 转换为 W*m^-2

% 海洋长波辐射量 （负值）
file_str='STR2009_daily.nc'; 
file_strin=ncinfo (file_str); 
STR=ncread(file_str,'str');
STR= STR(:,:,215);%读取这一年某天的STR

STR=STR./3600; %转换单位 从 J*m^-2 转换为 W*m^-2

% 太阳短波辐射量  （正值）
file_ssr='SSR2009_daily.nc'; 
file_ssrin=ncinfo (file_ssr);
SSR=ncread(file_ssr,'ssr');
SSR= SSR(:,:,215);  %读取这一年某天的SSR

SSR=SSR./3600; %转换单位 从 J*m^-2 转换为 W*m^-2

Qnet=SSR+SLHF+STR+SSHF;

%% 热通量数据：拼接 + 经纬度插值到0.25°
clc; clear; close all;

data_path = 'G:\第二篇\1998-2024月(daily)平均热通量和降水';
years = 1998:2024;

% ===== 目标网格（0.25°，来自SST）=====
load sst_Lon_Lat
x2 = find(Lon>=100 & Lon<=125);
y2 = find(Lat>=0   & Lat<=25);
Lon = Lon(x2);
Lat = Lat(y2);

[Lon_025, Lat_025] = meshgrid(Lon, Lat);

% ===== 初始化 =====
SSR  = [];
STR  = [];
SLHF = [];
SSHF = [];

for yr = years
    
    disp(['Processing year: ', num2str(yr)]);
    
    % ===== 文件名 =====
    file_ssr  = fullfile(data_path, ['SSR',  num2str(yr), '_daily.nc']);
    file_str  = fullfile(data_path, ['STR',  num2str(yr), '_daily.nc']);
    file_slhf = fullfile(data_path, ['SLHF', num2str(yr), '_daily.nc']);
    file_sshf = fullfile(data_path, ['SSHF', num2str(yr), '_daily.nc']);
    
    % ===== 读取数据 =====
    ssr  = ncread(file_ssr,  'ssr');
    str  = ncread(file_str,  'str');
    slhf = ncread(file_slhf, 'slhf');
    sshf = ncread(file_sshf, 'sshf');
    
    % ===== 经纬度（只需第一年读取一次）=====
    if yr == years(1)
        lat01 = ncread(file_ssr, 'latitude');
        lon01 = ncread(file_ssr, 'longitude');
        
        x1 = find(lon01>=100 & lon01<=125);
        y1 = find(lat01>=0   & lat01<=25);
        
        lon01 = lon01(x1);
        lat01 = lat01(y1);
        
        [Lon_021, Lat_021] = meshgrid(lon01, lat01);
    end
    
    % ===== 裁剪原始数据 =====
    ssr  = ssr(x1,y1,:);
    str  = str(x1,y1,:);
    slhf = slhf(x1,y1,:);
    sshf = sshf(x1,y1,:);
    
    % ===== 单位转换：J/m² → W/m² =====
    ssr  = ssr  / 3600;
    str  = str  / 3600;
    slhf = slhf / 3600;
    sshf = sshf / 3600;
    
    % ===== 转单精度 =====
    ssr  = single(ssr);
    str  = single(str);
    slhf = single(slhf);
    sshf = single(sshf);
    
    % ===== 初始化插值后的数组 =====
    nt = size(ssr,3);
    ssr_i  = nan(length(Lat), length(Lon), nt, 'single');
    str_i  = nan(length(Lat), length(Lon), nt, 'single');
    slhf_i = nan(length(Lat), length(Lon), nt, 'single');
    sshf_i = nan(length(Lat), length(Lon), nt, 'single');
    
    % ===== 每天做经纬度插值 =====
    for t = 1:nt
        
        % ⚠️ 如果原始数据是 (lon,lat)，需要转置
        ssr_i(:,:,t)  = interp2(Lon_021, Lat_021, ssr(:,:,t)',  Lon_025, Lat_025, 'linear');
        str_i(:,:,t)  = interp2(Lon_021, Lat_021, str(:,:,t)',  Lon_025, Lat_025, 'linear');
        slhf_i(:,:,t) = interp2(Lon_021, Lat_021, slhf(:,:,t)', Lon_025, Lat_025, 'linear');
        sshf_i(:,:,t) = interp2(Lon_021, Lat_021, sshf(:,:,t)', Lon_025, Lat_025, 'linear');
        
    end
    
    % ===== 拼接 =====
    SSR  = cat(3, SSR,  ssr_i);
    STR  = cat(3, STR,  str_i);
    SLHF = cat(3, SLHF, slhf_i);
    SSHF = cat(3, SSHF, sshf_i);
    
end

% ===== 计算净热通量 =====
Qnet = SSR + STR + SLHF + SSHF;

% ===== 转单精度（经纬度）=====
Lon = single(Lon);
Lat = single(Lat);

% ===== 检查 =====
disp(size(Qnet))   % 应为 (lon × lat× 9862)

% ===== 保存 =====
save('Qnet_interp025_1998_2024.mat', ...
    'Qnet','SSR','STR','SLHF','SSHF','Lon','Lat','-v7.3');

disp('✅ 拼接 + 插值 + 统一网格 已完成');



%% 画图查看
clc,clear;close all;
load Qnet_interp025_1998_2024.mat

Qnetx=Qnet(:,:,225);

figure('pos',[100 100 900 800]);  %绘图框
m_proj('equidistant','lon',[105 123],'lat',[0 25]);%投影范围
m_contourf(Lon,Lat,Qnetx,50,'linestyle','none');  %这里画图不需要转置
shading interp;

hold on;
m_gshhs_i('linewidth',1,'color','k');%生成地形
m_gshhs_i('patch',[0.7,0.7,0.7]);%陆地颜色灰度
m_grid('linestyle','none','box','on','linewidth',1,'fontsize',18);%生成经纬度网格

hold on;
colormap(m_colmap('diverging',14));
cb=colorbar('southoutside');
caxis([-50  200]);
set(cb,'tickdir','in'); % 刻度线朝内
set(cb,'linewidth',1,'fontsize',18,'edgecolor','k','FontName','Times New Roman');%设置colorbar字体大小和框线宽度
set(get(cb,'ylabel'),'string','Qnet(W m^{-2})','FontSize',26,'FontName','Times New Roman','FontWeight','bold');%生成colorbar标题，单位


%%%% chla ，热通量的数据是lon*lat*time  经度从小到大排列  纬度从大到小排列，画图的时候需要矩阵顺序要匹配
%%%% 我的wind sst 以及用sst做经纬度插值的数据 其经纬度都是从小到大排列 chla经过经纬度降分辨率插值后已经和sst数据保持一致方向


SLHF1=SLHF(:,:,225);  
SSHF1=SSHF(:,:,225);
SSR1=SSR(:,:,225);
STR1=STR(:,:,225);
