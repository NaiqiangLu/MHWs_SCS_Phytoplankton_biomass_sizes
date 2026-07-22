%% 计算 mixed-layer averaged irradiance
clc; clear; close all;

% ================== 1. 读取数据 ==================
load sst_Lon_Lat.mat          % Lon (100x1), Lat (100x1)
load mld1998_2024_0.25.mat    % MLD (100x100x9862)
load PAR1998_2024_0.25.mat    % PAR (100x100x9862)
load KD1998_2024_0.25.mat     % Kd  (100x100x9862)

% ================== 2. 预分配 ==================
[ni, nj, nt] = size(MLD);
Ibar = nan(ni, nj, nt);   % mixed-layer averaged irradiance

% ================== 3. 参数处理 ==================
% 防止除0或数值爆炸
Kd(Kd <= 0) = NaN;
MLD(MLD <= 0) = NaN;

% ================== 4. 核心计算 ==================
% 公式：
% Ibar = (1 / (kd * zml)) * I0 * exp(-kd * zml) * (1 - exp(-kd * zml))

for t = 1:nt
    kd  = Kd(:,:,t);
    zml = MLD(:,:,t);
    I0  = PAR(:,:,t);
    
    % 计算 kd*zml
    kz = kd .* zml;
    
    % 为避免 kz 很小导致数值不稳定，做一个判断
    small_mask = kz < 1e-6;
    
    tmp = nan(ni,nj);
    
    % 正常情况
    tmp(~small_mask) = (1 ./ kz(~small_mask)) .* I0(~small_mask) .* ...
        exp(-kz(~small_mask)) .* (1 - exp(-kz(~small_mask)));
    
    % kz → 0 时的极限（泰勒展开）
    % Ibar ≈ I0 * (1 - kz/2)
    tmp(small_mask) = I0(small_mask) .* (1 - kz(small_mask)/2);
    
    Ibar(:,:,t) = tmp;
end

% ================== 5. 保存 ==================
save('Ibar1998_2024.mat','Ibar','Lon','Lat','-v7.3');

disp('Mixed-layer averaged irradiance 计算完成');



%% 查看ibar图
clc,clear;close all;

load Ibar1998_2024.mat

% 2003年1月平均
Ibar=nanmean(Ibar(:,:,200:365),3);

figure('Position',[100 100 900 800]);

m_proj('miller','lon',[105 125],'lat',[0 25]);

m_pcolor(Lon, Lat, Ibar);  %kd不用转置，原始mld为lon*lat，降低分辨率插值方向是lat*lon，但在插值之前mld矩阵已经提前转置过，因此画图不需要转置
shading interp;

m_gshhs_i('patch',[0.7 0.7 0.7]);
m_grid('linestyle','none','box','on','tickdir','in');

colormap(nclCM(156,20));       %叶绿素相关更常用
colorbar('location','southoutside');
caxis([0, 25]);   %范围0-20
title('Mixed-layer averaged irradiance');