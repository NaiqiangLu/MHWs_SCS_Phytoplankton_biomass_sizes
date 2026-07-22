%% 计算1998–2024多年日平均风场气候态（u, v, ws，含2月29日）
clc; clear; close all;
tic

% ================== 1. 数据路径 ==================
data_path = 'E:\第二篇\wind1998-2024';
load(fullfile(data_path,'CCMP_wind1998_2024.mat'));
% 变量：
% u (100×100×9862)
% v (100×100×9862)
% ws (100×100×9862)


[ni,nj,ndays] = size(u);
fprintf('Wind data size: %d x %d x %d\n',ni,nj,ndays);

% ================== 2. 异常值处理 ==================
% 风速不应为负（保险处理）
u(abs(u)>100)   = NaN;
v(abs(v)>100)   = NaN;
ws(ws<0)        = NaN;

% ================== 3. 构造时间 ==================
dates_all = datetime(1998,1,1):datetime(2024,12,31);

if length(dates_all) ~= ndays
    error('时间长度与数据天数不一致');
end

% ================== 4. 初始化 ==================
u_sum   = zeros(ni,nj,366);
v_sum   = zeros(ni,nj,366);
ws_sum  = zeros(ni,nj,366);

u_count  = zeros(ni,nj,366);
v_count  = zeros(ni,nj,366);
ws_count = zeros(ni,nj,366);

fprintf('Start accumulating climatology...\n');

% ================== 5. 主循环 ==================
for d = 1:ndays

    date_now = dates_all(d);

    y = year(date_now);
    m = month(date_now);

    % 判断闰年
    leap = (mod(y,4)==0 & mod(y,100)~=0) | mod(y,400)==0;

    % 原始DOY
    doy = day(date_now,'dayofyear');

    % 平年3月后 +1（保证366框架）
    if ~leap && m>=3
        doy = doy + 1;
    end

    % ===== 取数据 =====
    u_data  = u(:,:,d);
    v_data  = v(:,:,d);
    ws_data = ws(:,:,d);

    % ===== 有效值 =====
    valid_u  = ~isnan(u_data);
    valid_v  = ~isnan(v_data);
    valid_ws = ~isnan(ws_data);

    % ===== NaN置0用于累加 =====
    u_tmp = u_data;   u_tmp(~valid_u) = 0;
    v_tmp = v_data;   v_tmp(~valid_v) = 0;
    ws_tmp = ws_data; ws_tmp(~valid_ws) = 0;

    % ===== 累加 =====
    u_sum(:,:,doy)  = u_sum(:,:,doy)  + u_tmp;
    v_sum(:,:,doy)  = v_sum(:,:,doy)  + v_tmp;
    ws_sum(:,:,doy) = ws_sum(:,:,doy) + ws_tmp;

    u_count(:,:,doy)  = u_count(:,:,doy)  + valid_u;
    v_count(:,:,doy)  = v_count(:,:,doy)  + valid_v;
    ws_count(:,:,doy) = ws_count(:,:,doy) + valid_ws;

    if mod(d,365)==0
        fprintf('Processed %d / %d\n',d,ndays);
    end
end

% ================== 6. 计算气候态 ==================
u_clim  = u_sum ./ u_count;
v_clim  = v_sum ./ v_count;
ws_clim = ws_sum ./ ws_count;

u_clim(u_count==0)   = NaN;
v_clim(v_count==0)   = NaN;
ws_clim(ws_count==0) = NaN;

% ================== 7. 保存 ==================
save(fullfile(data_path,'climWind_1998_2024.mat'), ...
    'u_clim','v_clim','ws_clim','Lon','Lat','-v7.3');

fprintf('Wind climatology saved successfully!\n');

toc
