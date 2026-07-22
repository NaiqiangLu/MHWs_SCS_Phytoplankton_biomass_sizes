function clim = daily_climatology(data, start_date, end_date)
% =========================================================
% 通用函数：计算1998–2024逐日气候态（366天框架）
%
% 输入：
%   data        : (ni × nj × nt) 数据
%   start_date  : datetime，如 datetime(1998,1,1)
%   end_date    : datetime，如 datetime(2024,12,31)
%
% 输出：
%   clim        : (ni × nj × 366) 每日气候态
%
% 特点：
%   ✔ 自动处理闰年
%   ✔ 自动跳过 NaN
%   ✔ 适用于任意变量（MLD / Chl / PAR / Wind等）
% =========================================================

[ni,nj,nt] = size(data);

% ===== 构造时间 =====
dates_all = start_date:end_date;

if length(dates_all) ~= nt
    error('时间长度与数据不一致');
end

% ===== 初始化 =====
data_sum   = zeros(ni,nj,366);
data_count = zeros(ni,nj,366);

fprintf('Start calculating daily climatology...\n');

% ===== 主循环 =====
for d = 1:nt
    
    date_now = dates_all(d);
    
    y = year(date_now);
    m = month(date_now);

    % ===== 判断闰年 =====
    leap = (mod(y,4)==0 & mod(y,100)~=0) | mod(y,400)==0;

    % ===== DOY =====
    doy = day(date_now,'dayofyear');

    % ===== 平年补位（核心）=====
    if ~leap && m>=3
        doy = doy + 1;
    end

    % ===== 取数据 =====
    tmp = data(:,:,d);

    % ===== 有效值 =====
    valid = ~isnan(tmp);

    % ===== NaN置0 =====
    tmp(~valid) = 0;

    % ===== 累加 =====
    data_sum(:,:,doy)   = data_sum(:,:,doy)   + tmp;
    data_count(:,:,doy) = data_count(:,:,doy) + valid;

    % ===== 进度 =====
    if mod(d,365)==0
        fprintf('Processed %d / %d\n',d,nt);
    end
end

% ===== 计算平均 =====
clim = data_sum ./ data_count;
clim(data_count==0) = NaN;

fprintf('Climatology done!\n');

end
