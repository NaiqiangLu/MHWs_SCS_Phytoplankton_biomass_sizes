function result = mk_test(x, alpha)
% MANN-KENDALL TREND TEST + SEN'S SLOPE
%
% 输入：
%   x     - 时间序列（列向量或行向量）
%   alpha - 显著性水平（默认 0.05）
%
% 输出（结构体 result）：
%   result.trend      - 趋势方向 ('increasing','decreasing','no trend')
%   result.h          - 是否显著 (1=显著, 0=不显著)
%   result.p          - p值
%   result.Z          - 标准统计量
%   result.S          - MK统计量
%   result.slope      - Sen's slope (每年变化率)
%   result.slope_decade - decade趋势
%
% 作者：适用于SCI绘图分析

if nargin < 2
    alpha = 0.05;
end

x = x(:);
n = length(x);

%% ===== 1. 计算 S =====
S = 0;
for i = 1:n-1
    for j = i+1:n
        S = S + sign(x(j) - x(i));
    end
end

%% ===== 2. 方差（无重复值简化版）=====
varS = n*(n-1)*(2*n+5)/18;

%% ===== 3. Z值 =====
if S > 0
    Z = (S - 1)/sqrt(varS);
elseif S < 0
    Z = (S + 1)/sqrt(varS);
else
    Z = 0;
end

%% ===== 4. p值 =====
p = 2 * (1 - normcdf(abs(Z),0,1));
h = p < alpha;

%% ===== 5. Sen's slope =====
slopes = [];
for i = 1:n-1
    for j = i+1:n
        slopes(end+1) = (x(j) - x(i)) / (j - i);
    end
end

sen_slope = median(slopes);

%% ===== 6. 趋势判断 =====
if h == 1
    if sen_slope > 0
        trend = 'increasing';
    else
        trend = 'decreasing';
    end
else
    trend = 'no trend';
end

%% ===== 7. 输出 =====
result.trend = trend;
result.h = h;
result.p = p;
result.Z = Z;
result.S = S;
result.slope = sen_slope;               % 每年变化
result.slope_decade = sen_slope * 10;   % decade
end
