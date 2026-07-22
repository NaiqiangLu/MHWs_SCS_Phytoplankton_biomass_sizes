%% 1. Loading data
clc,clear; close all;
% Load NOAA OI SST V2 data
sst_full=NaN(100,100,datenum(2025,12,31)-datenum(1982,1,1)+1);  %我拥有的sst数据长度
for i=1982:2025;
    file_here=['sst_' num2str(i)];
    load(file_here);
    eval(['data_here=sst_' num2str(i) ';'])  %eval:执行字符串中的代码并输出到工作区  这里是把每一年的sst数据按序读取
    sst_full(:,:,(datenum(i,1,1):datenum(i,12,31))-datenum(1982,1,1)+1)=data_here;  %这里把1982-2025所有年份的天数拼接在第三维后 把日期转换为天数 相减得到16071天
end

sst_full(sst_full < 0 | sst_full > 45) = NaN;  %把温度FillValue去掉

% save('sst_19822025.mat',"sst_full");
% 1982 to 2025.

load('sst_Lon_Lat'); 
size(sst_full); %size of data
datenum(2025,12,31)-datenum(1982,1,1)+1; % The temporal length is 44 years.

%% 2. Detecting MHWs and MCSs

% Here we detect marine heatwaves SCS based on the
% traditional definition of MHWs (Hobday et al. 2016). We detected MHWs
% during 1998 to 2024 for climatologies and thresholds in 1995 to 2024.
tic
[MHW,mclim,m90,mhw_ts]=detect(sst_full,datenum(1982,1,1):datenum(2025,12,31),datenum(1995,1,1),datenum(2024,12,31),datenum(1982,1,1),datenum(2024,12,31)); 
toc   %第一个时间是海温数据时间   第二个时间是气候态基期的时间  第三个时间是要检测的热浪事件时间
%这里调用detect函数，注意输入量，输入量是按需要是可变的，不需要全部的变量

% [MCS,~,m10,mcs_ts]=detect(sst_full,datenum(1982,1,1):datenum(2016,12,31),datenum(1982,1,1),datenum(2005,12,31),datenum(1993,1,1),datenum(2016,12,31),'Event','MCS','Threshold',0.1);

% Have a look of these two data.看前5行的事件
MHW(1:5,:);
% MCS(1:5,:);

% You could see that the properties `mhw_onset` and `mhw_end` are in a
% strange format. This is due to the fact that they are originally
% constructed in numeric format. We could change it to date format by
% following steps.日期转换，方便查看

datevec(num2str(MHW{1:5,:}),'yyyymmdd') % vector
datestr(datevec(num2str(MHW{1:5,:}),'yyyymmdd')) % string
datenum(num2str(MHW{1:5,:}),'yyyymmdd') % number

%% 3. Visualizing MHW/MCS time series

% Have a look of MHW events in grid (1,2) from Sep 2015 to Apr 2016
% figure('pos',[10 10 1000 1000]);
% event_line(sst_full,MHW,mclim,m90,[1 14],1982,[2023 1 1],[2024 12 31]);
% 数据长度，mhw事件表，气候态基期，90%阈值，格点位置，检测的事件长度
% % Have a look of MCS events in grid (1,2) during 1994
% figure('pos',[10 10 1000 1000]);
% event_line(sst_full,MCS,mclim,m10,[1 2],1982,[1994 1 1],[1994 12 31],'Event','MCS','Color',[0.5 0.5 1]);

%% 4. Mean states and trends

% Now we would like to know the mean states and annual trends of MHW
% frequency, i.e. how many MHW events would be detected per year and how it
% changes with time.

% [mean_freq,annual_freq,trend_freq,p_freq]=mean_and_trend(MHW,mhw_ts,1982,'Metric','Frequency');  %只计算mhw频率方面的多年平均和年度平均值
% 调用mean_and_trend函数要注意，年份：datastart是指你设定的mhw检测时期的起始年份，而不是sst数据源的起始年份
%这句代码是调用mean_and_trend函数仅计算MHWs频率方面的多年均值，每年频率，趋势和显著性

% These four outputs separately represent the total mean, annual mean,
% annual trend and associated p value of frequency.

% This function could detect mean states and trends for six different
% variables (Frequency, mean intensity, max intensity, duration and total
% MHW/MCs days). 

metric_used={'Frequency','MeanInt','MaxInt','CumInt','Duration','Days'}; %选取mhw所有参数

for i=1:6
    eval(['[mean_' metric_used{i} ',annual_' metric_used{i} ',trend_' metric_used{i} ',p_' metric_used{i} ']=mean_and_trend(MHW,mhw_ts,1982,' '''' 'Metric' '''' ',' 'metric_used{i}' ');'])
end   
%计算6个参数的mean_metric：多年平均指标，研究时段内（如 1998–2024）的 平均 MHW 指标。 这里的年份是指mhw_ts的起始时间，检测热浪的起始时间
% annual_metric：逐年指标  每一年计算一次 MHW 指标。


%% ================== 0. 数据整理  画图自定义属性 ==================
% plot mean and trend

for i = 1:6
    mean_data{i}  = eval(['mean_'  metric_used{i}]);
    trend_data{i} = eval(['trend_' metric_used{i}]);
    p_data{i}     = eval(['p_'     metric_used{i}]);
end


%================== Figure 1: Mean（tiledlayout版） ==================
figure('Position',[50 50 1100 900],'Color','w');

%  紧凑布局
tl = tiledlayout(2,3,'TileSpacing','tight','Padding','loose');

% ===== 每个子图 colorbar 范围 =====
clim_mean = [
    1   3.3;
    0.5 2.1;
    1   5.8;
    7   26;
    7   15;
    10   35 ];

% ===== 子图标题 =====
title_mean = {'(a) Frequency (Count)','(b) MeanInt (℃)','(c) MaxInt (℃)',...
    '(d) CumInt (℃ days)','(e) Duration (Days)','(f) Total Days (Days)'};

for i = 1:6
    
    nexttile   
    
    % ===== 数据 =====
    mean_here = mean_data{i};
    
    %  m_map 投影必须每个子图单独设置
    m_proj('mercator','lon',[105 123],'lat',[0 25]);
    
    % ===== 作图 =====
    m_pcolor(Lon,Lat,mean_here');
    shading interp
    hold on
    
    m_gshhs_i('patch',[0.7 0.7 0.7]);
    m_grid('linestyle','none','box','on','tickdir','in','fontsize',14,'xtick',[106:5:123],'Linewidth',1.1);
    
    colormap(m_colmap('diverging',32));
    
    %  每个子图独立范围
    caxis(clim_mean(i,:));

    %  colorbar（更紧凑）
    cb = colorbar('southoutside');
    cb.FontSize = 14;
    cb.LineWidth =1.2;

    %wscs上升流系统区域
    hold on
    % ===== 红色区域框 =====
    lon_box = [107.5 113.5 113.5 107.5 107.5];
    lat_box = [8.5   8.5   13.5   13.5   8.5];
    m_line(lon_box, lat_box, 'color','red','linewidth',1.3);

    %luzon上升流系统区域
    hold on
    % ===== 红色区域框 =====
    lon_box = [117   121    121 117  117];
    lat_box = [16.5  16.5   21  21  16.5];
    m_line(lon_box, lat_box, 'color','k','linewidth',1.3);
    % 标题
    title(title_mean{i},'fontsize',16,'fontname','Times New Roman');
    
    hold off
end


%% ================== Figure 2: Trend（tiledlayout版） ==================
[Lon2 , Lat2]=meshgrid(Lon, Lat);
figure('Position',[50 50 1100 900],'Color','w');

%  紧凑布局（核心）
tl = tiledlayout(2,3,'TileSpacing','tight','Padding','loose');

% ===== colorbar范围 =====
clim_trend = [
    -0.12 0.23;
    -0.03 0.025;
    -0.04 0.05;
    -0.81 1.22;
    -0.45 1.07;
    -0.91 3.42
];

% ===== 标题 =====
title_trend = {'(a) Frequency (Count/year)','(b) MeanInt (℃/year)','(c) MaxInt (℃/year)',...
    '(d) CumInt (℃ days/year)','(e) Duration (Days/year)','(f) Total days (Days/year)'};

for i = 1:6
    
    nexttile   
    
    % ===== 数据 =====
    trend_here = trend_data{i};
    p_here     = p_data{i};
    
    %  m_map必须在每个子图重新设置
    m_proj('mercator','lon',[105 123],'lat',[0 25]);
    
    % ===== 作图 =====
    m_pcolor(Lon,Lat,trend_here');
    shading interp
    hold on
    
    m_gshhs_i('patch',[0.7 0.7 0.7]);
    m_grid('linestyle','none','box','on','tickdir','in','fontsize',14,'xtick',[106:5:123],'linewidth',1.1);
    
    colormap(m_colmap('diverging',32));
    
    %  每个子图独立范围
    caxis(clim_trend(i,:))
    
    %  colorbar（更紧凑）
    cb = colorbar('southoutside');
    cb.FontSize = 14;
    cb.LineWidth =1.2;

    %wscs上升流系统区域
    hold on
    % ===== 红色区域框 =====
    lon_box = [107.5 113.5 113.5 107.5 107.5];
    lat_box = [8.5   8.5   13.5   13.5   8.5];
    m_line(lon_box, lat_box, 'color','red','linewidth',1.3);

    %luzon上升流系统区域
    hold on
    % ===== 红色区域框 =====
    lon_box = [117   121    121 117  117];
    lat_box = [16.5  16.5   21  21  16.5];
    m_line(lon_box, lat_box, 'color','k','linewidth',1.3);

    %  标题
    title(title_trend{i},'fontsize',16,'fontname','Times New Roman');
    
    % ===== 显著性 =====
    sig = (p_here' == 1);
    m_plot(Lon2(sig),Lat2(sig),'.k','markersize',1.5);
    
    hold off
end



%% 画南海mhws时间序列趋势分析

% ===== 区域范围 =====
x = find(Lon>=105 & Lon<=122);
y = find(Lat>=1  &  Lat<=23);

years = (1982:2024)';

% ===== 所有变量 =====
vars = {
    annual_Frequency, 'Frequency (Count)',  '(a)';
    annual_Days,      'Total Days (Days)',  '(b)';
    annual_Duration,  'Duration (Days)',    '(c)';
    annual_CumInt,    'CumInt (℃ days)',   '(d)';
    annual_MaxInt,    'MaxInt (℃)',        '(e)';
    annual_MeanInt,   'MeanInt (℃)',       '(f)';
    };

% ===== bootstrap次数 =====
nboot = 1000;

figure('color','w','position',[100 100 1600 900]);
tiledlayout(2,3,'TileSpacing','tight','Padding','loose');

for k = 1:6
    
    data = vars{k,1};
    ylabel_str = vars{k,2};
    label_str = vars{k,3};
    
    % ===== 区域平均 =====
    tmp = data(x,y,:);
    ts = squeeze(nanmean(tmp, [1 2]));
    
    % ===== 线性拟合 =====
    p = polyfit(years, ts, 1);
    y_fit = polyval(p, years);
    
    % ===== MK检验 =====
    result = mk_test(ts);
    trend = result.slope_decade;
    p_value = result.p;
    
    % ===== Bootstrap 计算趋势误差 =====
    slopes = zeros(nboot,1);
    n = length(ts);
    
    for b = 1:nboot
        idx = randi(n,n,1);   % 重采样
        ts_boot = ts(idx);
        years_boot = years(idx);
        
        p_boot = polyfit(years_boot, ts_boot, 1);
        slopes(b) = p_boot(1)*10;  % decade
    end
    
    trend_std = std(slopes);   % 不确定性
    
    % ===== 子图 =====
    nexttile; hold on
    
    plot(years, ts, 'k-o','LineWidth',1.2,'MarkerSize',4)
    plot(years, y_fit, 'r--','LineWidth',1.5)
    
    grid on; box on
    title([label_str ' ' ylabel_str],'FontWeight','normal','FontSize',18,'FontName','times')
    
    % ===== 标注 =====
    y_max = max(ts);
    y_min = min(ts);
    
    % trend ± std
    text(years(2), y_max - 0.01*(y_max-y_min), ...
        sprintf('Trend = %.2f \\pm %.2f /decade', trend, trend_std), ...
        'fontsize',16,'color','r','FontName','times');
    
    % p值规则
    if p_value < 0.01
        p_str = 'p < 0.01';
    elseif p_value < 0.05
        p_str = 'p < 0.05';
    else
        p_str = sprintf('p = %.3f', p_value);
    end
    
    text(years(2), y_max - 0.2*(y_max-y_min), ...
        p_str, 'fontsize',16,'FontName','times');
    
    if k > 3
        xlabel('Year')
    end
    
    ylabel(ylabel_str)
    
end

%% 提取南海年际平均mhw指标
frequency = squeeze(nanmean(annual_Frequency, [1 2]));
duration= squeeze(nanmean(annual_Duration, [1 2]));
days= squeeze(nanmean(annual_Days, [1 2]));
cumint= squeeze(nanmean(annual_CumInt, [1 2]));
maxint= squeeze(nanmean(annual_MaxInt, [1 2]));
meanint= squeeze(nanmean(annual_MeanInt, [1 2]));
