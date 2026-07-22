%% rainfall 在整个mhws期间的变化异常

clear;clc;close all;

load Rain_19982024_interp025.mat   %加载降水数据
Rain = permute(Rain,[2 1 3]);  %转置数据 交换维度，让数据矩阵正常投影

load Clim_precipitation.mat    %加载rainafall每日气候态数据
Precip_clim=permute(Precip_clim,[2,1,3]);

studytime=datenum(1998,1,1):datenum(2024,12,31);    %csd的研究时间为2003-2024年
period_plot_v=datevec(studytime);    %对齐 daily climatology
period_unique=datevec(datenum(2016,1,1):datenum(2016,12,31));    %构建如闰年DOY序列
[~,loc_plot]=ismember(period_plot_v(:,2:3),period_unique(:,2:3),'rows');   %取研究时间的月日然后把他投射到闰年doy序列
mPrecip = squeeze(Precip_clim(:,:,loc_plot));    %1998–2024 每一天对应的 climatology rainfall,同时用于去掉长度为 1 的维度。

%% SST 读取海温数据和热浪
load sst_19822025.mat
tic
[MHW,smclim,m90,mhw_ts]=detect(sst_full,datenum(1982,1,1):datenum(2025,12,31),...
    datenum(1995,1,1),datenum(2024,12,31),...
    datenum(1998,1,1),datenum(2024,12,31),'Threshold',0.9);
toc

%% 所有MHW事件下 rainfall anomaly（不分季节）
tic

mhw = MHW{:,:};   % 转为矩阵

% ===== 输出变量初始化 =====
VPrecip = nan(size(Rain,1),size(Rain,2));
RPrecip = nan(size(Rain,1),size(Rain,2));
SigMask_val   = nan(size(Rain,1),size(Rain,2));
SigMask_per   = nan(size(Rain,1),size(Rain,2));
% ===== 找所有发生过MHW的网格 =====
loc_full = unique(mhw(:,8:9),'rows');

for m = 1:size(loc_full,1)

    loc_here = loc_full(m,:);   % 当前网格点

    % ===== 当前网格所有MHW事件 =====
    mhw_here = mhw(mhw(:,8)==loc_here(1) & mhw(:,9)==loc_here(2),:);

    % ===== 每次MHW的时间段 =====
    period_mhw = [datenum(num2str(mhw_here(:,1)),'yyyymmdd'), ...
                  datenum(num2str(mhw_here(:,2)),'yyyymmdd')];

    vprecip = nan(size(period_mhw,1),1);
    rprecip = nan(size(period_mhw,1),1);  %一定要初始化矩阵，不然数据会被覆盖
    for loc = 1:size(period_mhw,1)

        % ===== 当前事件时间 =====
        mhw_time = period_mhw(loc,1):period_mhw(loc,2);

        % ===== 提取数据 =====
        mhw_precip = squeeze(Rain(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));

        clim_precip = squeeze(mPrecip(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));
        
        % ===== 单次事件 anomaly =====
        vprecip(loc) = nanmean(mhw_precip - clim_precip);  % 异常值

        rprecip(loc) = nanmean(((mhw_precip - clim_precip)./clim_precip).*100); %异常百分比值
    end

    % ===== 多事件平均 =====
    
    VPrecip(loc_here(1),loc_here(2)) =nanmean(vprecip);  %异常值

    RPrecip(loc_here(1),loc_here(2))= nanmean(rprecip);  %异常百分比
   
    % 1===== 显著性检验 value
    if sum(~isnan(vprecip)) > 1
        [h,p] = ttest(vprecip,0,'Alpha',0.05);

        if h==1
            SigMask_val(loc_here(1),loc_here(2)) = 1;
        else
            SigMask_val(loc_here(1),loc_here(2)) = 0;
        end
    end

    % 2===== 显著性检验 percent
    if sum(~isnan(rprecip)) > 1
        [h,p] = ttest(rprecip,0,'Alpha',0.05);

        if h==1
            SigMask_per(loc_here(1),loc_here(2)) = 1;
        else
            SigMask_per(loc_here(1),loc_here(2)) = 0;
        end
    end

end
toc


%% ===== Figure: all Rain anomaly 异常值=====
tic
figure('Position',[200,200,900,800]);
% ===== 投影 =====
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
% ===== 绘制rainfall异常 =====
m_contourf(Lon, Lat, VPrecip', 80, 'linestyle', 'none'); % 注意转置
shading interp

hold on;
% ===== 显著性黑点 =====
[i, j] = find(SigMask_val==1);
m_scatter(Lon(i), Lat(j), 3, 'k', 'filled');  %调节点大小
% ===== 海岸线 =====
m_gshhs_i('linewidth', 1.1, 'color', 'k');
m_gshhs_i('patch', [.7 .7 .7]);  % 灰色陆地填充
% ===== 网格线 =====
m_grid('linestyle','none','linewidth',1,'fontsize',16,...
    'xtick',[107:5:123],'ytick',[0:5:25]);

% ===== 配色方案 =====
colormap(m_colmap('diverging',12));   % 你原来的色带
caxis([-3.5,  3.5]);            % 设置颜色范围

% ===== Colorbar =====
cb = colorbar("southoutside");
set(cb,'linewidth',1,'fontsize',16,'edgecolor','k');
set(get(cb,'ylabel'),'string',' Precipitation anomalies (mm)','FontSize',18,'FontName','Times New Roman');
% title(cb,'(%)','FontSize',18,'FontName','Times New Roman'); %给colorbar加%号
title('Precipitation anomalies during MHWs.','fontsize',20,'FontName','times');
toc


% ===== Figure: all Rain anomaly 异常值百分比=====
tic
figure('Position',[200,200,900,800]);
% ===== 投影 =====
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
% ===== 绘制rainfall异常 =====
m_contourf(Lon, Lat, RPrecip', 80, 'linestyle', 'none'); % 注意转置
shading interp

hold on;
% ===== 显著性黑点 =====
[i, j] = find(SigMask_val==1);
m_scatter(Lon(i), Lat(j), 3, 'k', 'filled');  %调节点大小
% ===== 海岸线 =====
m_gshhs_i('linewidth', 1.1, 'color', 'k');
m_gshhs_i('patch', [.7 .7 .7]);  % 灰色陆地填充
% ===== 网格线 =====
m_grid('linestyle','none','linewidth',1,'fontsize',16,...
    'xtick',[106:5:123],'ytick',[0:5:25]);

% ===== 配色方案 =====
colormap(m_colmap('diverging',12));   % 色带
caxis([-50,  50]);            % 设置颜色范围

% ===== Colorbar =====
cb = colorbar("southoutside");
set(cb,'linewidth',1,'fontsize',16,'edgecolor','k');
set(get(cb,'ylabel'),'string',' Precipitation anomalies (mm)','FontSize',18,'FontName','Times New Roman');
% title(cb,'(%)','FontSize',18,'FontName','Times New Roman'); %给colorbar加%号
title('Precipitation anomalies during MHWs.','fontsize',20,'FontName','times');
toc
