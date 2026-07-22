%% ssw 海表风场在mhw期间的变化异常

clear;clc;close all;

load CCMP_wind1998_2024.mat   %加载数据
load climWind_1998_2024.mat    %加载wind每日气候态数据；

studytime=datenum(1998,1,1):datenum(2024,12,31);    %csd的研究时间为2003-2024年
period_plot_v=datevec(studytime);    %对齐 daily climatology
period_unique=datevec(datenum(2016,1,1):datenum(2016,12,31));    %构建如闰年DOY序列
[~,loc_plot]=ismember(period_plot_v(:,2:3),period_unique(:,2:3),'rows');   %取研究时间的月日然后把他投射到闰年doy序列

mWindS = squeeze(ws_clim(:,:,loc_plot));    %1998–2024 每一天对应的 climatology Chla,同时用于去掉长度为 1 的维度。
mU = squeeze(u_clim(:,:,loc_plot));
mV = squeeze(v_clim(:,:,loc_plot));

%% SST 读取海温数据和热浪
load sst_19822025.mat
tic
[MHW,smclim,m90,mhw_ts]=detect(sst_full,datenum(1982,1,1):datenum(2025,12,31),...
    datenum(1995,1,1),datenum(2024,12,31),...
    datenum(1998,1,1),datenum(2024,12,31),'Threshold',0.9);
toc

%% 所有MHW事件下 CSD anomaly（不分季节）
tic

mhw = MHW{:,:};   % 转为矩阵

% ===== 输出变量初始化 =====
ValWindS = nan(size(ws,1),size(ws,2));
RWindS = nan(size(ws,1),size(ws,2));
UWinds = nan(size(ws,1),size(ws,2));
VWinds = nan(size(ws,1),size(ws,2));
SigMask_val   = nan(size(ws,1),size(ws,2));
SigMask_per   = nan(size(ws,1),size(ws,2));

% ===== 找所有发生过MHW的网格 =====
loc_full = unique(mhw(:,8:9),'rows');  %筛选唯一发生mhw的网格，提取一次就行，某个网格可能会发生多次mhw

for m = 1:size(loc_full,1)

    loc_here = loc_full(m,:);   % 索引当前的mhw的网格点

    % ===== 当前网格所有MHW事件 =====
    mhw_here = mhw(mhw(:,8)==loc_here(1) & mhw(:,9)==loc_here(2),:);  %把某个网格的所有mhw事件筛选出来

    % ===== 每次MHW的时间段 =====
    period_mhw = [datenum(num2str(mhw_here(:,1)),'yyyymmdd'), ...
                  datenum(num2str(mhw_here(:,2)),'yyyymmdd')];  %索引某个网格所有mhw事件的起止时间，某个网格的所有mhw事件的每次起止事件被提取出来

    valwinds = nan(size(period_mhw,1),1);  %为某个网格的所有mhw事件xx值留下赋值矩阵
    rwinds = nan(size(period_mhw,1),1);  %一定要初始化矩阵，不然数据会被覆盖
    uwinds=nan(size(period_mhw,1),1);  %初始化空矩阵
    vwinds=nan(size(period_mhw,1),1);  %初始化空矩阵
    
    for loc = 1:size(period_mhw,1)   %索引某网格mhw事件行，表示按顺序一个一个mhw来进行计算

        % ===== 当前事件时间 =====
        mhw_time = period_mhw(loc,1):period_mhw(loc,2);  %提取每次mhw事件的时间

        % ===== 提取数据 =====
        mhw_ws = squeeze(ws(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1)); %把某个发生mhw的网格的数据提取出来，mhw_time - datenum(1998,1,1) + 1)的目的是精准定位发生mhw是在（100*100*9862天）的第几天到第几天

        clim_ws = squeeze(mWindS(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1)); %同理，与mhw事件对应的网格被提取出来，然后把相应的mhw期间对应所在气候态日期的映射，因为366日气候态数据已经事先把 月日 → 映射成闰年的 DOY 编号
        %因为数据是 3 维矩阵 (100×100×9862)，当固定了前两个经纬度网格后，数据会变成 1×1×N 的三维矩阵！用 squeeze 把两个长度为 1 的维度删掉，变成 1 维向量 N，才能用
        % ===== 去异常值 =====
        mhw_ws(mhw_ws< 0)  = NaN;
        clim_ws(clim_ws< 0)= NaN;


        mhw_u = squeeze(u(loc_here(1),loc_here(2),...
            mhw_time-datenum(1998,1,1)+1));

        clim_u = squeeze(mU(loc_here(1),loc_here(2),...
            mhw_time-datenum(1998,1,1)+1));

        mhw_v = squeeze(v(loc_here(1),loc_here(2),...
            mhw_time-datenum(1998,1,1)+1));

        clim_v = squeeze(mV(loc_here(1),loc_here(2),...
            mhw_time-datenum(1998,1,1)+1));

        uwinds(loc) = nanmean(mhw_u - clim_u);
        vwinds(loc) = nanmean(mhw_v - clim_v);  %用某次mhw事件期间的uv分量减去气候态uv分量 uv分量正负号表示方向，数值表示大小，是矢量，u为纬向（东西风，东为正），v是经向(南北风，北为正），这句表示风向异常
        % uwinds(loc) = nanmean(mhw_u);
        % vwinds(loc) = nanmean(mhw_v);   %不减去气候态 uv分量，表示求所有mhw期间的平均风向

        % ===== 单次事件 anomaly =====
        valwinds(loc) = nanmean(mhw_ws - clim_ws);  %把每次mhw事件的xx值减去相应气候态doy的值，得到mhw期间某指标相对应与气候态的异常值

        rwinds(loc) = nanmean(((mhw_ws - clim_ws)./clim_ws).*100);
    end

    % ===== 多事件平均 =====

    ValWindS(loc_here(1),loc_here(2)) =nanmean(valwinds);

    RWindS(loc_here(1),loc_here(2))= nanmean(rwinds);

    UWinds(loc_here(1),loc_here(2)) = nanmean(uwinds);
    VWinds(loc_here(1),loc_here(2)) = nanmean(vwinds);   %对当前格点把所有mhw事件求得的异常值索引到当前计算的网格并赋值进去

    % 1===== 显著性检验 value
    if sum(~isnan(valwinds)) > 1
        [h,p] = ttest(valwinds,0,'Alpha',0.05);

        if h==1
            SigMask_val(loc_here(1),loc_here(2)) = 1;
        else
            SigMask_val(loc_here(1),loc_here(2)) = 0;
        end
    end


    % 2===== 显著性检验 percent
    if sum(~isnan(rwinds)) > 1
        [h,p] = ttest(rwinds,0,'Alpha',0.05);

        if h==1
            SigMask_per(loc_here(1),loc_here(2)) = 1;
        else
            SigMask_per(loc_here(1),loc_here(2)) = 0;
        end
    end

end
toc


%% ===== Figure: all MHW winds anomaly 百分比=====
[Lon2, Lat2] = meshgrid(Lon, Lat);
figure('Position',[200,200,1000,900],'Resize','off','Color','w');
% ===== 投影 =====
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
% ===== 绘制winds异常 =====
m_contourf(Lon, Lat, RWindS', 50, 'linestyle', 'none'); % 注意VCSDslope转置
shading interp

hold on
skip = 4;  % 控制箭头稀疏程度
m_quiver(Lon2(1:skip:end,1:skip:end), Lat2(1:skip:end,1:skip:end), ...
    UWinds(1:skip:end,1:skip:end)', ...
    VWinds(1:skip:end,1:skip:end)', ...
    'k');

hold on
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',16,...
    'xtick',[106:5:123],'ytick',[0:5:25]);

% ===== 配色方案 =====
colormap(nclCM(156,20));
caxis([-40,  40]); %上下限
cb=colorbar('southoutside');
set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k');
set(get(cb,'ylabel'),'string','Wind Speed anomaly(%)','FontSize',20,'FontName','Times New Roman');%生成colorbar标题，单位
title('Wind Speed anomaly during MHWs','FontSize',16,'FontWeight','bold');


%% ===== Figure: all mhws wind speed anomaly 异常值=====

[Lon2, Lat2] = meshgrid(Lon, Lat);
tic
figure('Position',[200,200,1000,900],'Resize','on','Color','w');
% ===== 投影 =====
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
% ===== 绘制winds异常 =====
m_contourf(Lon, Lat, ValWindS', 50, 'linestyle', 'none'); % 注意VCSDslope转置
shading interp
hold on
skip = 4;  % 控制箭头稀疏程度
m_quiver(Lon2(1:skip:end,1:skip:end), Lat2(1:skip:end,1:skip:end), ...
    UWinds(1:skip:end,1:skip:end)', ...
    VWinds(1:skip:end,1:skip:end)', ...
    'k');

hold on;
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',16,...
    'xtick',[106:5:123],'ytick',[0:5:25]);

% ===== 配色方案 =====
colormap(nclCM(156,20));
caxis([-3.5,  3.5]); %上下限
cb=colorbar('southoutside');
set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k');
set(get(cb,'ylabel'),'string','Wind Speed anomaly(m/s)','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位
title('Wind Speed anomaly during MHWs','FontSize',16,'FontWeight','bold');





