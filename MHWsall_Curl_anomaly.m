%% curl 风应力旋度 在mhw期间的变化异常

clear;clc;close all;

load UI&curl1998-2024.mat   %加载curl数据
load UI_Curl_clim1998_2024.mat     %加载curl的每日气候态数据 Curl_clim

studytime=datenum(1998,1,1):datenum(2024,12,31);    %csd的研究时间为2003-2024年
period_plot_v=datevec(studytime);    %对齐 daily climatology
period_unique=datevec(datenum(2016,1,1):datenum(2016,12,31));    %构建如闰年DOY序列
[~,loc_plot]=ismember(period_plot_v(:,2:3),period_unique(:,2:3),'rows');   %取研究时间的月日然后把他投射到闰年doy序列

mCurl = squeeze(Curl_clim(:,:,loc_plot));    %2003–2024 每一天对应的 climatology Chla,同时用于去掉长度为 1 的维度。

%% SST 读取海温数据和热浪
load sst_19822025.mat
tic
[MHW,smclim,m90,mhw_ts]=detect(sst_full,datenum(1982,1,1):datenum(2025,12,31),...
    datenum(1995,1,1),datenum(2024,12,31),...
    datenum(1998,1,1),datenum(2024,12,31),'Threshold',0.9);
toc

%% 所有MHW事件下 UI anomaly（不分季节）
tic

mhw = MHW{:,:};   % 表格转为矩阵

% ===== 输出变量初始化 =====
VCurl = nan(size(curl,1),size(curl,2));
SigMask_val   = nan(size(curl,1),size(curl,2));

% ===== 找所有发生过MHW的网格 =====
loc_full = unique(mhw(:,8:9),'rows');

for m = 1:size(loc_full,1)

    loc_here = loc_full(m,:);   % 当前网格点

    % ===== 当前网格所有MHW事件 =====
    mhw_here = mhw(mhw(:,8)==loc_here(1) & mhw(:,9)==loc_here(2),:);

    % ===== 每次MHW的时间段 =====
    period_mhw = [datenum(num2str(mhw_here(:,1)),'yyyymmdd'), ...
                  datenum(num2str(mhw_here(:,2)),'yyyymmdd')];

    vcurl = nan(size(period_mhw,1),1);
   
    for loc = 1:size(period_mhw,1)

        % ===== 当前事件时间 =====
        mhw_time = period_mhw(loc,1):period_mhw(loc,2);

        % ===== 提取数据 =====
        mhw_curl = squeeze(curl(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));

        clim_curl = squeeze(mCurl(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));

        % ===== 单次MHW事件 anomaly =====
        vcurl(loc) = nanmean(mhw_curl - clim_curl);   %单次mhw事件的异常值

    end

    % ===== 多事件平均 =====

    VCurl(loc_here(1),loc_here(2)) =nanmean(vcurl);
    
    % 1===== 显著性检验 value
    if sum(~isnan(vcurl)) > 1
        [h,p] = ttest(vcurl,0,'Alpha',0.05);

        if h==1
            SigMask_val(loc_here(1),loc_here(2)) = 1;
        else
            SigMask_val(loc_here(1),loc_here(2)) = 0;
        end
    end

end
toc

VCurl=VCurl./1e-6;  %风应力旋度 N/m^3 （10^-6数量级）正值：气旋式风应力旋度（海洋中引发上升流）负值：反气旋式风应力旋度（海洋中引发下沉流）
%除于数量级方便画图的colorbar控制

%% ===== Figure: curl in vietnam outlier during all mhws events 异常值=====
tic
figure('Position',[200,200,1000,900],'Color','w');
% ===== 投影 =====
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
% ===== 绘制curl异常 =====

m_contourf(Lon, Lat, VCurl', 80, 'linestyle', 'none'); % 注意数据转置画图方向
shading interp

% ===== 显著性黑点 =====
hold on;
[i, j] = find(SigMask_val==1);
m_scatter(Lon(i), Lat(j), 2.5, 'k', 'filled');  %调节点大小

% ===== 海岸线 =====
hold on
m_gshhs_i('linewidth', 1.5, 'color', 'k');
m_gshhs_i('patch', [.7 .7 .7]);  % 灰色陆地填充
% ===== 网格线 =====
m_grid('linestyle','none','linewidth',1,'fontsize',16,...
    'xtick',106:3:123,'ytick',0:3:25);

% ===== 配色方案 =====
colormap(m_colmap('diverging',20));   % 色带颜色
caxis([-0.35, 0.35]);                 % 设置颜色范围

% ===== Colorbar =====
cb = colorbar("southoutside");
set(cb,'linewidth',1,'fontsize',16,'edgecolor','k');
set(get(cb,'ylabel'),'string','Wind stress curl anomalies (10^{-6} N m^{-3})','FontSize',18,'FontName','Times New Roman');

title('Wind stress curl anomalies during MHWs.','fontname','times','fontsize',18,'FontWeight','bold');
toc

%% ===== Figure: curl in Luzon outlier during all mhws events 异常值=====
tic
figure('Position',[200,200,1000,900],'Color','w');
% ===== 投影 =====
m_proj('equidistant','lon',[116 123],'lat',[15 24]);
% ===== 绘制curl异常 =====
m_contourf(Lon, Lat, VCurl', 80, 'linestyle', 'none'); % 注意数据的转置画图
shading interp

% ===== 显著性黑点 =====
hold on;
[i, j] = find(SigMask_val==1);
m_scatter(Lon(i), Lat(j), 2.5, 'k', 'filled');  %调节点大小

% ===== 海岸线 =====
hold on
m_gshhs_i('linewidth', 1.5, 'color', 'k');
m_gshhs_i('patch', [.7 .7 .7]);  % 灰色陆地填充
% ===== 网格线 =====
m_grid('linestyle','none','linewidth',1,'fontsize',16,...
    'xtick',106:3:123,'ytick',0:3:25);

% ===== 配色方案 =====
colormap(m_colmap('diverging',20));   % 色带颜色
caxis([-0.6, 0.6]);                   % 设置颜色范围

% ===== Colorbar =====
cb = colorbar("southoutside");
set(cb,'linewidth',1,'fontsize',16,'edgecolor','k');
set(get(cb,'ylabel'),'string','Wind stress curl anomalies (10^{-6} N m^{-3})','FontSize',18,'FontName','Times New Roman');

title('Wind stress curl anomalies during MHWs.','fontname','times','fontsize',18,'FontWeight','bold');
toc


