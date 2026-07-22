%% UI上升流强度 在mhw期间的变化异常

clear;clc;close all;

load UI&curl1998-2024.mat   %加载UI数据
load UI_Curl_clim1998_2024.mat     %加载ui_vn的每日气候态数据 clim_ui_vn

load CCMP_wind1998_2024.mat   %加载数据
load climWind_1998_2024.mat    %加载wind每日气候态数据；

studytime=datenum(1998,1,1):datenum(2024,12,31);    %csd的研究时间为2003-2024年
period_plot_v=datevec(studytime);    %对齐 daily climatology
period_unique=datevec(datenum(2016,1,1):datenum(2016,12,31));    %构建如闰年DOY序列
[~,loc_plot]=ismember(period_plot_v(:,2:3),period_unique(:,2:3),'rows');   %取研究时间的月日然后把他投射到闰年doy序列

mUI_VN = squeeze(UI_VNclim(:,:,loc_plot));    %2003–2024 每一天对应的 climatology Chla,同时用于去掉长度为 1 的维度。
mUI_LZ =squeeze(UI_LZclim(:,:,loc_plot));
mU = squeeze(u_clim(:,:,loc_plot));
mV = squeeze(v_clim(:,:,loc_plot));
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
VUI_VN = nan(size(UI_VN,1),size(UI_VN,2));
RUI_VN = nan(size(UI_VN,1),size(UI_VN,2));
VUI_LZ = nan(size(UI_LZ,1),size(UI_LZ,2));
UWinds = nan(size(ws,1),size(ws,2));
VWinds = nan(size(ws,1),size(ws,2));
% SigMaskVN_val   = nan(size(UI_VN,1),size(UI_VN,2));
% SigMaskVN_per   = nan(size(UI_VN,1),size(UI_VN,2));
% SigMaskLZ_val   = nan(size(UI_LZ,1),size(UI_LZ,2));

% ===== 找所有发生过MHW的网格 =====
loc_full = unique(mhw(:,8:9),'rows');

for m = 1:size(loc_full,1)

    loc_here = loc_full(m,:);   % 索引当前mhws网格点

    % ===== 当前网格所有MHW事件 =====
    mhw_here = mhw(mhw(:,8)==loc_here(1) & mhw(:,9)==loc_here(2),:);

    % ===== 每次MHW的时间段 =====
    period_mhw = [datenum(num2str(mhw_here(:,1)),'yyyymmdd'), ...
                  datenum(num2str(mhw_here(:,2)),'yyyymmdd')];

    vui_vn = nan(size(period_mhw,1),1);
    rui_vn = nan(size(period_mhw,1),1);  %一定要初始化矩阵，不然数据会被覆盖
    vui_lz = nan(size(period_mhw,1),1);
    uwinds=nan(size(period_mhw,1),1);  %初始化空矩阵
    vwinds=nan(size(period_mhw,1),1);  %初始化空矩阵
    for loc = 1:size(period_mhw,1)

        % ===== 当前事件时间 =====
        mhw_time = period_mhw(loc,1):period_mhw(loc,2);

        % ===== 提取数据 =====
        mhw_ui_vn = squeeze(UI_VN(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));

        clim_ui_vn = squeeze(mUI_VN(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));

        mhw_ui_lz = squeeze(UI_LZ(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));

        clim_ui_lz = squeeze(mUI_LZ(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));


        mhw_u = squeeze(u(loc_here(1),loc_here(2),...
            mhw_time-datenum(1998,1,1)+1));

        clim_u = squeeze(mU(loc_here(1),loc_here(2),...
            mhw_time-datenum(1998,1,1)+1));

        mhw_v = squeeze(v(loc_here(1),loc_here(2),...
            mhw_time-datenum(1998,1,1)+1));

        clim_v = squeeze(mV(loc_here(1),loc_here(2),...
            mhw_time-datenum(1998,1,1)+1));

        uwinds(loc) = nanmean(mhw_u - clim_u);
        vwinds(loc) = nanmean(mhw_v - clim_v);  %这两句表示风向异常


        % ===== 单次MHW事件 anomaly =====
        vui_vn(loc) = nanmean(mhw_ui_vn - clim_ui_vn);   %单次mhw事件的异常值

        rui_vn(loc) = nanmean(((mhw_ui_vn - clim_ui_vn)./clim_ui_vn).*100);  %单次mhw事件的异常程度（%）

        vui_lz(loc) = nanmean(mhw_ui_lz - clim_ui_lz);   %单次mhw事件的异常值
    end

    % ===== 多事件平均 =====

    VUI_VN(loc_here(1),loc_here(2)) =nanmean(vui_vn);
    RUI_VN(loc_here(1),loc_here(2))= nanmean(rui_vn);

    VUI_LZ(loc_here(1),loc_here(2)) =nanmean(vui_lz);

    UWinds(loc_here(1),loc_here(2)) = nanmean(uwinds);
    VWinds(loc_here(1),loc_here(2)) = nanmean(vwinds);

    % % 1===== 显著性检验 value
    % if sum(~isnan(vui_vn)) > 1
    %     [h,p] = ttest(vui_vn,0,'Alpha',0.05);
    % 
    %     if h==1
    %         SigMaskVN_val(loc_here(1),loc_here(2)) = 1;
    %     else
    %         SigMaskVN_val(loc_here(1),loc_here(2)) = 0;
    %     end
    % end
    % 
    % if sum(~isnan(vui_lz)) > 1
    %     [h,p] = ttest(vui_lz,0,'Alpha',0.05);
    % 
    %     if h==1
    %         SigMaskLZ_val(loc_here(1),loc_here(2)) = 1;
    %     else
    %         SigMaskLZ_val(loc_here(1),loc_here(2)) = 0;
    %     end
    % end
    % 
    % % 2===== 显著性检验 percent
    % if sum(~isnan(rui_vn)) > 1
    %     [h,p] = ttest(rui_vn,0,'Alpha',0.05);
    % 
    %     if h==1
    %         SigMaskVN_per(loc_here(1),loc_here(2)) = 1;
    %     else
    %         SigMaskVN_per(loc_here(1),loc_here(2)) = 0;
    %     end
    % end

end
toc


%% ===== Figure:  ui in vietnam anomaly% during all mhws events 百分比=====
[Lon2, Lat2] = meshgrid(Lon, Lat); %初始化风矢量网格

figure('Position',[200,200,1000,900],'Color','w');
% ===== 投影 =====
m_proj('equidistant','lon',[105 115],'lat',[5 18]);
% ===== 绘制ibar异常 =====
m_contourf(Lon, Lat, RUI_VN', 50, 'linestyle', 'none'); % 注意VCSDslope转置
shading interp

% hold on;
% % ===== 显著性黑点 =====
% [i, j] = find(SigMaskVN_per==1);
% m_scatter(Lon(i), Lat(j), 2.5, 'k', 'filled');  %调节点大小

% ===== 海岸线 =====
hold on
m_gshhs_i('linewidth', 1.5, 'color', 'k');
m_gshhs_i('patch', [.7 .7 .7]);  % 灰色陆地填充
% ===== 网格线 =====
m_grid('linestyle','none','linewidth',1,'fontsize',16,...
    'xtick',106:3:123,'ytick',0:3:25);

% ===== 配色方案 =====
colormap(m_colmap('diverging',20));   % 你原来的色带
caxis([-60,  60]);        % 设置颜色范围

% ===== Colorbar =====
cb = colorbar("southoutside");
set(cb,'linewidth',1,'fontsize',16,'edgecolor','k','YTick',-60:20:60);
set(get(cb,'ylabel'),'string','UI anomaly(%)','FontSize',18,'FontName','Times New Roman');
% title(cb,'(%)','FontSize',14,'FontName','Times New Roman');
title('Upwelling index during MHWs.','fontname','times','FontSize',18,'FontWeight','bold');


%% ===== Figure: ui in vietnam outlier during all mhws events 异常值=====
tic
figure('Position',[200,200,1000,900],'Color','w');
% ===== 投影 =====
m_proj('equidistant','lon',[105 115],'lat',[5 18]);
% ===== 绘制ibar异常 =====
m_contourf(Lon, Lat, VUI_VN', 50, 'linestyle', 'none'); % 注意VCSDslope转置
shading interp

hold on
skip = 2;  % 控制箭头稀疏程度
m_quiver(Lon2(1:skip:end,1:skip:end), Lat2(1:skip:end,1:skip:end), ...
    UWinds(1:skip:end,1:skip:end)', ...
    VWinds(1:skip:end,1:skip:end)', ...
    'k');

% hold on;
% % ===== 显著性黑点 =====
% [i, j] = find(SigMaskVN_val==1);
% m_scatter(Lon(i), Lat(j), 2.5, 'k', 'filled');  %调节点大小

% ===== 海岸线 =====
hold on
m_gshhs_i('linewidth', 1.5, 'color', 'k');
m_gshhs_i('patch', [.7 .7 .7]);  % 灰色陆地填充
% ===== 网格线 =====
m_grid('linestyle','none','linewidth',1,'fontsize',16,...
    'xtick',106:3:123,'ytick',0:3:25);

% ===== 配色方案 =====
colormap(m_colmap('diverging',20));   % 你原来的色带
caxis([-2.5, 2.5]);            % 设置颜色范围

% ===== Colorbar =====
cb = colorbar("southoutside");
set(cb,'linewidth',1,'fontsize',16,'edgecolor','k','YTick',-2.5:1:2.5);
set(get(cb,'ylabel'),'string','UI (m^{2} s^{-1})','FontSize',18,'FontName','Times New Roman');

title('Upwelling index during MHWs.','fontname','times','fontsize',18,'FontWeight','bold');
toc

%% ===== Figure: ui in Luzon outlier during all mhws events 异常值=====
tic
figure('Position',[200,200,1000,900],'Color','w');
% ===== 投影 =====
m_proj('equidistant','lon',[116 123],'lat',[16 24]);
% ===== 绘制UI_LZ异常 =====
m_contourf(Lon, Lat, VUI_LZ', 50, 'linestyle', 'none'); % 注意数据的转置画图
shading interp

hold on
skip = 2;  % 控制箭头稀疏程度
m_quiver(Lon2(1:skip:end,1:skip:end), Lat2(1:skip:end,1:skip:end), ...
    UWinds(1:skip:end,1:skip:end)', ...
    VWinds(1:skip:end,1:skip:end)', ...
    'k');

% hold on;
% % ===== 显著性黑点 =====
% [i, j] = find(SigMaskLZ_val==1);
% m_scatter(Lon(i), Lat(j), 2.5, 'k', 'filled');  %调节点大小

% ===== 海岸线 =====
hold on
m_gshhs_i('linewidth', 1.5, 'color', 'k');
m_gshhs_i('patch', [.7 .7 .7]);  % 灰色陆地填充
% ===== 网格线 =====
m_grid('linestyle','none','linewidth',1,'fontsize',16,...
    'xtick',106:3:123,'ytick',0:3:25);

% ===== 配色方案 =====
colormap(m_colmap('diverging',20));   % 你原来的色带
caxis([-2, 2]);            % 设置颜色范围

% ===== Colorbar =====
cb = colorbar("southoutside");
set(cb,'linewidth',1,'fontsize',16,'edgecolor','k','YTick',-2.5:1:2.5);
set(get(cb,'ylabel'),'string','UI (m^{2} s^{-1})','FontSize',18,'FontName','Times New Roman');

title('Upwelling index during MHWs.','fontname','times','fontsize',18,'FontWeight','bold');
toc


