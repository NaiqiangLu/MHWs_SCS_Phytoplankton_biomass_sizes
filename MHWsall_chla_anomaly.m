%% chla 在全部mhws期间的异常

clear;clc;close all;

load  Chla_19982024_interp025.mat   %加载Chla数据
chla_full(chla_full< 0)=nan;   %变量去异常值
chla_full = permute(chla_full,[2 1 3]); %转置数据 交换维度，把叶绿素矩阵正常投影

load Chla_clim.mat      %加载csd每日气候态数据；clim_CSD为变量名
climchla=permute(climchla,[2,1,3]);

studytime=datenum(1998,1,1):datenum(2024,12,31);    %csd的研究时间为2003-2024年
period_plot_v=datevec(studytime);    %对齐 daily climatology
period_unique=datevec(datenum(2016,1,1):datenum(2016,12,31));    %构建如闰年DOY序列
[~,loc_plot]=ismember(period_plot_v(:,2:3),period_unique(:,2:3),'rows');   %取研究时间的月日然后把他投射到闰年doy序列
mChla = squeeze(climchla(:,:,loc_plot));    %2003–2024 每一天对应的 climatology Chla,同时用于去掉长度为 1 的维度。

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
VChla = nan(size(chla_full,1),size(chla_full,2));
RChla = nan(size(chla_full,1),size(chla_full,2));
SigMask   = nan(size(chla_full,1),size(chla_full,2));

% ===== 找所有发生过MHW的网格 =====
loc_full = unique(mhw(:,8:9),'rows');

for m = 1:size(loc_full,1)

    loc_here = loc_full(m,:);   % 当前网格点

    % ===== 当前网格所有MHW事件 =====
    mhw_here = mhw(mhw(:,8)==loc_here(1) & mhw(:,9)==loc_here(2),:);

    % ===== 每次MHW的时间段 =====
    period_mhw = [datenum(num2str(mhw_here(:,1)),'yyyymmdd'), ...
                  datenum(num2str(mhw_here(:,2)),'yyyymmdd')];

    vchla = nan(size(period_mhw,1),1);
    rchla = nan(size(period_mhw,1),1);  %一定要初始化矩阵，不然数据会被覆盖
    for loc = 1:size(period_mhw,1)

        % ===== 当前事件时间 =====
        mhw_time = period_mhw(loc,1):period_mhw(loc,2);

        % ===== 提取数据 =====
        mhw_chla = squeeze(chla_full(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));

        clim_chla = squeeze(mChla(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));

        % ===== 去异常值 =====
        mhw_chla(mhw_chla< 0)  = NaN;
        clim_chla(clim_chla< 0)= NaN;

        % ===== 单次事件 anomaly =====
        vchla(loc) = nanmean(log10(mhw_chla) - log10(clim_chla));  % log10异常浓度值

        rchla(loc) = nanmean(((mhw_chla - clim_chla)./clim_chla).*100);
    end

    % ===== 多事件平均 =====
    
    VChla(loc_here(1),loc_here(2)) =nanmean(vchla);

    RChla(loc_here(1),loc_here(2))= nanmean(rchla);


    % 1===== 显著性检验 value
    if sum(~isnan(vchla)) > 1
        [h,p] = ttest(vchla,0,'Alpha',0.05);

        if h==1
            SigMask_val(loc_here(1),loc_here(2)) = 1;
        else
            SigMask_val(loc_here(1),loc_here(2)) = 0;
        end
    end


    % 2===== 显著性检验 percent
    if sum(~isnan(rchla)) > 1
        [h,p] = ttest(rchla,0,'Alpha',0.05);

        if h==1
            SigMask_per(loc_here(1),loc_here(2)) = 1;
        else
            SigMask_per(loc_here(1),loc_here(2)) = 0;
        end
    end

end

toc



%% ===== Figure: all MHW fpico anomaly 百分比=====
tic
figure('Position',[200,200,1000,900],'Resize','off','Color','w');
% ===== 投影 =====
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
% ===== 绘制异常 =====
m_contourf(Lon, Lat, RChla', 100, 'linestyle', 'none'); % 注意VCSDslope转置
shading interp
hold on;
% ===== 显著性黑点 =====
[i, j] = find(SigMask_per==1);
m_scatter(Lon(i), Lat(j), 3, 'k', 'filled');  %调节点大小
% ===== 海岸线 =====
m_gshhs_i('linewidth', 1.5, 'color', 'k');
m_gshhs_i('patch', [.7 .7 .7]);  % 灰色陆地填充
% ===== 网格线 =====
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');

% ===== 配色方案 =====
colormap(mycols(3));  % 你原来的色带
caxis([-50,  50]);        % 设置颜色范围

% ===== Colorbar =====
% cb = colorbar;
% set(cb,'YTick',-50:10:50);  
% set(cb,'location','eastoutside','linewidth',1.5,'fontsize',16,'edgecolor','k');
% set(get(cb,'ylabel'),'string','Chl-a anomaly','FontSize',18,'FontName','Times New Roman');
% title(cb,'(%)','FontSize',18,'FontName','Times New Roman');
% title('Chl-a Anomaly during MHWs','FontSize',18,'FontWeight','bold');
toc

%% ===== Figure: all MHW fpico anomaly log10异常值=====
tic
figure('Position',[200,200,1000,900],'Resize','off','Color','w');
% ===== 投影 =====
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
% ===== 绘制异常 =====
m_contourf(Lon, Lat, VChla', 80, 'linestyle', 'none'); % 注意VCSDslope转置
shading interp
hold on;
% ===== 显著性黑点 =====
[i, j] = find(SigMask_val==1);
m_scatter(Lon(i), Lat(j), 3, 'k', 'filled');  %调节点大小
% ===== 海岸线 =====
m_gshhs_i('linewidth', 1.5, 'color', 'k');
m_gshhs_i('patch', [.7 .7 .7]);  % 灰色陆地填充
% ===== 网格线 =====
m_grid('linestyle','none','linewidth',1.5,'fontsize',16,...
    'xtick',106:5:123,'ytick',0:5:25,'fontname','times');

% ===== 配色方案 =====
colormap(m_colmap('diverging',20));   % 你原来的色带
caxis([-0.3,  0.3]);            % 设置颜色范围

% ===== Colorbar =====
cb = colorbar;
set(cb,'YTick',-0.3:0.1:0.3);  
set(cb,'location','eastoutside','linewidth',1.5,'fontsize',16,'edgecolor','k');
set(get(cb,'ylabel'),'string','log_{10}mg m^{-3}','FontSize',28,'FontName','Times New Roman');
% title(cb,'(%)','FontSize',18,'FontName','Times New Roman');
title('Chl-a Anomaly during MHWs','FontSize',18,'FontWeight','bold');
toc




