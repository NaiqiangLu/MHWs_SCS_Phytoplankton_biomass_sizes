%% cellpigmentation 细胞着色度的季节性异常分布

clear;clc;close all;

load  cellpigment_025.mat     %加载数据
cellpig = permute(cellpig,[2 1 3]); %转置数据 交换维度，把叶绿素矩阵正常投影

load  climCellpig.mat   %加载csd每日气候态数据；clim_CSD为变量名
cellpig_clim=permute(cellpig_clim,[2,1,3]);

studytime=datenum(1998,1,1):datenum(2024,12,31);  %csd的研究时间为2003-2024年
period_plot_v=datevec(studytime);   %对齐 daily climatology
period_unique=datevec(datenum(2016,1,1):datenum(2016,12,31));  %构建如闰年DOY序列
[~,loc_plot]=ismember(period_plot_v(:,2:3),period_unique(:,2:3),'rows');  %取研究时间的月日然后把他投射到闰年doy序列
mCellpig = squeeze(cellpig_clim(:,:,loc_plot));  %2003–2024 每一天对应的 climatology Chla,同时用于去掉长度为 1 的维度。

%% SST 读取海温数据和热浪
load sst_19822025.mat
tic
[MHW,smclim,m90,mhw_ts]=detect(sst_full,datenum(1982,1,1):datenum(2025,12,31),...
    datenum(1995,1,1),datenum(2024,12,31),...
    datenum(1998,1,1),datenum(2024,12,31),'Threshold',0.9);
toc

%% 所有MHW事件下 Cellpigment anomaly（不分季节）
tic

mhw = MHW{:,:};   % 转为矩阵

% ===== 输出变量初始化 =====
VCellpig = nan(size(cellpig,1),size(cellpig,2));
RCellpig = nan(size(cellpig,1),size(cellpig,2));
SigMask = nan(size(cellpig,1),size(cellpig,2));

% ===== 找所有发生过MHW的网格 =====
loc_full = unique(mhw(:,8:9),'rows');

for m = 1:size(loc_full,1)

    loc_here = loc_full(m,:);   % 当前网格点

    % ===== 当前网格所有MHW事件 =====
    mhw_here = mhw(mhw(:,8)==loc_here(1) & mhw(:,9)==loc_here(2),:);

    % ===== 每次MHW的时间段 =====
    period_mhw = [datenum(num2str(mhw_here(:,1)),'yyyymmdd'), ...
                  datenum(num2str(mhw_here(:,2)),'yyyymmdd')];

    vcellpig = nan(size(period_mhw,1),1);
    rcellpig = nan(size(period_mhw,1),1);  %一定要初始化矩阵，不然数据会被覆盖
    for loc = 1:size(period_mhw,1)

        % ===== 当前事件时间 =====
        mhw_time = period_mhw(loc,1):period_mhw(loc,2);

        % ===== 提取数据 =====
        mhw_cellpig = squeeze(cellpig(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));

        clim_cellpig = squeeze(mCellpig(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));

        % ===== 单次事件 anomaly =====
        vcellpig(loc) = nanmean(mhw_cellpig - clim_cellpig);
        rcellpig(loc) = nanmean(((mhw_cellpig - clim_cellpig)./clim_cellpig).*100);
    end

    % ===== 多事件平均 =====
    
    VCellpig(loc_here(1),loc_here(2)) =nanmean(vcellpig);

    RCellpig(loc_here(1),loc_here(2))= nanmean(rcellpig);


    % ===== 显著性检验 =====
    if sum(~isnan(vcellpig)) > 1
        [h,p] = ttest(vcellpig,0,'Alpha',0.05);

        if h==1
            SigMask(loc_here(1),loc_here(2)) = 1;
        else
            SigMask(loc_here(1),loc_here(2)) = 0;
        end
    end

end
toc

%% ===== Figure: all MHW fpico anomaly 异常值=====

figure('position',[200,200,1000,900],'Color','w');
% ===== 投影 =====
m_proj('equidistant','lon',[105 123],'lat',[0 25]);

% ===== 绘制异常 =====
m_contourf(Lon, Lat, VCellpig', 100, 'linestyle', 'none'); % 注意VCSDslope转置
shading interp;
hold on;
% ===== 显著性黑点 =====
[i, j] = find(SigMask==1);
m_scatter(Lon(i), Lat(j), 3, 'k', 'filled');  %调节点大小
% ===== 海岸线 =====
m_gshhs_i('linewidth', 1.5, 'color', 'k');
m_gshhs_i('patch', [.7 .7 .7]);  % 灰色陆地填充
% ===== 网格线 =====
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
% ===== 配色方案 =====
colormap(nclCM(160, 18));   %色带
caxis([-25 , 25]);          %设置颜色范围
%===== Colorbar =====
% cb = colorbar('southoutside');
% set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k');
% set(get(cb,'ylabel'),'string','θ anomaly (mgC (mg Chl-a)^{-1})','FontSize',24,'FontName','Times New Roman');
% title(cb,'(%)','FontSize',18,'FontName','Times New Roman');
% title('cellular pigmentation(θ) Anomaly during MHWs','FontSize',18);





