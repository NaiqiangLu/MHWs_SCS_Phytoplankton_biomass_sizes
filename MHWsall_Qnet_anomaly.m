%% Qnet 净海气热通量在整个mhws期间的变化异常

clear;clc;close all;

load Qnet_interp025_1998_2024.mat %加载Qnet数据
Qnet = permute(Qnet,[2 1 3]);  %转置数据 交换维度，把叶绿素矩阵正常投影

load Clim_Qnet_SSR_STR_SLHF_SSHF.mat    %加载Qnet每日气候态数据
Qnet_clim=permute(Qnet_clim,[2,1,3]);

studytime=datenum(1998,1,1):datenum(2024,12,31);    %csd的研究时间为2003-2024年
period_plot_v=datevec(studytime);    %对齐 daily climatology
period_unique=datevec(datenum(2016,1,1):datenum(2016,12,31));    %构建如闰年DOY序列
[~,loc_plot]=ismember(period_plot_v(:,2:3),period_unique(:,2:3),'rows');   %取研究时间的月日然后把他投射到闰年doy序列
mQnet = squeeze(Qnet_clim(:,:,loc_plot));    %1998–2024 每一天对应的 climatology Qnet,同时用于去掉长度为 1 的维度。

%% SST 读取海温数据和热浪
load sst_19822025.mat
tic
[MHW,smclim,m90,mhw_ts]=detect(sst_full,datenum(1982,1,1):datenum(2025,12,31),...
    datenum(1995,1,1),datenum(2024,12,31),...
    datenum(1998,1,1),datenum(2024,12,31),'Threshold',0.9);
toc

%% 所有MHW事件下 Qnet anomaly（不分季节）
tic

mhw = MHW{:,:};   % 转为矩阵

% ===== 输出变量初始化 =====
VQnet = nan(size(Qnet,1),size(Qnet,2));
SigMask_val   = nan(size(Qnet,1),size(Qnet,2));
% ===== 找所有发生过MHW的网格 =====
loc_full = unique(mhw(:,8:9),'rows');

for m = 1:size(loc_full,1)

    loc_here = loc_full(m,:);   % 当前网格点

    % ===== 当前网格所有MHW事件 =====
    mhw_here = mhw(mhw(:,8)==loc_here(1) & mhw(:,9)==loc_here(2),:);

    % ===== 每次MHW的时间段 =====
    period_mhw = [datenum(num2str(mhw_here(:,1)),'yyyymmdd'), ...
                  datenum(num2str(mhw_here(:,2)),'yyyymmdd')];

    vqnet = nan(size(period_mhw,1),1);
      %一定要初始化矩阵，不然数据会被覆盖
    for loc = 1:size(period_mhw,1)

        % ===== 当前事件时间 =====
        mhw_time = period_mhw(loc,1):period_mhw(loc,2);

        % ===== 提取数据 =====
        mhw_qnet = squeeze(Qnet(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));

        clim_qnet = squeeze(mQnet(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));
        
        % ===== 单次事件 anomaly =====
        vqnet(loc) = nanmean(mhw_qnet - clim_qnet);  % 异常值

    end

    % ===== 多事件平均 =====   
    VQnet(loc_here(1),loc_here(2)) =nanmean(vqnet);

    % 1===== 显著性检验 value
    if sum(~isnan(vqnet)) > 1
        [h,p] = ttest(vqnet,0,'Alpha',0.05);

        if h==1
            SigMask_val(loc_here(1),loc_here(2)) = 1;
        else
            SigMask_val(loc_here(1),loc_here(2)) = 0;
        end
    end

end
toc


%% ===== Figure: Qnet anomaly =====
tic
figure('Position',[200,200,900,800],'Color','w');

% ===== 投影 =====
m_proj('equidistant','lon',[105 123],'lat',[0 25]);

% ===== 填色图 =====
m_contourf(Lon,Lat,VQnet',80,'linestyle','none');
shading flat      % 建议用flat，方便保持色阶分明

hold on

% ===== 添加等值线（无数字）=====
contour_levels = -25:25:50;   % 等值线间隔可调整

m_contour(Lon,Lat,VQnet',contour_levels,...
    'LineColor','k',...
    'LineWidth',0.8);

% ===== 显著性黑点 =====
[i,j]=find(SigMask_val==1);
m_scatter(Lon(i),Lat(j),3,'k','filled');

% ===== 海岸线 =====
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[0.7 0.7 0.7]);

% ===== 网格 =====
m_grid('linestyle','none',...
       'linewidth',1.5,...
       'fontsize',20,...
       'xtick',107:5:123,...
       'ytick',0:5:25,...
       'fontname','times');

% ===== 14级离散色带 =====
colormap(nclCM(130,14));

caxis([-25  55]);

% ===== Colorbar =====
cb=colorbar('southoutside');

set(cb,'linewidth',1.5,...
       'fontsize',16,...
       'edgecolor','k');

ylabel(cb,'Qnet (W m^{-2})',...
       'FontSize',18,...
       'FontName','Times New Roman');

% ===== colorbar分割黑线 =====
drawnow

pos = cb.Position;

for k = 2:14
    x = pos(1)+(k-1)/14*pos(3);

    annotation('line',...
        [x x],...
        [pos(2) pos(2)+pos(4)],...
        'Color','k',...
        'LineWidth',0.8);
end


