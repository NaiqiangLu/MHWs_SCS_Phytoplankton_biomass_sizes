%% MLD混合层深度在mhw期间的变化异常

clear;clc;close all;

load mld1998_2024_0.25.mat   %加载CSD数据

MLD = permute(MLD,[2 1 3]); %转置数据 交换维度，把叶绿素矩阵正常投影

load clim_MLD1998_2024.mat    %加载ibar每日气候态数据；clim_CSD为变量名
MLDclim=permute(MLDclim,[2,1,3]);

studytime=datenum(1998,1,1):datenum(2024,12,31);    %csd的研究时间为2003-2024年
period_plot_v=datevec(studytime);    %对齐 daily climatology
period_unique=datevec(datenum(2016,1,1):datenum(2016,12,31));    %构建如闰年DOY序列
[~,loc_plot]=ismember(period_plot_v(:,2:3),period_unique(:,2:3),'rows');   %取研究时间的月日然后把他投射到闰年doy序列
mMLD = squeeze(MLDclim(:,:,loc_plot));    %2003–2024 每一天对应的 climatology Chla,同时用于去掉长度为 1 的维度。

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
VMLD = nan(size(MLD,1),size(MLD,2));
RMLD = nan(size(MLD,1),size(MLD,2));
SigMask_val   = nan(size(MLD,1),size(MLD,2));
SigMask_per   = nan(size(MLD,1),size(MLD,2));

% ===== 找所有发生过MHW的网格 =====
loc_full = unique(mhw(:,8:9),'rows');

for m = 1:size(loc_full,1)

    loc_here = loc_full(m,:);   % 当前网格点

    % ===== 当前网格所有MHW事件 =====
    mhw_here = mhw(mhw(:,8)==loc_here(1) & mhw(:,9)==loc_here(2),:);

    % ===== 每次MHW的时间段 =====
    period_mhw = [datenum(num2str(mhw_here(:,1)),'yyyymmdd'), ...
                  datenum(num2str(mhw_here(:,2)),'yyyymmdd')];

    vmld = nan(size(period_mhw,1),1);
    rmld = nan(size(period_mhw,1),1);  %一定要初始化矩阵，不然数据会被覆盖
    for loc = 1:size(period_mhw,1)

        % ===== 当前事件时间 =====
        mhw_time = period_mhw(loc,1):period_mhw(loc,2);

        % ===== 提取数据 =====
        mhw_mld = squeeze(MLD(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));

        clim_mld = squeeze(mMLD(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));

        % ===== 单次事件 anomaly =====
        vmld(loc) = nanmean(mhw_mld - clim_mld);  % 异常值

        rmld(loc) = nanmean(((mhw_mld - clim_mld)./clim_mld).*100); %异常百分比值
    end

    % ===== 多事件平均 =====
    
    VMLD(loc_here(1),loc_here(2)) =nanmean(vmld);

    RMLD(loc_here(1),loc_here(2))= nanmean(rmld);


    % 1===== 显著性检验 value
    if sum(~isnan(vmld)) > 1
        [h,p] = ttest(vmld,0,'Alpha',0.05);

        if h==1
            SigMask_val(loc_here(1),loc_here(2)) = 1;
        else
            SigMask_val(loc_here(1),loc_here(2)) = 0;
        end
    end


    % 2===== 显著性检验 percent
    if sum(~isnan(rmld)) > 1
        [h,p] = ttest(rmld,0,'Alpha',0.05);

        if h==1
            SigMask_per(loc_here(1),loc_here(2)) = 1;
        else
            SigMask_per(loc_here(1),loc_here(2)) = 0;
        end
    end

end
toc


%% ===== Figure: all MHW Mixed layer-averaged irradiance anomaly 百分比=====
tic
figure('Position',[200,200,900,800]);
% ===== 投影 =====
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
% ===== 绘制ibar异常 =====
m_contourf(Lon, Lat, RMLD', 500, 'linestyle', 'none'); % 注意VCSDslope转置
shading interp
hold on;
% ===== 显著性黑点 =====
[i, j] = find(SigMask_per==1);
m_scatter(Lon(i), Lat(j), 3, 'k', 'filled');  %调节点大小
% ===== 海岸线 =====
m_gshhs_i('linewidth', 1.5, 'color', 'k');
m_gshhs_i('patch', [.7 .7 .7]);  % 灰色陆地填充
% ===== 网格线 =====
m_grid('linestyle','none','linewidth',1,'fontsize',16,...
    'xtick',106:5:123,'ytick',0:5:25);

% ===== 配色方案 =====
colormap(m_colmap('diverging',20));   % 你原来的色带
caxis([-40,  40]);        % 设置颜色范围

% ===== Colorbar =====
cb = colorbar("southoutside");
set(cb,'YTick',-50:10:50);  
set(cb,'linewidth',1,'fontsize',12,'edgecolor','k');
set(get(cb,'ylabel'),'string','Mixed layer depth anomaly(%)','FontSize',18,'FontName','Times New Roman');
% title(cb,'(%)','FontSize',16,'FontName','Times New Roman');
title('Mixed layer-averaged irradiance during MHWs.','FontSize',12,'FontWeight','bold');
toc

%% ===== Figure: all Mixed layer-averaged irradiance anomaly 异常值=====
tic
figure('Position',[200,200,900,800]);
% ===== 投影 =====
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
% ===== 绘制ibar异常 =====
m_contourf(Lon, Lat, VMLD', 500, 'linestyle', 'none'); % 注意VCSDslope转置
shading interp
hold on;
% ===== 显著性黑点 =====
[i, j] = find(SigMask_val==1);
m_scatter(Lon(i), Lat(j), 3, 'k', 'filled');  %调节点大小
% ===== 海岸线 =====
m_gshhs_i('linewidth', 1.5, 'color', 'k');
m_gshhs_i('patch', [.7 .7 .7]);  % 灰色陆地填充
% ===== 网格线 =====
m_grid('linestyle','none','linewidth',1,'fontsize',16,...
    'xtick',106:5:123,'ytick',0:5:25);

% ===== 配色方案 =====
colormap(m_colmap('diverging',20));   % 你原来的色带
caxis([-10,  10]);            % 设置颜色范围

% ===== Colorbar =====
cb = colorbar("southoutside");
set(cb,'YTick',-15:5:15);  
set(cb,'linewidth',1,'fontsize',16,'edgecolor','k');
set(get(cb,'ylabel'),'string','Mixed layer depth anomaly(m)','FontSize',18,'FontName','Times New Roman');
% title(cb,'(%)','FontSize',18,'FontName','Times New Roman'); %给colorbar加%号
title('Mixed layer-averaged irradiance during MHWs.',12,'FontWeight','bold');
toc




