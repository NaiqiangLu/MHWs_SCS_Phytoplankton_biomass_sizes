%% 浮游植物粒径比例异常 热浪相对于气候态高出了 多了 xx比例/百分比
clear;clc;close all;

load  Fphyto1998-2024.mat   %加载CSD数据
Fnano(Fnano< 0)=nan;   %变量去异常值
Fnano = permute(Fnano,[2 1 3]); %转置数据 交换维度，把叶绿素矩阵正常投影

load ClimFphyto1998-2024.mat   %加载csd每日气候态数据；clim_CSD为变量名
climFnano=permute(climFnano,[2,1,3]);

studytime=datenum(1998,1,1):datenum(2024,12,31);  %csd的研究时间为2003-2024年
period_plot_v=datevec(studytime);   %对齐 daily climatology
period_unique=datevec(datenum(2016,1,1):datenum(2016,12,31));  %构建如闰年DOY序列
[~,loc_plot]=ismember(period_plot_v(:,2:3),period_unique(:,2:3),'rows');  %取研究时间的月日然后把他投射到闰年doy序列
mFnano = squeeze(climFnano(:,:,loc_plot));  %2003–2024 每一天对应的 climatology Chla,同时用于去掉长度为 1 的维度。

%% SST 读取海温数据和热浪
load sst_19822025.mat
tic
[MHW,smclim,m90,mhw_ts]=detect(sst_full,datenum(1982,1,1):datenum(2025,12,31),...
    datenum(1995,1,1),datenum(2024,12,31),...
    datenum(1998,1,1),datenum(2024,12,31),'Threshold',0.9);
toc

%% 所有MHW事件下 anomaly（不分季节）
tic

mhw = MHW{:,:};   % 转为矩阵

% ===== 输出变量初始化 =====

RFnano = nan(size(Fnano,1),size(Fnano,2));
SigMask   = nan(size(Fnano,1),size(Fnano,2));

% ===== 找所有发生过MHW的网格 =====
loc_full = unique(mhw(:,8:9),'rows');

for m = 1:size(loc_full,1)

    loc_here = loc_full(m,:);   % 当前网格点

    % ===== 当前网格所有MHW事件 =====
    mhw_here = mhw(mhw(:,8)==loc_here(1) & mhw(:,9)==loc_here(2),:);

    % ===== 每次MHW的时间段 =====
    period_mhw = [datenum(num2str(mhw_here(:,1)),'yyyymmdd'), ...
                  datenum(num2str(mhw_here(:,2)),'yyyymmdd')];

    
    rfnano = nan(size(period_mhw,1),1);  %一定要初始化矩阵，不然数据会被覆盖
    for loc = 1:size(period_mhw,1)

        % ===== 当前事件时间 =====
        mhw_time = period_mhw(loc,1):period_mhw(loc,2);

        % ===== 提取数据 =====
        mhw_nano = squeeze(Fnano(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));

        clim_nano = squeeze(mFnano(loc_here(1),loc_here(2), ...
            mhw_time - datenum(1998,1,1) + 1));

        % ===== 单次事件 anomaly =====
        
        rfnano(loc) = nanmean((mhw_nano - clim_nano).*100);
    end

    % ===== 多事件平均 =====
    
    RFnano(loc_here(1),loc_here(2))= nanmean(rfnano);


    % ===== 显著性检验 =====
    if sum(~isnan(rfnano)) > 1
        [h,p] = ttest(rfnano,0,'Alpha',0.05);

        if h==1
            SigMask(loc_here(1),loc_here(2)) = 1;
        else
            SigMask(loc_here(1),loc_here(2)) = 0;
        end
    end

end
toc

%% ===== Figure: all MHW fnano anomaly =====
tic
figure('Position',[200,200,900,800],'Color','w');

% ===== 投影 =====
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
% ===== 绘制异常 =====
m_contourf(Lon, Lat, RFnano', 80, 'linestyle', 'none'); % 注意转置
shading interp

hold on;
% ===== 显著性黑点 =====
[i, j] = find(SigMask==1);
m_scatter(Lon(i), Lat(j), 3, 'k', 'filled');  %调节点大小

% ===== 海岸线 =====
m_gshhs_i('linewidth', 1.5, 'color', 'k');
m_gshhs_i('patch', [.7 .7 .7]);  % 灰色陆地填充

% ===== 网格线 =====
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,'xtick',107:5:123,'ytick',0:5:25,'fontname','times');

% ===== 配色方案 =====
colormap(nclCM(227,18));  % 你原来的色带
caxis([-6   6]);        % 设置颜色范围

% ===== Colorbar =====
% cb = colorbar('southoutside');
% set(cb,'YTick',-6:2:6);  
% set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k');
% set(get(cb,'ylabel'),'string','Fnano anomaly (%)','FontSize',20,'FontName','Times New Roman');
% title(cb,'(%)','FontSize',16,'FontName','Times New Roman');
% title('Fnano Anomaly during MHWs','FontSize',18,'FontWeight','bold');

%wscs上升流系统区域
hold on
% ===== 红色区域框 =====
lon_box = [107.5 113.5 113.5 107.5  107.5];
lat_box = [8.5   8.5  13.5   13.5  8.5];
m_line(lon_box, lat_box, 'color','red','linewidth',2);
%luzon上升流系统区域
hold on
% ===== 红色区域框 =====
lon_box = [117   121   121  117  117];
lat_box = [16.5    16.5    21    21   16.5];
m_line(lon_box, lat_box, 'color','k','linewidth',2);


%% 提取数据

% -------- 区域1 --------
lon1_min = 109; lon1_max = 111;
lat1_min = 10;   lat1_max = 12;

ix1 = find(Lon >= lon1_min & Lon <= lon1_max);
iy1 = find(Lat >= lat1_min & Lat <= lat1_max);

RF_region1 = RFnano(ix1, iy1);
RF_region1 = RF_region1(:);          % 变列向量
RF_region1 = RF_region1(~isnan(RF_region1));  % 去除 NaN

mean1 = mean(RF_region1);
std1  = std(RF_region1,1);
min1  = min(RF_region1);

% -------- 区域2 --------
lon2_min = 118.5;   lon2_max = 121;
lat2_min = 18;  lat2_max = 20;

ix2 = find(Lon >= lon2_min & Lon <= lon2_max);
iy2 = find(Lat >= lat2_min & Lat <= lat2_max);

RF_region2 = RFnano(ix2, iy2);
RF_region2 = RF_region2(:);
RF_region2 = RF_region2(~isnan(RF_region2));

mean2 = mean(RF_region2);
std2  = std(RF_region2,1);   
min2  = min(RF_region2);


