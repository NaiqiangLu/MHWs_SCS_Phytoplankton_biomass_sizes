%%  Curl 在 MHWs 的季节性异常变化

clear;clc;close all;

load UI&curl1998-2024.mat       %加载curl数据
load UI_Curl_clim1998_2024.mat  %加载curl每日气候态数据；Curl_clim为变量名

studytime=datenum(1998,1,1):datenum(2024,12,31);  %研究时间
period_plot_v=datevec(studytime);   %对齐 daily climatology
period_unique=datevec(datenum(2016,1,1):datenum(2016,12,31));  %构建如闰年DOY序列
[~,loc_plot]=ismember(period_plot_v(:,2:3),period_unique(:,2:3),'rows');  %取研究时间的月日然后把他投射到闰年doy序列

mCurl = squeeze(Curl_clim(:,:,loc_plot));    %1998–2024 每一天对应的 climatology curl,同时用于去掉长度为 1 的维度。

%% SST 读取海温数据和热浪
tic
load sst_19822025.mat

[MHW,smclim,m90,mhw_ts]=detect(sst_full,datenum(1982,1,1):datenum(2025,12,31),...
    datenum(1995,1,1),datenum(2024,12,31),...
    datenum(1998,1,1),datenum(2024,12,31),'Threshold',0.9);
toc

%% 匹配数据计算异常值
tic
mhw = MHW{:,:};  %把mhw事件表转化为mhw矩阵
season = [3 4 5;6 7 8;9 10 11;12 1 2];  %设置好季节矩阵

for seasonal=1:size(season,1)        %seasons是一个4*3的矩阵，返回行数 第一行春季（3 4 5月）  按循环四季
    loc_plot = MHWs_Season_judge(mhw,season(seasonal,:));  %读取mhw矩阵，调用函数 把mhw事件所属季节位置索引出来
    mhw_season = mhw(loc_plot,:);   %把热浪发生时的月份所属的季节索引出来

    loc_full=unique(mhw_season(:,8:9),'rows');   %找到属于季节的mhw网格点，去掉重复的网格点

    VCurl = nan(size(curl,1),size(curl,2));   %初始化矩阵
    SigMask_Curl  = nan(size(curl,1),size(curl,2));

    MHW_Curl = nan(size(curl,1),size(curl,2));
    Clim_Curl = nan(size(curl,1),size(curl,2));

    for m=1:size(loc_full,1)             %逐个处理每个网格点
        loc_here=loc_full(m,:);          %当前网格位置
        % 网格点对应的mhw信息
        mhw_here = mhw_season(mhw_season(:,8)==loc_here(1) & mhw_season(:,9)==loc_here(2),:);
        % 找出这个网格所有 MHW
        % 提取每次mhw时间
        period_mhw = [datenum(num2str(mhw_here(:,1)),'yyyymmdd'),...
            datenum(num2str(mhw_here(:,2)),'yyyymmdd')];

        vcurl =nan(size(period_mhw,1),1);  %初始化空矩阵

        mhw_Curlz  = nan(size(period_mhw,1),1);  %提取四季热浪期间curl
        clim_Curlz = nan(size(period_mhw,1),1);  %提取四季气候态curl

        for loc = 1:size(period_mhw,1)
            % 时间段筛选
            mhw_time = period_mhw(loc,1):period_mhw(loc,2);  %网格是同一个网格，但是要遍历所有mhw period

            mhw_curl = squeeze(curl(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));  %提取mhw期间curl   提取单次mhw事件的curl

            clim_curl = squeeze(mCurl(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));  %提取气候态curl     提取单次mhw事件期间的气候态curl

            vcurl(loc) = nanmean(mhw_curl-clim_curl);  %计算异常值

            % MHW期间平均Curl
            mhw_Curlz(loc) = nanmean(mhw_curl);
            % 气候态平均curl
            clim_Curlz(loc) = nanmean(clim_curl);

        end

        % 对异常值进行平均
        VCurl(loc_here(1),loc_here(2)) = nanmean(vcurl);

        MHW_Curl(loc_here(1),loc_here(2)) = nanmean(mhw_Curlz);   %提取四季热浪期间curl
        Clim_Curl(loc_here(1),loc_here(2)) = nanmean(clim_Curlz);  %提取四季气候态curl

        % t-test
        if sum(~isnan(vcurl)) > 1
            [h,p] = ttest(vcurl,0,'Alpha',0.05);
            if h==1
                SigMask_Curl(loc_here(1),loc_here(2)) = 1;
            else
                SigMask_Curl(loc_here(1),loc_here(2)) = 0;
            end
        end

    end

    VCurl=VCurl./1e-6;  %风应力旋度 N/m^3 （10^-6数量级）正值：气旋式风应力旋度（海洋中引发上升流）负值：反气旋式风应力旋度（海洋中引发下沉流）
    % 除以数量级方便画图的colorbar控制

    MHWs_Curl_Season{seasonal} = VCurl;     %curl 的异常值
    SigMask_Curl_Season{seasonal} = SigMask_Curl;  

    MHW_curl_Season{seasonal} = MHW_Curl;  %热浪期间
    Clim_curl_Season{seasonal} = Clim_Curl;  %气候态
end
toc

%% figure Curl spring
tic
figure('Position',[200,200,900,800],'Resize','off','Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_Curl_Season{1}',50,'linestyle','none');
shading flat      % 建议用flat，方便保持色阶分明
hold on
% ===== 添加等值线（无数字）=====
contour_levels = -0.5:0.2:0.5;   % 等值线间隔可调整
m_contour(Lon,Lat,MHWs_Curl_Season{1}',contour_levels,'LineColor','k','LineWidth',0.8);
hold on
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Curl_Season{1}==1);
m_scatter(Lon(i),Lat(j),3,'k','filled');
% 地形
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(141,14));
caxis([-0.5,  0.5]);
cb=colorbar('southoutside');
set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k','ytick',[-0.5:0.25:0.5],'fontname','times');
set(get(cb,'ylabel'),'string','Wind stress curl anomaly (10^{-6} N m^{-3})','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位
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


%%%   Curl summer
figure('Position',[200,200,900,800],'Resize','off','Color','w');
m_proj('equidistant','lon',[106 116],'lat',[5 17]);  %上升流系统局部放大图
m_contourf(Lon,Lat,MHWs_Curl_Season{2}',50,'linestyle','none');
shading flat      % 建议用flat，方便保持色阶分明
hold on
% ===== 添加等值线（无数字）=====
contour_levels = -0.5:0.2:0.5;   % 等值线间隔可调整
m_contour(Lon,Lat,MHWs_Curl_Season{2}',contour_levels,'LineColor','k','LineWidth',0.8);
hold on
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Curl_Season{2}==1);
m_scatter(Lon(i),Lat(j),3,'k','filled');
% 地形
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:4:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(141,14));
caxis([-0.6,  0.6]);
cb=colorbar('southoutside');
set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k','ytick',[-0.5:0.25:0.5],'fontname','times');
set(get(cb,'ylabel'),'string','Wind stress curl anomaly (10^{-6} N m^{-3})','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位
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
%===== wscs上升流系统区域 ====
hold on
% ===== 红色区域框 =====
lon_box = [107.5 113.5 113.5 107.5 107.5];
lat_box = [8.5   8.5   13.5   13.5   8.5];
m_line(lon_box, lat_box, 'color','red','linewidth',1.5);


%%%  Curl autumn
figure('Position',[200,200,900,800],'Resize','off','Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_Curl_Season{3}',50,'linestyle','none');
shading flat      % 建议用flat，方便保持色阶分明
hold on
% ===== 添加等值线（无数字）=====
contour_levels = -0.5:0.2:0.5;   % 等值线间隔可调整
m_contour(Lon,Lat,MHWs_Curl_Season{3}',contour_levels,'LineColor','k','LineWidth',0.8);
hold on
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Curl_Season{3}==1);
m_scatter(Lon(i),Lat(j),3,'k','filled');
% 地形
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(141,14));
caxis([-0.5, 0.5]);
cb=colorbar('southoutside');
set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k','ytick',[-0.5:0.25:0.5],'fontname','times');
set(get(cb,'ylabel'),'string','Wind stress curl anomaly (10^{-6} N m^{-3})','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位
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


%%%  Curl winter
figure('Position',[200,200,900,800],'Resize','off','Color','w');
m_proj('equidistant','lon',[113 123],'lat',[13 25]);
m_contourf(Lon,Lat,MHWs_Curl_Season{4}',50,'linestyle','none');
shading flat      % 建议用flat，方便保持色阶分明
hold on
% ===== 添加等值线（无数字）=====
contour_levels = -0.5:0.2:0.5;   % 等值线间隔可调整
m_contour(Lon,Lat,MHWs_Curl_Season{4}',contour_levels,'LineColor','k','LineWidth',0.8);
hold on
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Curl_Season{4}==1);   
m_scatter(Lon(i),Lat(j),3,'k','filled');  
% 地形
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[106:4:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(141,14));
caxis([-0.5, 0.5]);
cb=colorbar('southoutside');
set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k','ytick',[-0.5:0.25:0.5],'fontname','times');
set(get(cb,'ylabel'),'string','Wind stress curl anomaly (10^{-6} N m^{-3})','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位
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
%====== luzon上升流系统区域 ======
hold on
% ===== 红色区域框 =====
lon_box = [117   121    121 117  117];
lat_box = [16.5  16.5   21  21  16.5];
m_line(lon_box, lat_box, 'color','k','linewidth',1.5);


%%  提取 curlz  数据

% -------- 区域A --------四季
lon1_min = 109; lon1_max = 111;
lat1_min = 10;   lat1_max = 12;
ix1 = find(Lon >= lon1_min & Lon <= lon1_max);
iy1 = find(Lat >= lat1_min & Lat <= lat1_max);

%夏季 MHW期间
MHW_A_region2 = MHW_curl_Season{2}(ix1, iy1);
MHW_A_region2 = MHW_A_region2(:)./1e-6;          % 变列向量
MHW_A_region2 = MHW_A_region2(~isnan(MHW_A_region2));  % 去除 NaN
%夏季气候态
Clim_A_region2 = Clim_curl_Season{2}(ix1, iy1);
Clim_A_region2 = Clim_A_region2(:)./1e-6;          % 变列向量
Clim_A_region2 = Clim_A_region2(~isnan(Clim_A_region2));  % 去除 NaN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- 区域B --------四季
lon2_min = 118.5;   lon2_max = 121;
lat2_min = 18;  lat2_max = 20;
ix2 = find(Lon >= lon2_min & Lon <= lon2_max);
iy2 = find(Lat >= lat2_min & Lat <= lat2_max);

%冬季 MHW期间
MHW_B_region4 = MHW_curl_Season{4}(ix2, iy2);
MHW_B_region4 = MHW_B_region4(:)./1e-6;
MHW_B_region4 = MHW_B_region4(~isnan(MHW_B_region4));
%冬季气候态
Clim_B_region4 = Clim_curl_Season{4}(ix2, iy2);
Clim_B_region4 = Clim_B_region4(:)./1e-6;
Clim_B_region4 = Clim_B_region4(~isnan(Clim_B_region4));

