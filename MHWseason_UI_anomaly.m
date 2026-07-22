%% ui在MHWs 的季节性异常变化

clear;clc;close all;

load UI&curl1998-2024.mat   %加载细胞着色度数据
load UI_Curl_clim1998_2024.mat  %加载叶绿素每日气候态数据；climchla为变量名

load CCMP_wind1998_2024.mat
load climWind_1998_2024.mat  %加载叶绿素每日气候态数据；climchla为变量名

studytime=datenum(1998,1,1):datenum(2024,12,31);  %研究时间
period_plot_v=datevec(studytime);   %对齐 daily climatology
period_unique=datevec(datenum(2016,1,1):datenum(2016,12,31));  %构建如闰年DOY序列
[~,loc_plot]=ismember(period_plot_v(:,2:3),period_unique(:,2:3),'rows');  %取研究时间的月日然后把他投射到闰年doy序列

mUI_VN = squeeze(UI_VNclim(:,:,loc_plot));    %1998–2024 每一天对应的 climatology UI,同时用于去掉长度为 1 的维度。
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

%% 匹配数据计算异常值
tic
mhw = MHW{:,:};  %把mhw事件表转化为mhw矩阵
season = [3 4 5;6 7 8;9 10 11;12 1 2];  %设置好季节矩阵

for seasonal=1:size(season,1)        %seasons是一个4*3的矩阵，返回行数 第一行春季（3 4 5月）  按循环四季
    loc_plot = MHWs_Season_judge(mhw,season(seasonal,:));  %读取mhw矩阵，调用函数 把mhw事件所属季节位置索引出来
    mhw_season = mhw(loc_plot,:);   %把热浪发生时的月份所属的季节索引出来

    loc_full=unique(mhw_season(:,8:9),'rows');   %找到属于季节的mhw网格点，去掉重复的网格点

    VUI_VN = nan(size(UI_VN,1),size(UI_VN,2));  
    VUI_LZ = nan(size(UI_VN,1),size(UI_LZ,2));
    SigMaskUI_VN  = nan(size(UI_VN,1),size(UI_VN,2));
    SigMaskUI_LZ  = nan(size(UI_LZ,1),size(UI_LZ,2));
    UWinds = nan(size(ws,1),size(ws,2));
    VWinds = nan(size(ws,1),size(ws,2));
    
    MHW_UI_VN = nan(size(UI_VN,1),size(UI_VN,2));  %VN 
    Clim_UI_VN = nan(size(UI_VN,1),size(UI_VN,2));   %VN 
    MHW_UI_LZ = nan(size(UI_LZ,1),size(UI_LZ,2));  %LZ
    Clim_UI_LZ = nan(size(UI_LZ,1),size(UI_LZ,2));   %LZ

    for m=1:size(loc_full,1)             %逐个处理每个网格点
        loc_here=loc_full(m,:);          %当前网格位置
        % 网格点对应的mhw信息
        mhw_here = mhw_season(mhw_season(:,8)==loc_here(1) & mhw_season(:,9)==loc_here(2),:);
        % 找出这个网格所有 MHW
        % 提取每次mhw时间
        period_mhw = [datenum(num2str(mhw_here(:,1)),'yyyymmdd'),...
            datenum(num2str(mhw_here(:,2)),'yyyymmdd')];

        vui_vn =nan(size(period_mhw,1),1);  %初始化空矩阵     
        vui_lz =nan(size(period_mhw,1),1);  %初始化空矩阵
        uwinds=nan(size(period_mhw,1),1);  %初始化空矩阵
        vwinds=nan(size(period_mhw,1),1);  %初始化空矩阵

        mhw_ui_VN  = nan(size(period_mhw,1),1);  %提取四季热浪期间ui VN
        clim_ui_VN = nan(size(period_mhw,1),1);  %提取四季气候态ui VN
        mhw_ui_LZ = nan(size(period_mhw,1),1);  %提取四季热浪期间ui  LZ
        clim_ui_LZ = nan(size(period_mhw,1),1);  %提取四季气候态ui  LZ

        for loc = 1:size(period_mhw,1)
            % 时间段筛选
            mhw_time = period_mhw(loc,1):period_mhw(loc,2);  %网格是同一个网格，但是要遍历所有mhw period

            mhw_ui_vn = squeeze(UI_VN(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));  %提取mhw期间chla   提取单次mhw事件的chla

            clim_ui_vn = squeeze(mUI_VN(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));  %提取气候态chla     提取单次mhw事件期间的气候态chla


            mhw_ui_lz = squeeze(UI_LZ(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));  %提取mhw期间chla   提取单次mhw事件的chla

            clim_ui_lz = squeeze(mUI_LZ(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));  %提取气候态chla     提取单次mhw事件期间的气候态chla


            mhw_u = squeeze(u(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));

            clim_u = squeeze(mU(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));

            mhw_v = squeeze(v(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));

            clim_v = squeeze(mV(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));

            uwinds(loc) = nanmean(mhw_u - clim_u);   %第 loc 次 MHW 事件的“平均 xx异常值”
            vwinds(loc) = nanmean(mhw_v - clim_v);   %求当前格点每次mhws事件的异常风向，接下来对所有mhw事件下的异常风向取平均得到该格点mhw期间平均的异常风向


            vui_vn(loc) = nanmean(mhw_ui_vn-clim_ui_vn);  %计算异常值 VN

            vui_lz(loc) = nanmean(mhw_ui_lz-clim_ui_lz);  %计算异常值 LZ

            % MHW期间平均UI
            mhw_ui_VN(loc) = nanmean(mhw_ui_vn); %vn
            mhw_ui_LZ(loc) = nanmean(mhw_ui_lz); %lz
            % 气候态平均UI
            clim_ui_VN(loc) = nanmean(clim_ui_vn);  %vn
            clim_ui_LZ(loc) = nanmean(clim_ui_lz);  %lz

        end 

        % 对异常值进行平均
        VUI_VN(loc_here(1),loc_here(2)) = nanmean(vui_vn);
        VUI_LZ(loc_here(1),loc_here(2)) = nanmean(vui_lz);

        UWinds(loc_here(1),loc_here(2)) = nanmean(uwinds);
        VWinds(loc_here(1),loc_here(2)) = nanmean(vwinds);

        MHW_UI_VN(loc_here(1),loc_here(2)) = nanmean(mhw_ui_VN);   %提取四季热浪期间UI VN
        Clim_UI_VN(loc_here(1),loc_here(2)) = nanmean(clim_ui_VN);  %提取四季气候态UI VN

        MHW_UI_LZ(loc_here(1),loc_here(2)) = nanmean(mhw_ui_LZ);   %提取四季热浪期间UI LZ
        Clim_UI_LZ(loc_here(1),loc_here(2)) = nanmean(clim_ui_LZ);  %提取四季气候态UI LZ

        % t-test  VN
        if sum(~isnan(vui_vn)) > 1
            [h,p] = ttest(vui_vn,0,'Alpha',0.05);
            if h==1
                SigMaskUI_VN(loc_here(1),loc_here(2)) = 1;
            else
                SigMaskUI_VN(loc_here(1),loc_here(2)) = 0;
            end
        end

        % t-test  LZ
        if sum(~isnan(vui_lz)) > 1
            [h,p] = ttest(vui_lz,0,'Alpha',0.05);
            if h==1
                SigMaskUI_LZ(loc_here(1),loc_here(2)) = 1;
            else
                SigMaskUI_LZ(loc_here(1),loc_here(2)) = 0;
            end
        end

    end

    MHWs_UI_VN_Season{seasonal} = VUI_VN;     %ui VN的异常值   
    MHWs_UI_LZ_Season{seasonal} = VUI_LZ;     %ui LZ的异常值

    SigMaskUI_VN_Season{seasonal} = SigMaskUI_VN;  %异常显著性 VN
    SigMaskUI_LZ_Season{seasonal} = SigMaskUI_LZ;  %异常显著性 LZ

    MHWs_U_Season{seasonal} = UWinds;   %u分量异常
    MHWs_V_Season{seasonal} = VWinds;  %v分量异常

    %VN 的ui和气候态ui
    MHW_ui_VN_Season{seasonal} = MHW_UI_VN;
    Clim_ui_VN_Season{seasonal} = Clim_UI_VN;
    %LZ的ui和气候态ui值
    MHW_ui_LZ_Season{seasonal} = MHW_UI_LZ;
    Clim_ui_LZ_Season{seasonal} = Clim_UI_LZ;

end
toc

% MHWs_UI_VN_Season 是 1x4 cell，每个元素是二维矩阵
for s = 1:4
    data = MHWs_UI_VN_Season{s};   % 取出该季节数据
    % 限制异常值：绝对值大于x的设为 NaN
    data(abs(data) > 5) = NaN;
    % 放回 cell
    MHWs_UI_VN_Season{s} = data;
end

%% figure UI VN spring
tic
figure('Position',[200,200,900,800],'Resize','off','Color','w');
m_proj('equidistant','lon',[105 123],'lat',[2 25]);
m_contourf(Lon,Lat,MHWs_UI_VN_Season{1}',50,'linestyle','none');
shading interp
% ===== 显著性黑点 =====
hold on
[i,j] = find(SigMaskUI_VN_Season{1}==1);
m_scatter(Lon(i),Lat(j),2,'k','filled');
%地形
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[0.8 0.8 0.8]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(141,14));
caxis([-4, 4]);
cb=colorbar('southoutside');
set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k','fontname','times');
set(get(cb,'ylabel'),'string','UI anomaly (m^{2} s^{-1})','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位


%%%  UI VN summer
[Lon2, Lat2]=meshgrid(Lon,Lat);
figure('Position',[200,200,900,800],'Resize','off','Color','w');
m_proj('equidistant','lon',[106 116],'lat',[5 17]);
m_contourf(Lon,Lat,MHWs_UI_VN_Season{2}',50,'linestyle','none');
shading interp
%风向
hold on
skip = 3;  % 控制箭头稀疏程度
m_quiver(Lon2(1:skip:end,1:skip:end), Lat2(1:skip:end,1:skip:end), ...
    MHWs_U_Season{2}(1:skip:end,1:skip:end)', ...
    MHWs_V_Season{2}(1:skip:end,1:skip:end)', ...
    'k','linewidth',1.5);
% ===== 显著性黑点 =====
hold on
[i,j] = find(SigMaskUI_VN_Season{2}==1);
m_scatter(Lon(i),Lat(j),3,'k','filled');
%地形
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[105:3:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(141,14));
caxis([-4, 4]);
cb=colorbar('southoutside');
set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k','fontname','times');
set(get(cb,'ylabel'),'string','UI anomaly (m^{2} s^{-1})','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位
%===== wscs上升流系统区域 ====
hold on
% ===== 红色区域框 =====
lon_box = [107.5 113.5 113.5 107.5 107.5];
lat_box = [8.5   8.5   13.5   13.5   8.5];
m_line(lon_box, lat_box, 'color','red','linewidth',1.5);


%%%  UI VN autumn
figure('Position',[200,200,900,800],'Resize','off','Color','w');
m_proj('equidistant','lon',[105 123],'lat',[2 25]);
m_contourf(Lon,Lat,MHWs_UI_VN_Season{3}',50,'linestyle','none');
shading interp
% ===== 显著性黑点 =====
hold on
[i,j] = find(SigMaskUI_VN_Season{3}==1);
m_scatter(Lon(i),Lat(j),2,'k','filled');
%地形
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(141,14));
caxis([-4,4]);
cb=colorbar('southoutside');
set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k','fontname','times');
set(get(cb,'ylabel'),'string','UI anomaly (m^{2} s^{-1})','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位


%%%  UI VN winter
figure('Position',[200,200,900,800],'Resize','off','Color','w');
m_proj('equidistant','lon',[105 123],'lat',[2 25]);
m_contourf(Lon,Lat,MHWs_UI_VN_Season{4}',50,'linestyle','none');
shading interp
% ===== 显著性黑点 =====
hold on
[i,j] = find(SigMaskUI_VN_Season{4}==1);   
m_scatter(Lon(i),Lat(j),2,'k','filled');  
%地形
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(141,14));
caxis([-4, 4]);
cb=colorbar('southoutside');
set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k','fontname','times');
set(get(cb,'ylabel'),'string','UI anomaly (m^{2} s^{-1})','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位
%====== luzon上升流系统区域 ======
hold on
% ===== 红色区域框 =====
lon_box = [117   121    121 117  117];
lat_box = [16.5  16.5   21  21  16.5];
m_line(lon_box, lat_box, 'color','k','linewidth',1.5);


%%  UI LZ winter
[Lon2, Lat2]=meshgrid(Lon,Lat);

figure('Position',[200,200,900,800],'Resize','off','Color','w');
m_proj('equidistant','lon',[113 123],'lat',[13 25]);
m_contourf(Lon,Lat,MHWs_UI_LZ_Season{4}',80,'linestyle','none');
shading interp
%风向
hold on
skip = 3;  % 控制箭头稀疏程度
m_quiver(Lon2(1:skip:end,1:skip:end), Lat2(1:skip:end,1:skip:end), ...
    MHWs_U_Season{4}(1:skip:end,1:skip:end)', ...
    MHWs_V_Season{4}(1:skip:end,1:skip:end)', ...
    'k','linewidth',1.5);

hold on
% ===== 显著性黑点 =====
[i,j] = find(SigMaskUI_LZ_Season{4}==1);   
m_scatter(Lon(i),Lat(j),3,'k','filled');  
%地形
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[105:4:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(141,14));
caxis([-4, 4]);
cb=colorbar('southoutside');
set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k','fontname','times');
set(get(cb,'ylabel'),'string','UI anomaly (m^{2} s^{-1})','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位
%====== luzon上升流系统区域 ======
hold on
% ===== 红色区域框 =====
lon_box = [117   121    121 117  117];
lat_box = [16.5  16.5   21  21  16.5];
m_line(lon_box, lat_box, 'color','k','linewidth',1.5);


%% 提取ui数据
% -------- 区域A --------四季
lon1_min = 109; lon1_max = 111;
lat1_min = 10;   lat1_max = 12;
ix1 = find(Lon >= lon1_min & Lon <= lon1_max);
iy1 = find(Lat >= lat1_min & Lat <= lat1_max);

%夏季 MHW期间
MHW_A_region2 = MHW_ui_VN_Season{2}(ix1, iy1);
MHW_A_region2 = MHW_A_region2(:);          % 变列向量
MHW_A_region2 = MHW_A_region2(~isnan(MHW_A_region2));  % 去除 NaN
%夏季气候态
Clim_A_region2 = Clim_ui_VN_Season{2}(ix1, iy1);
Clim_A_region2 = Clim_A_region2(:);          % 变列向量
Clim_A_region2 = Clim_A_region2(~isnan(Clim_A_region2));  % 去除 NaN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- 区域B --------四季
lon2_min = 118.5;   lon2_max = 121;
lat2_min = 18;  lat2_max = 20;
ix2 = find(Lon >= lon2_min & Lon <= lon2_max);
iy2 = find(Lat >= lat2_min & Lat <= lat2_max);

%冬季 MHW期间
MHW_B_region4 = MHW_ui_LZ_Season{4}(ix2, iy2);
MHW_B_region4 = MHW_B_region4(:);
MHW_B_region4 = MHW_B_region4(~isnan(MHW_B_region4));
%冬季气候态
Clim_B_region4 = Clim_ui_LZ_Season{4}(ix2, iy2);
Clim_B_region4 = Clim_B_region4(:);
Clim_B_region4 = Clim_B_region4(~isnan(Clim_B_region4));

