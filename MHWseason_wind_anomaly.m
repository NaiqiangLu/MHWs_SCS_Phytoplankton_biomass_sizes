%% WIND 季节异常
clear;clc;close all

load CCMP_wind1998_2024.mat
load climWind_1998_2024.mat      %加载叶绿素每日气候态数据；climchla为变量名

studytime=datenum(1998,1,1):datenum(2024,12,31);  %研究时间
period_plot_v = datevec(studytime);    %对齐 daily climatology
period_unique = datevec(datenum(2016,1,1):datenum(2016,12,31));  %构建如闰年DOY序列
[~,loc_plot]= ismember(period_plot_v(:,2:3),period_unique(:,2:3),'rows');  %取研究时间的月日然后把他投射到闰年doy序列

mWindS=squeeze(ws_clim(:,:,loc_plot));  %1998–2024 每一天对应的 climatology Chla,同时用于去掉长度为 1 的维度。
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
mhw = MHW{:,:};  %把mhw事件表转化为mhw矩阵
season = [3 4 5;6 7 8;9 10 11;12 1 2];  %设置好季节矩阵

tic
for seasonal=1:size(season,1)        %seasons是一个4*3的矩阵，返回行数 第一行春季（3 4 5月）  按循环四季
    loc_plot =MHWs_Season_judge(mhw,season(seasonal,:));  %读取mhw矩阵，调用函数 把mhw事件所属季节位置索引出来
    mhw_season = mhw(loc_plot,:);   %把热浪所属季节索引出来

    loc_full=unique(mhw_season(:,8:9),'rows');   %找到属于季节的mhw网格点，去掉重复的网格点

    RWinds = nan(size(ws,1),size(ws,2));  %初始化异常值矩阵
    ValWinds = nan(size(ws,1),size(ws,2));
    UWinds = nan(size(ws,1),size(ws,2));
    VWinds = nan(size(ws,1),size(ws,2));
    SigMask  = nan(size(ws,1),size(ws,2));

    MHWWindSpeed = nan(size(ws,1),size(ws,2));
    ClimWindSpeed = nan(size(ws,1),size(ws,2));

    for m=1:size(loc_full,1)             %逐个处理每个网格点
        loc_here=loc_full(m,:);          %当前网格位置
        % 网格点对应的mhw信息
        mhw_here = mhw_season(mhw_season(:,8)==loc_here(1) & mhw_season(:,9)==loc_here(2),:);
        % 找出这个网格所有 MHW
        % 提取每次mhw时间
        period_mhw = [datenum(num2str(mhw_here(:,1)),'yyyymmdd'),...
            datenum(num2str(mhw_here(:,2)),'yyyymmdd')];

        valwinds=nan(size(period_mhw,1),1);
        rwinds=nan(size(period_mhw,1),1);  %初始化空矩阵
        uwinds=nan(size(period_mhw,1),1);  %初始化空矩阵
        vwinds=nan(size(period_mhw,1),1);  %初始化空矩阵

        mhw_windspeed  = nan(size(period_mhw,1),1);  %提取四季热浪期间风速
        clim_windspeed = nan(size(period_mhw,1),1);  %提取四季气候态风速

        for loc = 1:size(period_mhw,1)
            % 时间段筛选
            mhw_time = period_mhw(loc,1):period_mhw(loc,2);

            mhw_ws = squeeze(ws(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));  %提取mhw期间chla   提取单次mhw事件的chla

            clim_ws = squeeze(mWindS(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));  %提取气候态chla     提取单次mhw事件期间的气候态chla
            mhw_ws(mhw_ws< 0) = NaN;
            clim_ws(clim_ws< 0) = NaN;

            % MHW期间平均风速
            mhw_windspeed(loc) = nanmean(mhw_ws);
            % 气候态平均风速
            clim_windspeed(loc) = nanmean(clim_ws);


            mhw_u = squeeze(u(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));

            clim_u = squeeze(mU(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));

            mhw_v = squeeze(v(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));

            clim_v = squeeze(mV(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));

            uwinds(loc) = nanmean(mhw_u - clim_u);
            vwinds(loc) = nanmean(mhw_v - clim_v);  %求mhw期间异常风向

            % uwinds(loc) = nanmean(mhw_u);
            % vwinds(loc) = nanmean(mhw_v);  %求mhw期间的平均风向

            rwinds(loc) = nanmean(((mhw_ws-clim_ws)./clim_ws)*100);  %计算异常百分比 计算单次热浪chla异常

            valwinds(loc) = nanmean(mhw_ws-clim_ws);  %计算异常值
        end

        %计算异常值的百分比
        RWinds(loc_here(1),loc_here(2)) = nanmean(rwinds);
        % 对异常值进行平均
        ValWinds(loc_here(1),loc_here(2)) = nanmean(valwinds);

        UWinds(loc_here(1),loc_here(2)) = nanmean(uwinds);
        VWinds(loc_here(1),loc_here(2)) = nanmean(vwinds);

        MHWWindSpeed(loc_here(1),loc_here(2)) = nanmean(mhw_windspeed);   %提取四季热浪期间风速
        ClimWindSpeed(loc_here(1),loc_here(2)) = nanmean(clim_windspeed);  %提取四季气候态风速

        % t-test
        if sum(~isnan(vwinds)) > 1
            [h,p] = ttest(vwinds,0,'Alpha',0.05);
            if h==1
                SigMask(loc_here(1),loc_here(2)) = 1;
            else
                SigMask(loc_here(1),loc_here(2)) = 0;
            end
        end

    end

    MHWs_WindS_val_Season{seasonal} = ValWinds;      %Cphyto的异常值
    MHWs_WindS_per_Season{seasonal} = RWinds;      %Cphyto异常的百分比

    MHWs_U_Season{seasonal} = UWinds;
    MHWs_V_Season{seasonal} = VWinds;

    SigMask_Season{seasonal} = SigMask;

    MHWs_WindSpeed_Season{seasonal} = MHWWindSpeed;
    Clim_WindSpeed_Season{seasonal} = ClimWindSpeed;

end
toc

%% figure windanomaly spring

[Lon2, Lat2] = meshgrid(Lon, Lat);

figure('Position',[200,200,900,800],'Resize','on','Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_WindS_val_Season{1}',50,'linestyle','none');

hold on
skip = 4;  % 控制箭头稀疏程度
m_quiver(Lon2(1:skip:end,1:skip:end), Lat2(1:skip:end,1:skip:end), ...
    MHWs_U_Season{1}(1:skip:end,1:skip:end)', ...
    MHWs_V_Season{1}(1:skip:end,1:skip:end)', ...
    'k','linewidth',1.5);

hold on;
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(325,14));
caxis([-3.5,3.5]);
cb=colorbar;
set(cb,'location','southoutside','linewidth',1.5,'fontsize',16,'edgecolor','k','fontname','times');
set(get(cb,'ylabel'),'string','Windspeed anomaly (m/s)','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位


%%% windanomaly summer
figure('Position',[200,200,900,800],'Resize','on','Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_WindS_val_Season{2}',50,'linestyle','none');

hold on
skip = 4;  % 控制箭头稀疏程度
m_quiver(Lon2(1:skip:end,1:skip:end), Lat2(1:skip:end,1:skip:end), ...
    MHWs_U_Season{2}(1:skip:end,1:skip:end)', ...
    MHWs_V_Season{2}(1:skip:end,1:skip:end)', ...
    'k','linewidth',1.5);

hold on
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(325,14));
caxis([-3.5,3.5]);
cb=colorbar;
set(cb,'location','southoutside','linewidth',1.5,'fontsize',16,'edgecolor','k','fontname','times');
set(get(cb,'ylabel'),'string','Windspeed anomaly (m/s)','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位
%wscs上升流系统区域
hold on
% ===== 红色区域框 =====
lon_box = [107.5 113.5 113.5 107.5 107.5];
lat_box = [8.5   8.5   13.5   13.5   8.5];
m_line(lon_box, lat_box, 'color','red','linewidth',1.5);


%%% windanomaly autumn
figure('Position',[200,200,900,800],'Resize','on','Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_WindS_val_Season{3}',50,'linestyle','none');

hold on
skip = 4;  % 控制箭头稀疏程度
m_quiver(Lon2(1:skip:end,1:skip:end), Lat2(1:skip:end,1:skip:end), ...
    MHWs_U_Season{3}(1:skip:end,1:skip:end)', ...
    MHWs_V_Season{3}(1:skip:end,1:skip:end)', ...
    'k','linewidth',1.5);

hold on
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(325,14));
caxis([-3.5,3.5]);
cb=colorbar;
set(cb,'location','southoutside','linewidth',1.5,'fontsize',16,'edgecolor','k','fontname','times');
set(get(cb,'ylabel'),'string','Windspeed anomaly (m/s)','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位


%%%  windanomaly winter
figure('Position',[200,200,900,800],'Resize','on','Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_WindS_val_Season{4}',50,'linestyle','none');

hold on
skip = 4;  % 控制箭头稀疏程度
m_quiver(Lon2(1:skip:end,1:skip:end), Lat2(1:skip:end,1:skip:end), ...
    MHWs_U_Season{4}(1:skip:end,1:skip:end)', ...
    MHWs_V_Season{4}(1:skip:end,1:skip:end)', ...
    'k','linewidth',1.5);

hold on
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(325,14));
caxis([-3.5,3.5]);
cb=colorbar;
set(cb,'location','southoutside','linewidth',1.5,'fontsize',16,'edgecolor','k','fontname','times');
set(get(cb,'ylabel'),'string','Windspeed anomaly (m/s)','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位
%luzon上升流系统区域
hold on
% ===== 红色区域框 =====
lon_box = [117   121    121 117  117];
lat_box = [16.5  16.5   21  21  16.5];
m_line(lon_box, lat_box, 'color','k','linewidth',1.5);


%% 提取风速数据 

% -------- 区域A --------四季
lon1_min = 109; lon1_max = 111;
lat1_min = 10;   lat1_max = 12;
ix1 = find(Lon >= lon1_min & Lon <= lon1_max);
iy1 = find(Lat >= lat1_min & Lat <= lat1_max);

%夏季 MHW期间
MHW_A_region2 = MHWs_WindSpeed_Season{2}(ix1, iy1);
MHW_A_region2 = MHW_A_region2(:);          % 变列向量
MHW_A_region2 = MHW_A_region2(~isnan(MHW_A_region2));  % 去除 NaN
%夏季气候态
Clim_A_region2 =Clim_WindSpeed_Season{2}(ix1, iy1);
Clim_A_region2 = Clim_A_region2(:);          % 变列向量
Clim_A_region2 = Clim_A_region2(~isnan(Clim_A_region2));  % 去除 NaN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- 区域B --------四季
lon2_min = 118.5;   lon2_max = 121;
lat2_min = 18;  lat2_max = 20;
ix2 = find(Lon >= lon2_min & Lon <= lon2_max);
iy2 = find(Lat >= lat2_min & Lat <= lat2_max);

%冬季 MHW期间
MHW_B_region4 = MHWs_WindSpeed_Season{4}(ix2, iy2);
MHW_B_region4 = MHW_B_region4(:);
MHW_B_region4 = MHW_B_region4(~isnan(MHW_B_region4));
%冬季气候态
Clim_B_region4 = Clim_WindSpeed_Season{4}(ix2, iy2);
Clim_B_region4 = Clim_B_region4(:);
Clim_B_region4 = Clim_B_region4(~isnan(Clim_B_region4));

