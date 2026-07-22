%% STR 在MHWs 的季节性异常变化

clear;clc;close all;

load Qnet_interp025_1998_2024.mat   %加载细胞着色度数据
STR = permute(STR,[2 1 3]);       %转置数据 交换维度，把叶绿素矩阵正常投影

load Clim_Qnet_SSR_STR_SLHF_SSHF.mat   %加载叶绿素每日气候态数据；climchla为变量名
STR_clim = permute(STR_clim,[2,1,3]);

studytime=datenum(1998,1,1):datenum(2024,12,31);  %研究时间
period_plot_v=datevec(studytime);   %对齐 daily climatology
period_unique=datevec(datenum(2016,1,1):datenum(2016,12,31));  %构建如闰年DOY序列
[~,loc_plot]=ismember(period_plot_v(:,2:3),period_unique(:,2:3),'rows');  %取研究时间的月日然后把他投射到闰年doy序列
mSTR=squeeze(STR_clim(:,:,loc_plot));  %1998–2024 每一天对应的 climatology Chla,同时用于去掉长度为 1 的维度。

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
    loc_plot = MHWs_Season_judge(mhw,season(seasonal,:));  %读取mhw矩阵，调用函数 把mhw事件所属季节位置索引出来
    mhw_season = mhw(loc_plot,:);   %把热浪所属季节索引出来

    loc_full=unique(mhw_season(:,8:9),'rows');   %找到属于季节的mhw网格点，去掉重复的网格点
    
    VSTR = nan(size(STR,1),size(STR,2));
    SigMask  = nan(size(STR,1),size(STR,2));
    for m=1:size(loc_full,1)             %逐个处理每个网格点
        loc_here=loc_full(m,:);          %当前网格位置
        % 网格点对应的mhw信息
        mhw_here = mhw_season(mhw_season(:,8)==loc_here(1) & mhw_season(:,9)==loc_here(2),:);
        % 找出这个网格所有 MHW
        % 提取每次mhw时间
        period_mhw = [datenum(num2str(mhw_here(:,1)),'yyyymmdd'),...
            datenum(num2str(mhw_here(:,2)),'yyyymmdd')];

        vstr =nan(size(period_mhw,1),1);    %初始化空矩阵    
        for loc = 1:size(period_mhw,1)
            % 时间段筛选
            mhw_time = period_mhw(loc,1):period_mhw(loc,2);

            mhw_str = squeeze(STR(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));  %提取mhw期间chla   提取单次mhw事件的chla

            clim_str = squeeze(mSTR(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));  %提取气候态chla     提取单次mhw事件期间的气候态chla
            
            vstr(loc) = nanmean(mhw_str-clim_str);  %计算异常值
        end

        % 对异常值进行平均 
        VSTR(loc_here(1),loc_here(2)) = nanmean(vstr);

        % t-test
        if sum(~isnan(vstr)) > 1
            [h,p] = ttest(vstr,0,'Alpha',0.05);
            if h==1
                SigMask(loc_here(1),loc_here(2)) = 1;
            else
                SigMask(loc_here(1),loc_here(2)) = 0;
            end
        end

    end

    MHWs_STR_Season{seasonal} = VSTR;      %异常值

    SigMask_Season{seasonal} = SigMask;

end
toc

%% figure  spring

figure('Position',[200,200,900,800],'Resize','off','Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_STR_Season{1}',80,'linestyle','none');
shading flat      % 建议用flat，方便保持色阶分明
hold on
% ===== 添加等值线（无数字）=====
contour_levels = -10:5:10;   % 等值线间隔可调整
m_contour(Lon,Lat,MHWs_STR_Season{1}',contour_levels,...
    'LineColor','k',...
    'LineWidth',0.8);

hold on
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Season{1}==1);
m_scatter(Lon(i),Lat(j),3,'k','filled');
%海岸线
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
%网格
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25]);
%colorbar
colormap(nclCM(130,14));
caxis([-10, 10]);
cb=colorbar('southoutside');
set(cb,'linewidth',1.5,'fontsize',18,'edgecolor','k','ytick',[-10:5:10]);
set(get(cb,'ylabel'),'string','STR (W m^{-2})','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位
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



%%%  summer
figure('Position',[200,200,900,800],'Resize','off','Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_STR_Season{2}',80,'linestyle','none');
shading flat      % 建议用flat，方便保持色阶分明
hold on
% ===== 添加等值线（无数字）=====
contour_levels = -10:5:10;   % 等值线间隔可调整
m_contour(Lon,Lat,MHWs_STR_Season{2}',contour_levels,...
    'LineColor','k',...
    'LineWidth',0.8);

hold on
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Season{2}==1);
m_scatter(Lon(i),Lat(j),3,'k','filled'); 
%地形和经纬网
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25]);
%colorbar
colormap(nclCM(130,14));
caxis([-10, 10]);
cb=colorbar('southoutside');
set(cb,'linewidth',1.5,'fontsize',18,'edgecolor','k','ytick',[-10:5:10]);
set(get(cb,'ylabel'),'string','STR (W m^{-2})','FontSize',18,'FontName','Times New Roman');  %生成colorbar标题，单位
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



%%%  autumn
figure('Position',[200,200,900,800],'Resize','off','Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_STR_Season{3}',80,'linestyle','none');
shading flat      % 建议用flat，方便保持色阶分明
hold on
% ===== 添加等值线（无数字）=====
contour_levels = -10:5:10;   % 等值线间隔可调整
m_contour(Lon,Lat,MHWs_STR_Season{3}',contour_levels,...
    'LineColor','k',...
    'LineWidth',0.8);

hold on
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Season{3}==1);
m_scatter(Lon(i),Lat(j),3,'k','filled');

m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25]);
%colorbar
colormap(nclCM(130,14));
caxis([-10, 10]);
cb=colorbar('southoutside');
set(cb,'linewidth',1.5,'fontsize',18,'edgecolor','k','ytick',[-10:5:10]);
set(get(cb,'ylabel'),'string','STR (W m^{-2})','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位
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



%%%   winter
figure('Position',[200,200,900,800],'Resize','off','Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_STR_Season{4}',80,'linestyle','none');
shading flat      % 建议用flat，方便保持色阶分明
hold on
% ===== 添加等值线（无数字）=====
contour_levels = -10:5:10;   % 等值线间隔可调整
m_contour(Lon,Lat,MHWs_STR_Season{4}',contour_levels,...
    'LineColor','k',...
    'LineWidth',0.8);

hold on
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Season{4}==1);   
m_scatter(Lon(i),Lat(j),3,'k','filled');  

m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25]);
%colorbar
colormap(nclCM(130,14));
caxis([-10, 10]);
cb=colorbar('southoutside');
set(cb,'linewidth',1.5,'fontsize',18,'edgecolor','k','ytick',[-10:5:10]);
set(get(cb,'ylabel'),'string','STR (W m^{-2})','FontSize',18,'FontName','Times New Roman');%生成colorbar标题，单位
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


%% 提取数据
% -------- 区域A --------四季
lon1_min = 109;  lon1_max = 111;
lat1_min = 10;   lat1_max = 12;
ix1 = find(Lon >= lon1_min & Lon <= lon1_max);
iy1 = find(Lat >= lat1_min & Lat <= lat1_max);
%春季
A_region1 = MHWs_STR_Season{1}(ix1, iy1);
A_region1 = A_region1(:);          % 变列向量
A_region1 = A_region1(~isnan(A_region1));  % 去除 NaN
%夏季
A_region2 = MHWs_STR_Season{2}(ix1, iy1);
A_region2 = A_region2(:);          % 变列向量
A_region2 = A_region2(~isnan(A_region2));  % 去除 NaN
%秋季
A_region3 = MHWs_STR_Season{3}(ix1, iy1);
A_region3 = A_region3(:);          % 变列向量
A_region3 = A_region3(~isnan(A_region3));  % 去除 NaN
%冬季
A_region4 = MHWs_STR_Season{4}(ix1, iy1);
A_region4 = A_region4(:);          % 变列向量
A_region4 = A_region4(~isnan(A_region4));  % 去除 NaN


% -------- 区域B --------四季
lon2_min = 118.5;  lon2_max = 121;
lat2_min = 18;    lat2_max = 20;
ix2 = find(Lon >= lon2_min & Lon <= lon2_max);
iy2 = find(Lat >= lat2_min & Lat <= lat2_max);
%春季
B_region1 = MHWs_STR_Season{1}(ix2, iy2);
B_region1 = B_region1(:);
B_region1 = B_region1(~isnan(B_region1));
%夏季
B_region2 = MHWs_STR_Season{2}(ix2, iy2);
B_region2 = B_region2(:);
B_region2 = B_region2(~isnan(B_region2));
%秋季
B_region3 = MHWs_STR_Season{3}(ix2, iy2);
B_region3 = B_region3(:);
B_region3 = B_region3(~isnan(B_region3));
%冬季
B_region4 = MHWs_STR_Season{4}(ix2, iy2);
B_region4 = B_region4(:);
B_region4 = B_region4(~isnan(B_region4));

