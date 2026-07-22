%% Cellular pigmentation 细胞着色度，值越高表示细胞内chla含量越低

clear;clc;close all;

load cellpigment_025.mat   %加载细胞着色度数据
cellpig = permute(cellpig,[2 1 3]); %转置数据 交换维度，把叶绿素矩阵正常投影

load climCellpig.mat       %加载叶绿素每日气候态数据；climchla为变量名
cellpig_clim=permute(cellpig_clim,[2,1,3]);

studytime=datenum(1998,1,1):datenum(2024,12,31);  %研究时间
period_plot_v=datevec(studytime);   %对齐 daily climatology
period_unique=datevec(datenum(2016,1,1):datenum(2016,12,31));  %构建如闰年DOY序列
[~,loc_plot]=ismember(period_plot_v(:,2:3),period_unique(:,2:3),'rows');  %取研究时间的月日然后把他投射到闰年doy序列
mCellpig=squeeze(cellpig_clim(:,:,loc_plot));  %1998–2024 每一天对应的 climatology Chla,同时用于去掉长度为 1 的维度。

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

for ss=1:size(season,1)        %seasons是一个4*3的矩阵，返回行数 第一行春季（3 4 5月）  按循环四季
    loc_plot = MHWs_Season_judge(mhw,season(ss,:));  %读取mhw矩阵，调用函数 把mhw事件所属季节位置索引出来
    mhw_season = mhw(loc_plot,:);   %把热浪所属季节索引出来

    loc_full=unique(mhw_season(:,8:9),'rows');   %找到属于季节的mhw网格点，去掉重复的网格点

    RCellpig = nan(size(cellpig,1),size(cellpig,2));  %初始化异常值矩阵
    VCellpig = nan(size(cellpig,1),size(cellpig,2));
    SigMask  = nan(size(cellpig,1),size(cellpig,2));

    for m=1:size(loc_full,1)             %逐个处理每个网格点
        loc_here=loc_full(m,:);          %当前网格位置
        % 网格点对应的mhw信息
        mhw_here = mhw_season(mhw_season(:,8)==loc_here(1) & mhw_season(:,9)==loc_here(2),:);
        % 找出这个网格所有 MHW
        % 提取每次mhw时间
        period_mhw = [datenum(num2str(mhw_here(:,1)),'yyyymmdd'),...
            datenum(num2str(mhw_here(:,2)),'yyyymmdd')];

        vcellpig=nan(size(period_mhw,1),1);  %初始化空矩阵
        rcellpig=nan(size(period_mhw,1),1);  %初始化空矩阵
        for loc = 1:size(period_mhw,1)
            % 时间段筛选
            mhw_time = period_mhw(loc,1):period_mhw(loc,2);

            mhw_cellpig = squeeze(cellpig(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));  %提取mhw期间chla   提取单次mhw事件的chla

            clim_cellpig = squeeze(mCellpig(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));  %提取气候态chla     提取单次mhw事件期间的气候态chla
        
            rcellpig(loc) = nanmean(((mhw_cellpig-clim_cellpig)./clim_cellpig).*100);  %计算异常百分比 计算单次热浪chla异常

            vcellpig(loc) = nanmean(mhw_cellpig-clim_cellpig);  %计算异常值
        end

        %计算异常值的百分比
        RCellpig(loc_here(1),loc_here(2)) = nanmean(rcellpig);

        % 对异常值进行平均 
      
        VCellpig(loc_here(1),loc_here(2)) = nanmean(vcellpig);

        % t-test
        if sum(~isnan(vcellpig)) > 1
            [h,p] = ttest(vcellpig,0,'Alpha',0.05);
            if h==1
                SigMask(loc_here(1),loc_here(2)) = 1;
            else
                SigMask(loc_here(1),loc_here(2)) = 0;
            end
        end

    end

    MHWs_Cellpig_Season{ss} = VCellpig;      %Cellular pigment的异常值

    MHWs_Cellpig_per_Season{ss} = RCellpig;  %Cellular pigment异常的百分比

    SigMask_Season{ss} = SigMask;

end
toc

%% figure Cellular pigmentation spring

figure('position',[200,200,1000,900],'Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_Cellpig_Season{1}',100,'linestyle','none');
shading interp
hold on;
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Season{1}==1);
m_scatter(Lon(i),Lat(j),3,'k','filled');
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(160, 18));
caxis([-25,25]);
% cb=colorbar("southoutside");
% set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k','ytick',[-25:10:25]);
% set(get(cb,'ylabel'),'string','θ anomaly (mgC (mg Chl-a)^{-1})','FontSize',24,'FontName','Times New Roman');  %生成colorbar标题，单位

%%% Cellular pigmentation summer
figure('position',[200,200,1000,900],'Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_Cellpig_Season{2}',100,'linestyle','none');
shading interp
hold on
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Season{2}==1);
m_scatter(Lon(i),Lat(j),3,'k','filled');  
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(160, 18));
caxis([-25,25]);
% cb=colorbar("southoutside");
% set(cb,'location','eastoutside','linewidth',1.5,'fontsize',16,'edgecolor','k','ytick',[-25:10:25]);
% set(get(cb,'ylabel'),'string','θ anomaly (mgC (mg Chl-a)^{-1})','FontSize',24,'FontName','Times New Roman');  %生成colorbar标题，单位
%wscs上升流系统区域
hold on
% ===== 红色区域框 =====
lon_box = [107.5 113.5 113.5 107.5  107.5];
lat_box = [8.5   8.5    13.5   13.5  8.5];
m_line(lon_box, lat_box, 'color','red','linewidth',2);


%%% Cellular pigmentation autumn
figure('position',[200,200,1000,900],'Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_Cellpig_Season{3}',100,'linestyle','none');
shading interp
hold on
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Season{3}==1);
m_scatter(Lon(i),Lat(j),3,'k','filled');
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(160, 18));
caxis([-25,25]);
% cb=colorbar("southoutside");
% set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k','ytick',[-25:10:25]);
% set(get(cb,'ylabel'),'string','θ anomaly (mgC (mg Chl-a)^{-1})','FontSize',24,'FontName','Times New Roman');  %生成colorbar标题，单位


%%% Cellular pigmentation winter
figure('position',[200,200,1000,900],'Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_Cellpig_Season{4}',100,'linestyle','none');
shading interp
hold on
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Season{4}==1);   
m_scatter(Lon(i),Lat(j),3,'k','filled');  
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(160, 18));
caxis([-25,25]);
cb=colorbar("southoutside");
set(cb,'linewidth',1.5,'fontsize',16,'edgecolor','k');
set(get(cb,'ylabel'),'string','θ anomaly (mgC (mg Chl-a)^{-1})','FontSize',22,'FontName','Times New Roman');  %生成colorbar标题，单位
%luzon上升流系统区域
hold on
% ===== 红色区域框 =====
lon_box = [117   121   121  117  117];
lat_box = [16.5  16.5  21    21   16.5];
m_line(lon_box, lat_box, 'color','k','linewidth',2);


%% 提取数据
% -------- 区域A --------四季
lon1_min = 109; lon1_max = 111;
lat1_min = 10;   lat1_max = 12;
ix1 = find(Lon >= lon1_min & Lon <= lon1_max);
iy1 = find(Lat >= lat1_min & Lat <= lat1_max);

%夏季
A_region2 = MHWs_Cellpig_Season{2}(ix1, iy1);
A_region2 = A_region2(:);          % 变列向量
A_region2 = A_region2(~isnan(A_region2));  % 去除 NaN


% -------- 区域B --------四季
lon2_min = 118.5;   lon2_max = 121;
lat2_min = 18;  lat2_max = 20;
ix2 = find(Lon >= lon2_min & Lon <= lon2_max);
iy2 = find(Lat >= lat2_min & Lat <= lat2_max);

%冬季
B_region4 = MHWs_Cellpig_Season{4}(ix2, iy2);
B_region4 = B_region4(:);
B_region4 = B_region4(~isnan(B_region4));

%% ===== A_region2 =====
nA = numel(A_region2);    % 总有效数据个数

A_pos = sum(A_region2 > 0);   % 正值个数


A_pos_pct = A_pos / nA * 100;


fprintf('A_region2:\n');
fprintf('总样本数 = %d\n', nA);
fprintf('正值: %d (%.2f%%)\n', A_pos, A_pos_pct);

% ===== B_region4 =====
nB = numel(B_region4);
B_neg = sum(B_region4 < 0);

B_neg_pct = B_neg / nB * 100;

fprintf('B_region4:\n');
fprintf('总样本数 = %d\n', nB);
fprintf('负值: %d (%.2f%%)\n', B_neg, B_neg_pct);

