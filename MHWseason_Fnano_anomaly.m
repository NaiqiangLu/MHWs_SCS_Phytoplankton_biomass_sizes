%% PFTs浮游植物功能类型 分粒径浮游植物比例季节性异常
clear;clc;close all

load Fphyto1998-2024.mat   %三组分浮游植物比例数据
Fnano(Fnano<0)=nan;   %变量去异常值
Fnano = permute(Fnano,[2 1 3]); %转置数据 交换维度，把叶绿素矩阵正常投影

load ClimFphyto1998-2024.mat     %加载三组分叶绿素每日气候态数据；climchla为变量名
climFnano=permute(climFnano,[2,1,3]);

studytime=datenum(1998,1,1):datenum(2024,12,31);  %研究时间
period_plot_v=datevec(studytime);   %对齐 daily climatology
period_unique=datevec(datenum(2016,1,1):datenum(2016,12,31));  %构建如闰年DOY序列
[~,loc_plot]=ismember(period_plot_v(:,2:3),period_unique(:,2:3),'rows');  %取研究时间的月日然后把他投射到闰年doy序列
mFnano=squeeze(climFnano(:,:,loc_plot));  %1998–2024 每一天对应的 climatology Chla,同时用于去掉长度为 1 的维度。

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

    RFNano = nan(size(Fnano,1),size(Fnano,2));  %初始化异常值矩阵
    SigMask  = nan(size(Fnano,1),size(Fnano,2));
    for m=1:size(loc_full,1)             %逐个处理每个网格点
        loc_here=loc_full(m,:);          %当前网格位置
        % 网格点对应的mhw信息
        mhw_here = mhw_season(mhw_season(:,8)==loc_here(1) & mhw_season(:,9)==loc_here(2),:);
        % 找出这个网格所有 MHW
        % 提取每次mhw时间
        period_mhw = [datenum(num2str(mhw_here(:,1)),'yyyymmdd'),...
            datenum(num2str(mhw_here(:,2)),'yyyymmdd')];

        
        rfnano  = nan(size(period_mhw,1),1);
        for loc = 1:size(period_mhw,1)

            mhw_time = period_mhw(loc,1):period_mhw(loc,2);
            tidx = mhw_time - studytime(1) + 1;

            mhw_fnano  = squeeze(Fnano(loc_here(1),loc_here(2),tidx));
            clim_fnano = squeeze(mFnano(loc_here(1),loc_here(2),tidx));

            rfnano(loc)  = nanmean ((mhw_fnano-clim_fnano).*100);
            
        end

        % 至少1个mhw事件的值才能做平均
        if sum(~isnan(rfnano)) >= 1
            RFNano(loc_here(1),loc_here(2)) = nanmean(rfnano);         
        end

        % t-test
        if sum(~isnan(rfnano)) > 1    %至少有两个有效值来做t检验
            [h,p] = ttest(rfnano,0,'Alpha',0.05);
            if h==1
                SigMask(loc_here(1),loc_here(2)) = 1;
            else
                SigMask(loc_here(1),loc_here(2)) = 0;
            end
        end

    end

    SigMask_Season{seasonal} = SigMask;
    MHWs_FNano_Season{seasonal} = RFNano;  %PFTs异常的百分比
end
toc

%% figure  PFTs spring
tic
figure('Position',[200,200,900,800],'Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_FNano_Season{1}',80,'linestyle','none');
shading interp
hold on;
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Season{1}==1);
m_scatter(Lon(i),Lat(j),3,'k','filled');
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(227,18));
caxis([-7, 7]);
% cb=colorbar;
% set(cb,'location','eastoutside','linewidth',1.5,'fontsize',16,'edgecolor','k');
% set(get(cb,'ylabel'),'string','Fpico (Xi et al.2021)','FontSize',20,'FontName','Times New Roman');%生成colorbar标题，单位
% title(cb,'(%)','FontSize',18,'FontName','Times New Roman');
% title('Fnano (Xi et al.2021)','FontSize',16,'FontName','Times New Roman');

%%% PFTs summer
figure('Position',[200,200,900,800],'Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_FNano_Season{2}',80,'linestyle','none');
shading interp
hold on
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Season{2}==1);
m_scatter(Lon(i),Lat(j),3,'k','filled');
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(227,18));
caxis([-7, 7]);
% cb=colorbar;
% set(cb,'location','eastoutside','linewidth',1.5,'fontsize',16,'edgecolor','k');
% set(get(cb,'ylabel'),'string','Fpico (Xi et al.2021)','FontSize',20,'FontName','Times New Roman');%生成colorbar标题，单位
% title(cb,'(%)','FontSize',18,'FontName','Times New Roman');
% title('Fnano (Xi et al.2021)','FontSize',16,'FontName','Times New Roman');
%wscs上升流系统区域
hold on
% ===== 红色区域框 =====
lon_box = [107.5 113.5 113.5 107.5  107.5];
lat_box = [8.5   8.5  13.5   13.5  8.5];
m_line(lon_box, lat_box, 'color','red','linewidth',2);


%%% PFTs autumn
figure('Position',[200,200,900,800],'Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_FNano_Season{3}',80,'linestyle','none');
shading interp
hold on
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Season{3}==1);
m_scatter(Lon(i),Lat(j),3,'k','filled');
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(227,18));
caxis([-7, 7]);
% cb=colorbar;
% set(cb,'location','eastoutside','linewidth',1.5,'fontsize',16,'edgecolor','k');
% set(get(cb,'ylabel'),'string','Fpico (Xi et al.2021)','FontSize',20,'FontName','Times New Roman');%生成colorbar标题，单位
% title(cb,'(%)','FontSize',18,'FontName','Times New Roman');
% title('Fnano (Xi et al.2021)','FontSize',16,'FontName','Times New Roman');

%%% PFTs winter
figure('Position',[200,200,900,800],'Color','w');
m_proj('equidistant','lon',[105 123],'lat',[0 25]);
m_contourf(Lon,Lat,MHWs_FNano_Season{4}',80,'linestyle','none');
shading interp;
hold on
% ===== 显著性黑点 =====
[i,j] = find(SigMask_Season{4}==1);
m_scatter(Lon(i),Lat(j),3,'k','filled');
m_gshhs_i('linewidth',1.5,'color','k');
m_gshhs_i('patch',[.7 .7 .7]);
m_grid('linestyle','none','linewidth',1.5,'fontsize',20,...
    'xtick',[107:5:123],'ytick',[0:5:25],'fontname','times');
colormap(nclCM(227,18));
caxis([-7, 7]);
% cb=colorbar;
% set(cb,'location','eastoutside','linewidth',1.5,'fontsize',16,'edgecolor','k');
% set(get(cb,'ylabel'),'string','Fpico (Xi et al.2021)','FontSize',20,'FontName','Times New Roman');%生成colorbar标题，单位
% title(cb,'(%)','FontSize',18,'FontName','Times New Roman');
% title('Fnano (Xi et al.2021)','FontSize',16,'FontName','Times New Roman');
%luzon上升流系统区域
hold on
% ===== 红色区域框 =====
lon_box = [117   121   121  117  117];
lat_box = [16.5    16.5    21    21   16.5];
m_line(lon_box, lat_box, 'color','k','linewidth',2);


%% 提取数据  Fnano
% -------- 区域1 --------四季
lon1_min = 109; lon1_max = 111;
lat1_min = 10;   lat1_max = 12;

ix1 = find(Lon >= lon1_min & Lon <= lon1_max);
iy1 = find(Lat >= lat1_min & Lat <= lat1_max);
%春季 
A_region1 = MHWs_FNano_Season{1}(ix1, iy1);
A_region1 = A_region1(:);          % 变列向量
A_region1 = A_region1(~isnan(A_region1));  % 去除 NaN
%夏季
A_region2 = MHWs_FNano_Season{2}(ix1, iy1);
A_region2 = A_region2(:);          % 变列向量
A_region2 = A_region2(~isnan(A_region2));  % 去除 NaN
%秋季
A_region3 = MHWs_FNano_Season{3}(ix1, iy1);
A_region3 = A_region3(:);          % 变列向量
A_region3 = A_region3(~isnan(A_region3));  % 去除 NaN
%冬季
A_region4 = MHWs_FNano_Season{4}(ix1, iy1);
A_region4 = A_region4(:);          % 变列向量
A_region4 = A_region4(~isnan(A_region4));  % 去除 NaN


% -------- 区域2 --------四季
lon2_min = 118.5;   lon2_max = 121;
lat2_min = 18;  lat2_max = 20;

ix2 = find(Lon >= lon2_min & Lon <= lon2_max);
iy2 = find(Lat >= lat2_min & Lat <= lat2_max);
%春季
B_region1 = MHWs_FNano_Season{1}(ix2, iy2);
B_region1 = B_region1(:);
B_region1 = B_region1(~isnan(B_region1));
%夏季
B_region2 = MHWs_FNano_Season{2}(ix2, iy2);
B_region2 = B_region2(:);
B_region2 = B_region2(~isnan(B_region2));
%秋季
B_region3 = MHWs_FNano_Season{3}(ix2, iy2);
B_region3 = B_region3(:);
B_region3 = B_region3(~isnan(B_region3));
%冬季
B_region4 = MHWs_FNano_Season{4}(ix2, iy2);
B_region4 = B_region4(:);
B_region4 = B_region4(~isnan(B_region4));
