%% An example analysing seasonality of MHWs
% Here we provide an example about seasonality of MHWs
%在此模板中，我们提供了一个示例，其中包含对海洋热浪（MHWs）月度和季节性变化的简单分析。
%% 1. Loading data
clc,clear;close all;

% % Load NOAA OI SST V2 data
% sst_full=NaN(100,100,datenum(2025,12,31)-datenum(1982,1,1)+1);
% for i=1982:2025;
%     file_here=['sst_' num2str(i)];
%     load(file_here);
%     eval(['data_here=sst_' num2str(i) ';'])
%     sst_full(:,:,(datenum(i,1,1):datenum(i,12,31))-datenum(1982,1,1)+1)=data_here;
% end
% 
% sst_full(sst_full < 0 | sst_full > 45) = NaN;  %把温度FillValue去掉
% This data includes SST in [100-125E, 0-25N] in resolution of 0.25 from
% 1982 to 2025.

load sst_19822025.mat;
load sst_Lon_Lat.mat;

%% 2. Detecting MHWs and MCSs

% Here we detect marine heatwaves off eastern Tasmania based on the
% traditional definition of MHWs (Hobday et al. 2016). We detected MHWs
% during 2003 to 2025 for climatologies and thresholds in 1996 to 2025.
tic
[MHW,mclim,m90,mhw_ts]=detect(sst_full,datenum(1982,1,1):datenum(2025,12,31),datenum(1995,1,1),datenum(2024,12,31),datenum(1998,1,1),datenum(2024,12,31)); %take about 30 seconds.
toc
%检测热浪事件，sst数据源时间为1982-2025，气候基期为1995-2024，检测热浪时间范围是1998-2024.
%% 3. Generating monthly and seasonal MHW metrics
% Here we calculate monthly and seasonal MHW metrics including numbers of
% MHW days and mean MHW intensity

% Generating date matrix
date_used=datevec(datenum(1998,1,1):datenum(2024,12,31));   %热浪检测时间范围

% Determining land index
land_index=isnan(nanmean(mhw_ts,3));  %把第三维平均nanmean处理后得到二维，然后把陆地的部分筛选出来

% %% 4.Monthly
% mhwday_month=NaN(size(mhw_ts,1),size(mhw_ts,2),12);      % lon-lat-month
% mhwday_month_sum=NaN(size(mhw_ts,1),size(mhw_ts,2),12);  % lon-lat-month
% mhwint_month=NaN(size(mhw_ts,1),size(mhw_ts,2),12);      % lon-lat-month
% for i=1:12
%     index_used=date_used(:,2)==i;  %判断月份这一列是否和i相等，就是把同一个月份的全部取出来 
%     mhwday_month(:,:,i)=sum(~isnan(mhw_ts(:,:,index_used)),3,'omitnan')./(2024-1998+1); % ~isnan取出所有年份的某个月的逐日数据,是 MHW 的日子 → 1,对时间维求和，在1998–2024所有年份中，得到该月的 MHW 的总天数。符合条件的在某个月为true，不符合的为0;除以年份得到月平均mhw天数
%     mhwday_month_sum(:,:,i)=sum(~isnan(mhw_ts(:,:,index_used)),3,'omitnan'); %这里是指所有年份1-12个月的mhw天数各个月份的mhw天数之和
%     mhwint_month(:,:,i)=mean(mhw_ts(:,:,index_used),3,'omitnan'); %月平均 MHW 强度，这里直接把所有年份某个月份mhw强度求平均，所以已经得到月平均mhw强度(该月所有 MHW 日的平均强度)。 
% 
% 
% end                                 
% mhwday_month(repmat(land_index,1,1,12))=nan; %MATLAB不允许直接用二维逻辑索引去索引三维数组，用repmat在第三维复制12份，把陆地设置为nan

%% 5.Seasonal
% 北半球春季3 4 5；夏季6 7 8；秋季9 10 11；冬季12 1 2

seas=[3 4 5;...
    6 7 8;...
    9 10 11;...
    12 1 2];

mhwday_season4=NaN(size(mhw_ts,1),size(mhw_ts,2),4);
mhwint_seas=NaN(size(mhw_ts,1),size(mhw_ts,2),4); % lon-lat-seasons
nyear=2024-1998+1;
for i=1:4
    index_used=ismember(date_used(:,2),seas(i,:)); %把所有月份（的所有天）判断和归类是否在某一个季节，seas(i,:)按i行所有列索引月份的所属季节，例如：第一行属于春季（3 4 5月）

    mhwday_season4(:,:,i)=sum(~isnan(mhw_ts(:,:,index_used)),3,'omitnan')./nyear; %得到四个季节的平均mhw天数（1998-2024）;如果去掉./(2024-1998+1)则为1998-2024年分别总计四个季节的mhw天数

    mhwint_seas(:,:,i)=mean(mhw_ts(:,:,index_used),3,'omitnan');  %这里是把所有年份的按季节对mhw强度进行平均，得到季节平均mhw强度
  
end

mhwday_season4(repmat(land_index,1,1,4))=nan;  %把陆地部分赋值nan
% mhwint_seas (^{o}C) is the average intensity of MHW days in each season
% during 1998-2024


%% 6. Visualizing seasonal MHW metrics
m_proj('miller','lon',[nanmin(Lon(:))  nanmax(Lon(:))],'lat',[nanmin(Lat(:))   nanmax(Lat(:))]);
name_used={'SPR','SUM','AUT','WIN'};

for i=1:4
    subplot(2,4,i);
    m_pcolor(Lon,Lat,(mhwint_seas(:,:,i))');
    shading interp
    m_gshhs_l('patch',[0.7 0.7 0.7]);
    m_grid('linestyle','none','box','on','tickdir','in');
    colormap(m_colmap('diverging',22));
    caxis([0 2.5]);
    s=colorbar('location','southoutside');
    title(['MHW intensity (^{o}C)-' name_used{i}]);
end

for i=1:4
    subplot(2,4,i+4);
    m_pcolor(Lon,Lat,(mhwday_season4(:,:,i))');
    shading interp
    m_gshhs_l('patch',[0.7 0.7 0.7]);
    m_grid('linestyle','none','box','on','tickdir','in');
    colormap(m_colmap('diverging',22));
    caxis([0 12]);
    s=colorbar('location','southoutside');
    title(['MHW days/Season-',name_used{i}]);
end

%%%%%注意！！这个代码是画了南海四季的mhw平均强度和天数，其结果和mhw_seasonanalysis.m这份代码画出来的相同，
%以后用mhw_seasonanalysis.m这个代码来画热浪特征的图


