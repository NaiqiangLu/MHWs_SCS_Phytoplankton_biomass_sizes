%海洋热浪（MHWs）进行季节性变化的简单分析。
%% 1. Loading data
clc,clear;close all;

load sst_19822025.mat;
load sst_Lon_Lat.mat;

%% 2. Detecting MHWs and MCSs

% traditional definition of MHWs (Hobday et al. 2016). We detected MHWs
% during 1998 to 2024 for climatologies and thresholds in 1995 to 2024.
tic
[MHW,mclim,m90,mhw_ts]=detect(sst_full,datenum(1982,1,1):datenum(2025,12,31),datenum(1995,1,1),datenum(2024,12,31),datenum(1998,1,1),datenum(2024,12,31)); 
toc
%检测热浪事件，sst数据源时间为1982-2025，气候基期为1995-2024，热浪检测时间范围是1998-2024.

%% 3.0 Generating  seasonal MHW metrics
% Here we calculate seasonal MHW metrics including numbers of 
% MHW days and mean MHW intensity

seas = [3 4 5;    % SPR
    6 7 8;        % SUM
    9 10 11;      % AUT
    12 1 2];      % WIN
name_used = {'Spring','Summer','Autumn','Winter'}; %分为四个季节并命名

date_used = datevec(datenum(1998,1,1):datenum(2024,12,31));
nyears = 2024 - 1998 + 1;

land_index = isnan(nanmean(mhw_ts,3));

mhwday_seas = NaN(size(mhw_ts,1),size(mhw_ts,2),4);

for i = 1:4
    index_used = ismember(date_used(:,2), seas(i,:));
    mhwday_seas(:,:,i) = ...
        sum(~isnan(mhw_ts(:,:,index_used)),3,'omitnan')./nyears;
end

mhwday_seas(repmat(land_index,1,1,4)) = NaN;

mhwint_seas = NaN(size(mhw_ts,1),size(mhw_ts,2),4);

for i = 1:4
    index_used = ismember(date_used(:,2), seas(i,:));
    mhwint_seas(:,:,i) = mean(mhw_ts(:,:,index_used),3,'omitnan');
end

mhwint_seas(repmat(land_index,1,1,4)) = NaN;

%% 3.1 Generating seasonal MHW event metrics

nx = size(mhw_ts,1);
ny = size(mhw_ts,2);

% 初始化
freq_seas   = NaN(nx,ny,4);
dur_seas    = NaN(nx,ny,4);
cumint_seas = NaN(nx,ny,4);
maxint_seas = NaN(nx,ny,4);

% 事件总数
nevents = height(MHW);

% 用 cell 暂存
freq_tmp   = zeros(nx,ny,4);
dur_tmp    = cell(nx,ny,4);
cumint_tmp = cell(nx,ny,4);
maxint_tmp = cell(nx,ny,4);

for i = 1:nevents

    % 事件开始日期  “Seasonal event metrics are assigned based on the onset month
    % of each MHW event.”  mhw平均频率 持续时间 累计强度 最大强度 用mhw onset来进行统计是科学合理 主流的做法
    t = datenum(num2str(MHW.mhw_onset(i)),'yyyymmdd');   
    m = month(datetime(t,'ConvertFrom','datenum'));

    % 判断属于哪个季节
    if ismember(m,[3 4 5])
        s = 1;
    elseif ismember(m,[6 7 8])
        s = 2;
    elseif ismember(m,[9 10 11])
        s = 3;
    else
        s = 4;
    end

    ix = MHW.xloc(i);
    iy = MHW.yloc(i);

    % frequency
    freq_tmp(ix,iy,s) = freq_tmp(ix,iy,s) + 1;

    % duration
    dur_tmp{ix,iy,s}(end+1) = MHW.mhw_dur(i);

    % cumulative intensity
    cumint_tmp{ix,iy,s}(end+1) = MHW.int_cum(i);

    % max intensity
    maxint_tmp{ix,iy,s}(end+1) = MHW.int_max(i);

end

% 计算平均值
for i=1:nx
    for j=1:ny
        for s=1:4

            freq_seas(i,j,s) = freq_tmp(i,j,s)./nyears;   % freq_seas: seasonal MHW frequency (events per season)
%“Seasonal MHW frequency is defined as the number of events per year in each season (equivalent to events per season).”
            if ~isempty(dur_tmp{i,j,s})
                dur_seas(i,j,s) = mean(dur_tmp{i,j,s});
            end

            if ~isempty(cumint_tmp{i,j,s})
                cumint_seas(i,j,s) = mean(cumint_tmp{i,j,s});
            end

            if ~isempty(maxint_tmp{i,j,s})
                maxint_seas(i,j,s) = mean(maxint_tmp{i,j,s});
            end

        end
    end
end

% 去掉陆地
freq_seas(repmat(land_index,1,1,4))   = NaN;
dur_seas(repmat(land_index,1,1,4))    = NaN;
cumint_seas(repmat(land_index,1,1,4)) = NaN;
maxint_seas(repmat(land_index,1,1,4)) = NaN;


%% ==========================================================
%  Extract every MHW event using MHWs_Season_judge
% ==========================================================

mhw = MHW{:,:};

%% -------- Region index --------

x1 = find(Lon>=109 & Lon<=111);
y1 = find(Lat>=10  & Lat<=12);

x2 = find(Lon>=118.5 & Lon<=121);
y2 = find(Lat>=18    & Lat<=20);

nx1 = length(x1);
ny1 = length(y1);

nx2 = length(x2);
ny2 = length(y2);

season = [3 4 5;
          6 7 8;
          9 10 11;
          12 1 2];

%% ==========================================================
% Region1  Summer
% ==========================================================

loc_plot = MHWs_Season_judge(mhw,season(2,:));

summer_event = [];

for i = 1:length(loc_plot)

    k = loc_plot(i);

    if ismember(MHW.xloc(k),x1) && ismember(MHW.yloc(k),y1)

        summer_event(end+1)=k;

    end

end

Ns = length(summer_event);

Region1_Frequency = NaN(nx1,ny1,Ns);
Region1_Duration  = NaN(nx1,ny1,Ns);
Region1_CumInt    = NaN(nx1,ny1,Ns);
Region1_MaxInt    = NaN(nx1,ny1,Ns);
Region1_MeanInt   = NaN(nx1,ny1,Ns);
Region1_TotalDays = NaN(nx1,ny1,Ns);

for n=1:Ns

    k=summer_event(n);

    ix=find(x1==MHW.xloc(k));
    iy=find(y1==MHW.yloc(k));

    Region1_Frequency(ix,iy,n)=1;

    Region1_Duration(ix,iy,n)=MHW.mhw_dur(k);

    Region1_CumInt(ix,iy,n)=MHW.int_cum(k);

    Region1_MaxInt(ix,iy,n)=MHW.int_max(k);

    Region1_MeanInt(ix,iy,n)=MHW.int_mean(k);

    Region1_TotalDays(ix,iy,n)=MHW.mhw_dur(k);

end

%% 展开

Region1_Frequency_col = Region1_Frequency(:);
Region1_Frequency_col(isnan(Region1_Frequency_col))=[];

Region1_TotalDays_col = Region1_TotalDays(:);
Region1_TotalDays_col(isnan(Region1_TotalDays_col))=[];

Region1_MeanInt_col = Region1_MeanInt(:);
Region1_MeanInt_col(isnan(Region1_MeanInt_col))=[];

Region1_Duration_col = Region1_Duration(:);
Region1_Duration_col(isnan(Region1_Duration_col))=[];

Region1_CumInt_col = Region1_CumInt(:);
Region1_CumInt_col(isnan(Region1_CumInt_col))=[];

Region1_MaxInt_col = Region1_MaxInt(:);
Region1_MaxInt_col(isnan(Region1_MaxInt_col))=[];


%% ==========================================================
% Region2  Winter
% ==========================================================

loc_plot = MHWs_Season_judge(mhw,season(4,:));

winter_event = [];

for i=1:length(loc_plot)

    k = loc_plot(i);

    if ismember(MHW.xloc(k),x2) && ismember(MHW.yloc(k),y2)

        winter_event(end+1)=k;

    end

end

Nw = length(winter_event);

Region2_Frequency = NaN(nx2,ny2,Nw);
Region2_Duration  = NaN(nx2,ny2,Nw);
Region2_CumInt    = NaN(nx2,ny2,Nw);
Region2_MaxInt    = NaN(nx2,ny2,Nw);
Region2_MeanInt   = NaN(nx2,ny2,Nw);
Region2_TotalDays = NaN(nx2,ny2,Nw);

for n=1:Nw

    k=winter_event(n);

    ix=find(x2==MHW.xloc(k));
    iy=find(y2==MHW.yloc(k));

    Region2_Frequency(ix,iy,n)=1;

    Region2_Duration(ix,iy,n)=MHW.mhw_dur(k);

    Region2_CumInt(ix,iy,n)=MHW.int_cum(k);

    Region2_MaxInt(ix,iy,n)=MHW.int_max(k);

    Region2_MeanInt(ix,iy,n)=MHW.int_mean(k);

    Region2_TotalDays(ix,iy,n)=MHW.mhw_dur(k);

end

%% 展开

Region2_Frequency_col = Region2_Frequency(:);
Region2_Frequency_col(isnan(Region2_Frequency_col))=[];

Region2_TotalDays_col = Region2_TotalDays(:);
Region2_TotalDays_col(isnan(Region2_TotalDays_col))=[];

Region2_MeanInt_col = Region2_MeanInt(:);
Region2_MeanInt_col(isnan(Region2_MeanInt_col))=[];

Region2_Duration_col = Region2_Duration(:);
Region2_Duration_col(isnan(Region2_Duration_col))=[];

Region2_CumInt_col = Region2_CumInt(:);
Region2_CumInt_col(isnan(Region2_CumInt_col))=[];

Region2_MaxInt_col = Region2_MaxInt(:);
Region2_MaxInt_col(isnan(Region2_MaxInt_col))=[];


