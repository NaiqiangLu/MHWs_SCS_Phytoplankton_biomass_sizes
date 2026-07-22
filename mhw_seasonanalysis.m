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

%% 3. Generating  seasonal MHW metrics
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


%% ================= FIGURE 1 =================
figure('Position',[50 50 1500 1000],'Color','w');
tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

panel_label = 'a';

for row = 1:3
    for col = 1:4

        nexttile

        m_proj('miller','lon',[100 123],'lat',[0 25]);

        % ===== 变量选择 =====
        if row==1
            data = freq_seas(:,:,col);
            clim = [0 1];
            varname = 'Frequency (count)';
        elseif row==2
            data = mhwint_seas(:,:,col);
            clim = [0 2.2];
            varname = 'Mean Intensity (°C)';
        else
            data = mhwday_seas(:,:,col);
            clim = [2 10];
            varname = 'Total Days (Days)';
        end

        m_pcolor(Lon,Lat,data');
        shading interp
        hold on

        m_gshhs_i('patch',[0.7 0.7 0.7],'LineWidth',1);
        m_grid('linestyle','none','box','on',...
               'tickdir','in','fontsize',13,...
               'fontname','times','linewidth',1);

        colormap(nclCM(96,22));
        caxis(clim);

        % ===== 字母标注 =====
        text(0.02,0.95,['(' panel_label ')'],...
            'Units','normalized',...
            'FontSize',19,'FontWeight','bold',...
            'FontName','times');

        panel_label = char(panel_label + 1);

        if row==1
            title(name_used{col},'fontsize',16,'fontname','times');
        end

        hold off

        % ===== 每行最后一个图加 colorbar =====
        if col==4
            cb = colorbar;
            cb.Layout.Tile = 'east';   % 放在右侧
            cb.FontSize = 14;
            cb.FontName ='times';
            cb.LineWidth=1;
            ylabel(cb,varname,'FontSize',16,...
                   'FontName','times');
        end

    end
end

%% ================= FIGURE 2 =================
figure('Position',[50 50 1500 1000],'Color','w');
tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

panel_label = 'a';

for row = 1:3
    for col = 1:4

        nexttile

        m_proj('miller','lon',[100 123],'lat',[0 25]);

        if row==1
            data = maxint_seas(:,:,col);
            clim = [0.5 2.5];
            varname = 'Maximum Intensity (°C)';
        elseif row==2
            data = cumint_seas(:,:,col);
            clim = [5 35];
            varname = 'Cumulative Intensity (°C·days)';
        else
            data = dur_seas(:,:,col);
            clim = [5 22];
            varname = 'Duration (Days)';
        end

        m_pcolor(Lon,Lat,data');
        shading interp
        hold on

        m_gshhs_i('patch',[0.7 0.7 0.7],'LineWidth',1);
        m_grid('linestyle','none','box','on',...
               'tickdir','in','fontsize',13,...
               'fontname','times','linewidth',1);

        colormap(nclCM(96,22));
        caxis(clim);

        text(0.02,0.95,['(' panel_label ')'],...
            'Units','normalized',...
            'FontSize',19,'FontWeight','bold',...
            'FontName','times');

        panel_label = char(panel_label + 1);

        if row==1
            title(name_used{col},'fontsize',16,'fontname','times');
        end

        hold off

        % ===== 每行右侧 colorbar =====
        if col==4
            cb = colorbar;
            cb.Layout.Tile = 'east';
            cb.FontSize = 14;
            cb.FontName ='times';
            cb.LineWidth=1;
            ylabel(cb,varname,'FontSize',16,...
                'FontName','times');
        end

    end
end


%% ===== 提取两个小区域四季MHW指标，因为要跟浮游植物采样的小区域相同才有可比性（不做区域平均）=====

% -------- 区域1 --------
lon1_min = 109; lon1_max = 111;
lat1_min = 10;  lat1_max = 12;

% -------- 区域2 --------
lon2_min = 118.5; lon2_max = 121;
lat2_min = 18;    lat2_max = 20;

% 经纬度索引

x1 = find(Lon>=lon1_min & Lon<=lon1_max);
y1 = find(Lat>=lat1_min & Lat<=lat1_max);

x2 = find(Lon>=lon2_min & Lon<=lon2_max);
y2 = find(Lat>=lat2_min & Lat<=lat2_max);

% ===== 区域1提取 =====

Region1_Frequency = freq_seas(x1,y1,:);
Region1_TotalDays = mhwday_seas(x1,y1,:);
Region1_MeanInt   = mhwint_seas(x1,y1,:);

Region1_Duration  = dur_seas(x1,y1,:);
Region1_CumInt    = cumint_seas(x1,y1,:);
Region1_MaxInt    = maxint_seas(x1,y1,:);

% ===== 区域2提取 =====

Region2_Frequency = freq_seas(x2,y2,:);
Region2_TotalDays = mhwday_seas(x2,y2,:);
Region2_MeanInt   = mhwint_seas(x2,y2,:);

Region2_Duration  = dur_seas(x2,y2,:);
Region2_CumInt    = cumint_seas(x2,y2,:);
Region2_MaxInt    = maxint_seas(x2,y2,:);

% ===== 展平为Excel易复制的一列数据 =====

season_names = {'Spring','Summer','Autumn','Winter'};

% ===== 区域1 =====

for s = 1:4

    % Frequency
    tmp = Region1_Frequency(:,:,s);
    Region1_Frequency_col{s} = tmp(:);
    Region1_Frequency_col{s}(isnan(Region1_Frequency_col{s}))=[];

    % Total Days
    tmp = Region1_TotalDays(:,:,s);
    Region1_TotalDays_col{s} = tmp(:);
    Region1_TotalDays_col{s}(isnan(Region1_TotalDays_col{s}))=[];

    % Mean Intensity
    tmp = Region1_MeanInt(:,:,s);
    Region1_MeanInt_col{s} = tmp(:);
    Region1_MeanInt_col{s}(isnan(Region1_MeanInt_col{s}))=[];

    % Duration
    tmp = Region1_Duration(:,:,s);
    Region1_Duration_col{s} = tmp(:);
    Region1_Duration_col{s}(isnan(Region1_Duration_col{s}))=[];

    % CumInt
    tmp = Region1_CumInt(:,:,s);
    Region1_CumInt_col{s} = tmp(:);
    Region1_CumInt_col{s}(isnan(Region1_CumInt_col{s}))=[];

    % MaxInt
    tmp = Region1_MaxInt(:,:,s);
    Region1_MaxInt_col{s} = tmp(:);
    Region1_MaxInt_col{s}(isnan(Region1_MaxInt_col{s}))=[];

end

% ===== 区域2 =====

for s = 1:4

    % Frequency
    tmp = Region2_Frequency(:,:,s);
    Region2_Frequency_col{s} = tmp(:);
    Region2_Frequency_col{s}(isnan(Region2_Frequency_col{s}))=[];

    % Total Days
    tmp = Region2_TotalDays(:,:,s);
    Region2_TotalDays_col{s} = tmp(:);
    Region2_TotalDays_col{s}(isnan(Region2_TotalDays_col{s}))=[];

    % Mean Intensity
    tmp = Region2_MeanInt(:,:,s);
    Region2_MeanInt_col{s} = tmp(:);
    Region2_MeanInt_col{s}(isnan(Region2_MeanInt_col{s}))=[];

    % Duration
    tmp = Region2_Duration(:,:,s);
    Region2_Duration_col{s} = tmp(:);
    Region2_Duration_col{s}(isnan(Region2_Duration_col{s}))=[];

    % CumInt
    tmp = Region2_CumInt(:,:,s);
    Region2_CumInt_col{s} = tmp(:);
    Region2_CumInt_col{s}(isnan(Region2_CumInt_col{s}))=[];

    % MaxInt
    tmp = Region2_MaxInt(:,:,s);
    Region2_MaxInt_col{s} = tmp(:);
    Region2_MaxInt_col{s}(isnan(Region2_MaxInt_col{s}))=[];

end


