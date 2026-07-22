%% ============================================================
%  Extract Wind Speed anomaly for every MHW event
%  Region1 : Summer (109-111E,10-12N)
%  Region2 : Winter (118.5-121E,18-20N)
%
%  Output:
%  Region1_WindEvent
%  Region2_WindEvent
%
%  每一个元素对应一次MHW事件期间平均风速异常(m/s)
% =============================================================

clear;clc;close all

%% Wind

load CCMP_wind1998_2024.mat
load climWind_1998_2024.mat

studytime = datenum(1998,1,1):datenum(2024,12,31);

period_plot_v = datevec(studytime);

period_unique = datevec(datenum(2016,1,1):datenum(2016,12,31));

[~,loc_plot] = ismember(period_plot_v(:,2:3),...
                        period_unique(:,2:3),'rows');

mWindS = squeeze(ws_clim(:,:,loc_plot));

%% MHW

load sst_19822025.mat

[MHW,smclim,m90,mhw_ts] = detect(...
    sst_full,...
    datenum(1982,1,1):datenum(2025,12,31),...
    datenum(1995,1,1),...
    datenum(2024,12,31),...
    datenum(1998,1,1),...
    datenum(2024,12,31),...
    'Threshold',0.9);

mhw = MHW{:,:};

season = [3 4 5;
          6 7 8;
          9 10 11;
          12 1 2];

%% Region

load sst_Lon_Lat.mat

x1 = find(Lon>=109 & Lon<=111);
y1 = find(Lat>=10  & Lat<=12);

x2 = find(Lon>=118.5 & Lon<=121);
y2 = find(Lat>=18    & Lat<=20);

%% 保存所有事件

Region1_WindEvent = [];

Region2_WindEvent = [];

%% ============================================================
% Four seasons
%% ============================================================

for seasonal = 1:4

    loc_plot = MHWs_Season_judge(mhw,season(seasonal,:));

    mhw_season = mhw(loc_plot,:);

    loc_full = unique(mhw_season(:,8:9),'rows');

    for m = 1:size(loc_full,1)

        loc_here = loc_full(m,:);

        mhw_here = mhw_season(...
            mhw_season(:,8)==loc_here(1) & ...
            mhw_season(:,9)==loc_here(2),:);

        period_mhw = [...
            datenum(num2str(mhw_here(:,1)),'yyyymmdd'),...
            datenum(num2str(mhw_here(:,2)),'yyyymmdd')];

        valwinds = nan(size(period_mhw,1),1);

        for loc = 1:size(period_mhw,1)

            mhw_time = period_mhw(loc,1):period_mhw(loc,2);

            mhw_ws = squeeze(...
                ws(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));

            clim_ws = squeeze(...
                mWindS(loc_here(1),loc_here(2),...
                mhw_time-datenum(1998,1,1)+1));

            mhw_ws(mhw_ws<0)=NaN;
            clim_ws(clim_ws<0)=NaN;

            % 风速异常值 (m/s)

            valwinds(loc)=nanmean(mhw_ws-clim_ws);

        end

        %% Region1 Summer

        if seasonal==2

            if ismember(loc_here(1),x1) && ismember(loc_here(2),y1)

                Region1_WindEvent = ...
                    [Region1_WindEvent;
                     valwinds(:)];

            end

        end

        %% Region2 Winter

        if seasonal==4

            if ismember(loc_here(1),x2) && ismember(loc_here(2),y2)

                Region2_WindEvent = ...
                    [Region2_WindEvent;
                     valwinds(:)];

            end

        end

    end

end

%% 去NaN

Region1_WindEvent(isnan(Region1_WindEvent))=[];

Region2_WindEvent(isnan(Region2_WindEvent))=[];

%% 查看数量

length(Region1_WindEvent)

length(Region2_WindEvent)