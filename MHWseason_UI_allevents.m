%% ==========================================================
% Extract every MHW event UI anomaly
% Region1 (Summer) : UI_VN
% Region2 (Winter) : UI_LZ
%% ==========================================================

clear;clc;close all;

%% ------------------ Load data ------------------

load UI&curl1998-2024.mat
load UI_Curl_clim1998_2024.mat

load sst_19822025.mat
load sst_Lon_Lat.mat

%% Daily climatology

studytime = datenum(1998,1,1):datenum(2024,12,31);

period_plot_v = datevec(studytime);

period_unique = datevec(datenum(2016,1,1):datenum(2016,12,31));

[~,loc_plot] = ismember(period_plot_v(:,2:3),...
                        period_unique(:,2:3),'rows');

mUI_VN = UI_VNclim(:,:,loc_plot);

mUI_LZ = UI_LZclim(:,:,loc_plot);

%% Detect MHW

[MHW,~,~,~] = detect(...
    sst_full,...
    datenum(1982,1,1):datenum(2025,12,31),...
    datenum(1995,1,1),...
    datenum(2024,12,31),...
    datenum(1998,1,1),...
    datenum(2024,12,31),...
    'Threshold',0.9);

mhw = MHW{:,:};

%% ---------------- Region ----------------

x1 = find(Lon>=109 & Lon<=111);
y1 = find(Lat>=10  & Lat<=12);

x2 = find(Lon>=118.5 & Lon<=121);
y2 = find(Lat>=18    & Lat<=20);

%% 保存结果

Region1_UIEvent = [];

Region2_UIEvent = [];

%% Season

season = [3 4 5;
          6 7 8;
          9 10 11;
          12 1 2];

%% =======================================================
% Loop over seasons
%% =======================================================

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

        %% ---------------- Region1 Summer ----------------

        if seasonal==2 && ...
                ismember(loc_here(1),x1) && ...
                ismember(loc_here(2),y1)

            vui_vn = nan(size(period_mhw,1),1);

            for loc=1:size(period_mhw,1)

                mhw_time = period_mhw(loc,1):period_mhw(loc,2);

                mhw_ui = squeeze(UI_VN(...
                    loc_here(1),...
                    loc_here(2),...
                    mhw_time-datenum(1998,1,1)+1));

                clim_ui = squeeze(mUI_VN(...
                    loc_here(1),...
                    loc_here(2),...
                    mhw_time-datenum(1998,1,1)+1));

                vui_vn(loc) = nanmean(mhw_ui-clim_ui);

            end

            Region1_UIEvent = [Region1_UIEvent;vui_vn];

        end

        %% ---------------- Region2 Winter ----------------

        if seasonal==4 && ...
                ismember(loc_here(1),x2) && ...
                ismember(loc_here(2),y2)

            vui_lz = nan(size(period_mhw,1),1);

            for loc=1:size(period_mhw,1)

                mhw_time = period_mhw(loc,1):period_mhw(loc,2);

                mhw_ui = squeeze(UI_LZ(...
                    loc_here(1),...
                    loc_here(2),...
                    mhw_time-datenum(1998,1,1)+1));

                clim_ui = squeeze(mUI_LZ(...
                    loc_here(1),...
                    loc_here(2),...
                    mhw_time-datenum(1998,1,1)+1));

                vui_lz(loc) = nanmean(mhw_ui-clim_ui);

            end

            Region2_UIEvent = [Region2_UIEvent;vui_lz];

        end

    end

end

%% Remove NaN

Region1_UIEvent(isnan(Region1_UIEvent))=[];

Region2_UIEvent(isnan(Region2_UIEvent))=[];