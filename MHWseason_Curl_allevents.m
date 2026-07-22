%% ==========================================================
% Extract every MHW event Curl anomaly
% Region1 (Summer)
% Region2 (Winter)
%% ==========================================================

clear;clc;close all;

%% ---------------- Load data ----------------

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

mCurl = Curl_clim(:,:,loc_plot);

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

Region1_CurlEvent = [];

Region2_CurlEvent = [];

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

            vcurl = nan(size(period_mhw,1),1);

            for loc = 1:size(period_mhw,1)

                mhw_time = period_mhw(loc,1):period_mhw(loc,2);

                mhw_curl = squeeze(...
                    curl(loc_here(1),loc_here(2),...
                    mhw_time-datenum(1998,1,1)+1));

                clim_curl = squeeze(...
                    mCurl(loc_here(1),loc_here(2),...
                    mhw_time-datenum(1998,1,1)+1));

                vcurl(loc) = nanmean(mhw_curl-clim_curl);

            end

            Region1_CurlEvent = [Region1_CurlEvent; vcurl];

        end

        %% ---------------- Region2 Winter ----------------

        if seasonal==4 && ...
                ismember(loc_here(1),x2) && ...
                ismember(loc_here(2),y2)

            vcurl = nan(size(period_mhw,1),1);

            for loc = 1:size(period_mhw,1)

                mhw_time = period_mhw(loc,1):period_mhw(loc,2);

                mhw_curl = squeeze(...
                    curl(loc_here(1),loc_here(2),...
                    mhw_time-datenum(1998,1,1)+1));

                clim_curl = squeeze(...
                    mCurl(loc_here(1),loc_here(2),...
                    mhw_time-datenum(1998,1,1)+1));

                vcurl(loc) = nanmean(mhw_curl-clim_curl);

            end

            Region2_CurlEvent = [Region2_CurlEvent; vcurl];

        end

    end

end

%% Remove NaN

Region1_CurlEvent(isnan(Region1_CurlEvent)) = [];

Region2_CurlEvent(isnan(Region2_CurlEvent)) = [];

%% Unit conversion (与空间图一致)

Region1_CurlEvent = Region1_CurlEvent ./ 1e-6;

Region2_CurlEvent = Region2_CurlEvent ./ 1e-6;


