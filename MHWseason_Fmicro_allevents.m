%% 提取 Region1(夏季) 和 Region2(冬季) 每次MHW事件的Fmicro例异常
clear;clc;close all;

%% ---------------- Fmicro ----------------
load Fphyto1998-2024.mat

Fmicro(Fmicro< 0)=NaN;
Fmicro = permute(Fmicro,[2 1 3]);

load ClimFphyto1998-2024.mat

climFmicro = permute(climFmicro,[2 1 3]);

studytime = datenum(1998,1,1):datenum(2024,12,31);

period_plot_v = datevec(studytime);
period_unique = datevec(datenum(2016,1,1):datenum(2016,12,31));

[~,loc_plot] = ismember(period_plot_v(:,2:3),...
                        period_unique(:,2:3),'rows');

mFmicro = squeeze(climFmicro(:,:,loc_plot));

%% ---------------- SST & MHW ----------------
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

%% ---------------- Region ----------------

x1 = find(Lon>=109 & Lon<=111);
y1 = find(Lat>=10 & Lat<=12);

x2 = find(Lon>=118.5 & Lon<=121);
y2 = find(Lat>=18 & Lat<=20);

%% ---------------- 保存结果 ----------------

Region1_FMicro_Event = [];
Region2_FMicro_Event = [];

%% ==========================================================
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

        %% 当前网格所有MHW事件
        rfmicro = nan(size(period_mhw,1),1);

        for loc = 1:size(period_mhw,1)

            mhw_time = period_mhw(loc,1):period_mhw(loc,2);

            tidx = mhw_time - studytime(1) + 1;

            mhw_fmicro = squeeze(...
                Fmicro(loc_here(1),loc_here(2),tidx));

            clim_fmicro = squeeze(...
                mFmicro(loc_here(1),loc_here(2),tidx));

            %% 当前MHW事件Fpico比例异常（百分比）
            rfmicro(loc) = nanmean((mhw_fmicro-clim_fmicro).*100);

        end

        %% =====================================================
        %% Region1（夏季）

        if seasonal==2

            if ismember(loc_here(1),x1) && ismember(loc_here(2),y1)

                Region1_FMicro_Event = ...
                    [Region1_FMicro_Event; rfmicro(:)];

            end

        end

        %% =====================================================
        %% Region2（冬季）

        if seasonal==4

            if ismember(loc_here(1),x2) && ismember(loc_here(2),y2)

                Region2_FMicro_Event = ...
                    [Region2_FMicro_Event; rfmicro(:)];

            end

        end

    end

end

