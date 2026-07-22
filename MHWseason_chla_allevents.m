%% =====================================================
% Extract Chla anomaly for every MHW event
% Region1 : Summer
% Region2 : Winter
%% =====================================================

clear;clc;

%% ---------------- Chla ----------------

load Chla_19982024_interp025.mat

chla_full(chla_full< 0)=nan;

Chla = permute(chla_full,[2 1 3]);

load Chla_clim.mat

climchla = permute(climchla,[2 1 3]);

studytime = datenum(1998,1,1):datenum(2024,12,31);

period_plot_v = datevec(studytime);

period_unique = datevec(datenum(2016,1,1):datenum(2016,12,31));

[~,loc_plot] = ismember(period_plot_v(:,2:3),...
                        period_unique(:,2:3),'rows');

mChla = squeeze(climchla(:,:,loc_plot));

%% ---------------- SST ----------------

load sst_19822025.mat

[MHW,~,~,~] = detect(...
    sst_full,...
    datenum(1982,1,1):datenum(2025,12,31),...
    datenum(1995,1,1),...
    datenum(2024,12,31),...
    datenum(1998,1,1),...
    datenum(2024,12,31));

mhw = MHW{:,:};

%% ---------------- Region ----------------

load sst_Lon_Lat.mat

x1 = find(Lon>=109 & Lon<=111);
y1 = find(Lat>=10  & Lat<=12);

x2 = find(Lon>=118.5 & Lon<=121);
y2 = find(Lat>=18    & Lat<=20);

%% =====================================================
% Region1  Summer


loc_summer = MHWs_Season_judge(mhw,[6 7 8]);

mhw_summer = mhw(loc_summer,:);

Region1_ChlaEvent = [];

for k = 1:size(mhw_summer,1)

    ix = mhw_summer(k,8);
    iy = mhw_summer(k,9);

    if ~ismember(ix,x1) || ~ismember(iy,y1)
        continue
    end

    t1 = datenum(num2str(mhw_summer(k,1)),'yyyymmdd');
    t2 = datenum(num2str(mhw_summer(k,2)),'yyyymmdd');

    mhw_time = t1:t2;

    mhw_chla = squeeze(...
        Chla(ix,iy,mhw_time-datenum(1998,1,1)+1));

    clim_chla = squeeze(...
        mChla(ix,iy,mhw_time-datenum(1998,1,1)+1));

    mhw_chla(mhw_chla< 0)=NaN;
    clim_chla(clim_chla< 0)=NaN;

    rchla = nanmean(...
        ((mhw_chla-clim_chla)./clim_chla).*100);

    Region1_ChlaEvent(end+1,1)=rchla;

end

%% =====================================================
% Region2 Winter


loc_winter = MHWs_Season_judge(mhw,[12 1 2]);

mhw_winter = mhw(loc_winter,:);

Region2_ChlaEvent = [];

for k = 1:size(mhw_winter,1)

    ix = mhw_winter(k,8);
    iy = mhw_winter(k,9);

    if ~ismember(ix,x2) || ~ismember(iy,y2)
        continue
    end

    t1 = datenum(num2str(mhw_winter(k,1)),'yyyymmdd');
    t2 = datenum(num2str(mhw_winter(k,2)),'yyyymmdd');

    mhw_time = t1:t2;

    mhw_chla = squeeze(...
        Chla(ix,iy,mhw_time-datenum(1998,1,1)+1));

    clim_chla = squeeze(...
        mChla(ix,iy,mhw_time-datenum(1998,1,1)+1));

    mhw_chla(mhw_chla< 0)=NaN;
    clim_chla(clim_chla< 0)=NaN;

    rchla = nanmean(...
        ((mhw_chla-clim_chla)./clim_chla).*100);

    Region2_ChlaEvent(end+1,1)=rchla;

end

%% =====================================================
% check
%% =====================================================

size(Region1_ChlaEvent)
size(Region2_ChlaEvent)