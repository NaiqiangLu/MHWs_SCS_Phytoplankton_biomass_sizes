clc; clear; close all;

data_path = 'E:\第二篇\1998-2024月(daily)平均热通量\GPCPday 降水1998-2024';

load(fullfile(data_path,'Rain_19982024_interp025.mat'));
% 变量：降水(100×100×9862)

Precip_clim = daily_climatology(Rain, ...
                datetime(1998,1,1), ...
                datetime(2024,12,31));


save(fullfile(data_path,'Clim_precipitation.mat'), ...
     'Precip_clim','Lon','Lat','-v7.3');