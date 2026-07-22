clc; clear; close all;

data_path = 'E:\第二篇\ML_irradiance\混合层深度daily_1998_2024';

load(fullfile(data_path,'mld1998_2024_0.25.mat'));
% 变量：MLD (100×100×9862)

MLDclim = daily_climatology(MLD, ...
                datetime(1998,1,1), ...
                datetime(2024,12,31));

save(fullfile(data_path,'clim_MLD1998_2024.mat'), ...
     'MLDclim','Lon','Lat','-v7.3');
