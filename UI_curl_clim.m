clc; clear; close all;

data_path = 'E:\第二篇\wind1998-2024';

load(fullfile(data_path,'UI&curl1998-2024.mat'));
% 变量：ui & curl (100×100×9862)

UI_VNclim = daily_climatology(UI_VN, ...
    datetime(1998,1,1), ...
    datetime(2024,12,31));


UI_LZclim = daily_climatology(UI_LZ, ...
    datetime(1998,1,1), ...
    datetime(2024,12,31));


Curl_clim=daily_climatology(curl, ...
    datetime(1998,1,1), ...
    datetime(2024,12,31));


save(fullfile(data_path,'UI_Curl_clim1998_2024.mat'), ...
    'UI_VNclim','UI_LZclim','Curl_clim','Lon','Lat','-v7.3');