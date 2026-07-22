clc; clear; close all;

data_path = 'G:\第二篇\ML_irradiance';

load(fullfile(data_path,'Ibar1998_2024.mat'));


Ibarclim = daily_climatology(Ibar,datetime(1998,1,1),datetime(2024,12,31));


save(fullfile(data_path,'Ibarclim_1998_2024.mat'), ...
    'Ibarclim','Lon','Lat','-v7.3');