clc,clear;close all;
load CCMP_wind1998_2024.mat;

Clim_ws=mclim(ws,datenum(1998,1,1):datenum(2024,12,31),...
    datenum(1998,1,1),datenum(2024,12,31),31);  %数据:起止时间，气候基期，滑动窗口长度; 默认的窗口半宽为5


save('mhw_climway_WS.mat',"Clim_ws","Lon","Lat");




