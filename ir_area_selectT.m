function [sst,sstir] = ir_area_selectT(sst_full,lon,lat,waterline,lon_lat,lola)
%  Syntax
%  获取缩小区域减水深数据(sst)，以及不规则区域的数据(sstir)
%  ！！！只适用于珠江口海表温度！！！
%  需要自定义缩小区域
%
%  Input Arguments：
% 
%  SST_full - 处理的数据 sst_full(lon,lat,time)
%   
%  lon - 对应的经度 100 101 102……
%   
%  lat - 对应的纬度 0 1 2 3……
%   
%  waterline - 去除的水深 eg. waterline = -10；
%   
%  lon_lat - 不规则区域的顶点 eg. lon_lat = [111 21;112 20;118 22;117 23];
%
%  lola - 预处理 缩小的区域 eg. lola = [113,119,20.5,23.5];
% 
%  Output Arguments：
% 
%  sst - 缩小区域的Chla减去对应水深
%   
%  sstir - 缩小区域 减去 对应水深的 "不规则区域"的数据
%===========================================================%


x1 = find(lon>=lola(1) & lon<=lola(2)); y1 = find(lat>=lola(3) & lat<=lola(4));
lon=lon(x1);lat=lat(y1);SST = sst_full(x1,y1,:);

% 水深
lo = ncread('topo_20.1.nc','lon');
la = ncread('topo_20.1.nc','lat');
depth = ncread('topo_20.1.nc','z');

x3=find(lo>=lola(1) & lo<=lola(2)); y3=find(la>=lola(3) & la<=lola(4));
la = la(y3); lo = lo(x3); dep = depth(x3,y3);

% 去除水深
x=1:length(la);y=1:length(lo);
[X,Y]=meshgrid(x,y);
xi=1:(length(la)/length(lat)):length(la);
yi=1:(length(lo)/length(lon)):length(lo);
[Xi,Yi]=meshgrid(xi,yi);
z1=interp2(X,Y,dep,Xi,Yi,'spline');
z1(z1>=waterline)=nan;z1(z1<=waterline)=0;

sst = SST+z1;

% 选取不规则区域
[Lo,La]=meshgrid(lon,lat);
loc_in = inpolygon(Lo,La,lon_lat(:,1),lon_lat(:,2));  %把落到不规则多边形的网格圈出来

sstir=[];
for d1 = 1:size(sst,3)
    B = sst(:,:,d1);
    sstir(:,d1) = B(loc_in');
end

end