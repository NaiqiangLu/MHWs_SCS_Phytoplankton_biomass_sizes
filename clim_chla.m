clc; clear;close all;

% ===============================
%  Daily climatology of Chl-a
%  Period : 1998–2024
%  Output : 601 × 601 × 366
% ================================

% 数据路径
path_chla = 'E:\第二篇\chl-a_1998_2024\';

years = 1998:2024;

% 网格大小
nx = 601;   %chla_yyyy.mat 内的变量为Chla (601*601*365/366)
ny = 601;  

% 初始化
sumchl  = zeros(nx,ny,366);   % 累加
countchl = zeros(nx,ny,366);  % 有效数据计数

% ===============================
% 主循环

for y = 1:length(years)
    
    year = years(y);
    
    file = [path_chla,'chla_',num2str(year),'.mat'];
    
    load(file,'Chla');  % Chla: 601×601×365/366
    
    ndays = size(Chla,3);
    
    % 计算该年的日期
    dates = datenum(year,1,1):datenum(year,12,31);
    
    for d = 1:ndays
        
        % 当前日期
        date_current = dates(d);
        
        % 计算一年中的第几天 (1–366)
        doy = date_current - datenum(year,1,0);
        
        data = Chla(:,:,d);
        
        % 有效数据
        valid = ~isnan(data);
        
        % NaN -> 0 用于求和
        data(~valid) = 0;
        
        % 累加
        sumchl(:,:,doy) = sumchl(:,:,doy) + data;
        
        % 统计有效像元
        countchl(:,:,doy) = countchl(:,:,doy) + valid;
        
    end
    
    disp(['Finished year ',num2str(year)])

end

% 计算 climatology

climchl = sumchl ./ countchl;

% 无有效值位置设为 NaN
climchl(countchl == 0) = NaN;

% 保存
save('climchl_1998_2024.mat','climchl','-v7.3')

disp('Daily climatology calculation finished!')


%%  画图查看
clc,clear;
load climchl_1998_2024.mat
load chla_lat_lon.mat
climchl60=climchl(:,:,60);

figure('pos',[100 100 1000 1000]);  %绘图框
m_proj('equidistant','lon',[100 125],'lat',[0 25]);%投影范围
m_contourf(lon,lat,climchl60',500,'linestyle','none');
shading interp;
hold on;
m_gshhs_i('linewidth',1,'color','k');%生成地形
m_gshhs_i('patch',[0.7,0.7,0.7]);%陆地颜色灰度
m_grid('linestyle','none','box','on','linewidth',1,'fontsize',26);%生成经纬度网格

hold on
colormap(m_colmap('jet'));
cb=colorbar;
caxis([0  0.8]);
set(cb,'tickdir','in'); % 刻度线朝内
set(cb,'linewidth',1,'fontsize',28,'edgecolor','k','FontName','Times New Roman');%设置colorbar字体大小和框线宽度
set(cb,'YTick',0:0.2:0.8); %色标值范围及显示间隔
set(get(cb,'ylabel'),'string','Chla(​mg/m^{3})','FontSize',30,'FontName','Times New Roman','FontWeight','bold');%生成colorbar标题，单位