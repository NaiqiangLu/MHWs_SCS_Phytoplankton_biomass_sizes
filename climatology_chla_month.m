%% Seasonal climatology log10(Chla) maps
clear;clc;close all

load('Chla_clim.mat');   % climchla [100×100×366]

%% ===== 1. 366天daily climatology -> 12个月 climatology =====

month_days = [31 29 31 30 31 30 31 31 30 31 30 31];

clim_month = nan(size(climchla,1),size(climchla,2),12);

st = 1;

for m = 1:12

    ed = st + month_days(m)-1;

    clim_month(:,:,m) = mean(climchla(:,:,st:ed),3,'omitnan');

    st = ed+1;

end

%% ===== log10 Chla =====

clim_month(clim_month<=0)=nan;

logchl = log10(clim_month);

%% ===== 北半球季节 =====

season_months = {
    [12 1 2];      % DJF
    [3 4 5];       % MAM
    [6 7 8];       % JJA
    [9 10 11]};    % SON

season_names = {'DJF','MAM','JJA','SON'};

month_names = {'Jan','Feb','Mar','Apr','May','Jun',...
               'Jul','Aug','Sep','Oct','Nov','Dec'};

%% ===== colorbar范围（建议值）=====
% Chla一般log10范围

cmin = -1.5;
cmax = 0.8;

%% ===== 开始绘图 =====

for s = 1:4

    months = season_months{s};

    figure('Position',[100 200 1150 450],...
        'Color','w');

    tiledlayout(1,3,...
        'TileSpacing','compact',...
        'Padding','compact');

    for k = 1:3

        m = months(k);

        nexttile

        %% 投影
        m_proj('equidistant',...
            'lon',[105 123],...
            'lat',[0 25]);

        %% 主图
        m_contourf(Lon,Lat,...
            logchl(:,:,m),...
            80,...
            'linestyle','none');

        shading interp
        hold on

        %% 海岸线
        m_gshhs_i('patch',[0.8 0.8 0.8]);

        m_gshhs_i('color','k',...
            'linewidth',1.2);

        %% 网格
        m_grid('linestyle','none',...
            'box','on',...
            'linewidth',1.2,...
            'fontsize',12,...
            'fontname','times');

        %% colormap
        colormap(nclCM(160))

        caxis([cmin cmax])

        %% ===== WSCS区域框 =====
        lon_box = [107.5 113.5 113.5 107.5 107.5];
        lat_box = [8.5 8.5 13.5 13.5 8.5];

        m_line(lon_box,lat_box,...
            'color','red',...
            'linewidth',1.3);

        %% ===== Luzon区域框 =====
        lon_box = [117 121 121 117 117];
        lat_box = [16.5 16.5 21 21 16.5];

        m_line(lon_box,lat_box,...
            'color','k',...
            'linewidth',1.3);

        %% 标题
        title(month_names{m},...
            'fontsize',16,...
            'fontweight','bold',...
            'FontName','Times New Roman')

        %% ===== 仅第三张图加colorbar =====
        if k==3

            cb=colorbar('eastoutside');

            set(cb,...
                'tickdir','in',...
                'linewidth',1.2,...
                'fontsize',15,...
                'edgecolor','k',...
                'FontName','Times New Roman');

            ylabel(cb,...
                'log_{10}(Chl-a) (mg m^{-3})',...
                'fontsize',16,...
                'fontweight','bold',...
                'FontName','Times New Roman');

        end
    end

    %% 保存
    print('-dtiff','-r300',...
        ['Seasonal_logChla_' season_names{s}])

end

disp('All seasonal climatology figures saved.')