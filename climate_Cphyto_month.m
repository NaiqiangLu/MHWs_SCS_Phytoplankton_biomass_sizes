%% Seasonal mean log10(Cphyto) maps
clear;clc;close all

load('climCphyto.mat');

%% ===== 366天daily climatology -> 12个月 climatology =====

month_days=[31 29 31 30 31 30 31 31 30 31 30 31];

clim_month=nan(size(Cphyto_clim,1),...
    size(Cphyto_clim,2),12);

st=1;

for m=1:12

    ed=st+month_days(m)-1;

    clim_month(:,:,m)=mean(...
        Cphyto_clim(:,:,st:ed),...
        3,'omitnan');

    st=ed+1;

end

%% ===== 夏季(JJA) & 冬季(DJF) =====

summer=mean(clim_month(:,:,[6 7 8]),...
    3,'omitnan');

winter=mean(clim_month(:,:,[12 1 2]),...
    3,'omitnan');

%% ===== log10 =====

summer(summer<=0)=nan;
winter(winter<=0)=nan;

summer_log=log10(summer);
winter_log=log10(winter);

plot_data={summer_log,winter_log};

titles={'Summer (JJA)','Winter (DJF)'};

%% ===== Figure =====

figure('Position',[100 200 900 800],...
    'Color','w');

tl=tiledlayout(1,2,...
    'TileSpacing','loose',...
    'Padding','loose');

for k=1:2

    nexttile

    %% 投影
    m_proj('equidistant',...
        'lon',[105 123],...
        'lat',[0 25]);

    %% 主图
    m_contourf(Lon,Lat,...
        plot_data{k},...
        80,...
        'linestyle','none');

    shading interp
    hold on

    %% 海岸线
    m_gshhs_i('patch',[0.8 0.8 0.8]);
    m_gshhs_i('linewidth',1.2,...
        'color','k');

    %% 网格
    m_grid('linestyle','none',...
        'box','on',...
        'linewidth',1.5,...
        'fontsize',16,...
        'fontname','times');

    %% colormap
    colormap(nclCM(172,14))

    caxis([0.4   1.2])

    %% 区域框

    if k==1        % Summer (JJA)

        % WSCS框（红色）
        lon_box=[107.5 113.5 113.5 107.5 107.5];
        lat_box=[8.5 8.5 13.5 13.5 8.5];

        m_line(lon_box,lat_box,...
            'color','red',...
            'linewidth',2);

    elseif k==2    % Winter (DJF)

        % Luzon框（黑色）
        lon_box=[117 121 121 117 117];
        lat_box=[16.5 16.5 21 21 16.5];

        m_line(lon_box,lat_box,...
            'color','k',...
            'linewidth',2);
    end

end

%% ===== 手动放置共用colorbar =====

cb=colorbar('southoutside');

set(cb,...
    'Position',[0.25 0.10 0.50 0.03],...
    'tickdir','in',...
    'linewidth',1.2,...
    'fontsize',16,...
    'edgecolor','k',...
    'FontName','Times New Roman');

xlabel(cb,...
    'Cphyto (log_{10}(mg m^{-3}))',...
    'fontsize',18,'FontName','Times New Roman');

%% 保存

print('-dtiff','-r300',...
    'Seasonal_log10_Cphyto')

disp('Seasonal Cphyto figure saved.')