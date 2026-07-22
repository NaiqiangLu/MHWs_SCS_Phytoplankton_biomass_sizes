function [MHW,mclim,m90,mhw_ts,category_ts]=detect(temp,time,cli_start,cli_end,mhw_start,mhw_end,varargin)
%detect - Detecting spatial MHW/MCS  
%  Syntax
%
%  [MHW]=detect(temp,time,cli_start,cli_end,mhw_start,mhw_end)
%  [MHW,mclim,m90,mhw_ts]=detect(temp,time,cli_start,cli_end,mhw_start,mhw_end);
%  [MHW,mclim,m90,mhw_ts]=detect(temp,time,cli_start,cli_end,mhw_start,mhw_end,'Event','MCS','Threshold',0.1);
%
%  Description
%
%  [MHW]=detect(temp,time,cli_start,cli_end,mhw_start,mhw_end) returns
%  all detected MHW events for the m-by-n-by-t matrix TEMP starting in the
%  year DATA_START. m, n and t separately indicate two spatial dimensions
%  (m and n) and one temporal dimension (t). Climatologies used to
%  determine events are calculated based on TEMP from CLI_START to
%  CLI_END. MHW is a table where each row corresponds to a particular
%  event during MHW_START to MHW_END and each column indicates a metric.
%
%  [MHW,mclim,m90,mhw_ts]=detect(temp,time,cli_start,cli_end,mhw_start,mhw_end)
%  also return the spatial climatology MCLIM (m-by-n-by-366) and threshold
%  M90 (m-by-n-by-366) to calculate MHW events and resultant MHW time series
%  MHW_TS (m-by-n-by-t). In the condition that there is no missing value in
%  data TEMP, NaN in all outputs indicates lands and 0 in MHW_TS indicates
%  the corresponding day in that grid is not in a MHW event.
%
%  [MHW,mclim,m90,mhw_ts]=detect(temp,time,cli_start,cli_end,mhw_start,mhw_end,'Event','MCS','Threshold',0.1)
%  returns MCS events based on 10th percentile threshold.
%
%  Input Arguments
%
%   temp - 3D daily temperature to detect MHW/MCS events, specified as a
%   m-by-n-by-t matrix. m and n separately indicate two spatial dimensions
%   and t indicates temporal dimension. 
%
%   time - datenum(start_year,start_month,start_day):datenum(end_year,
%   end_month,end_day)
%
%   cli_start - A numeric value in format of datennum(yyyy,mm,dd), indicating the start date for the period
%   across which the spatial climatology and threshold are calculated. 
%
%   cli_end - A numeric value in format of datennum(yyyy,mm,dd) indicating the end year for the period across
%   which the spatial climatology and threshold are calculated. 
%
%   data_start - A numeric value in format of datennum(yyyy,mm,dd) indicating the start year of your input
%   data TEMP.
%
%   mhw_start - A numeric value in format of datennum(yyyy,mm,dd) indicating the start year for the period
%   across which MHW/MCS events are detected. 
%
%   mhw_end - A numeric value in format of datennum(yyyy,mm,dd) indicating the end year for the period across
%   which MHW/MCS events are detected.
%
%   'Event' - Default is 'MHW'.
%           - 'MHW' - detecting MHW events.
%           - 'MCS' - detecting MCS events.
%
%   'Threshold' - Default is 0.9. Threshold percentile to detect MHW/MCS
%   events.
%
%   'windowHalfWidth' - Default is 5. Width of sliding window to calculate
%   spatial climatology and threshold. 
%
%   'smoothPercentileWidth' - Default is 31. Width of moving mean window to smooth spatial
%   climatology and threshold.
%
%   'minDuration' - Default is 5. Minimum duration to accept a detection of MHW/MCS
%   event.
%
%   'maxGap' - Default is 2. Maximum gap accepting joining of MHW events.
%   
%   'ClimTemp' - Default is TEMP. The data used to calculate climatology
%   and thresholds.
%
%   'ClimTime' - A vector of datenum() corresponding to ClimTemp.
%   
%   'percentile' - Default is 'matlab'. Indicating the way to calculate
%   percentile, using either 'matlab way' or 'python way'.
%
%  Output Arguments
%   
%   MHW - A table containing all detected MHW/MCS events where each row
%   corresponds to a particular event and each column corresponds to a
%   metric. Specified metrics are:
%       - mhw_onset - onset date of each event.
%       - mhw_end - end date of each event.
%       - mhw_dur - duration of each event.
%       - int_max - maximum intensity of each event.
%       - int_mean - mean intensity of each event.
%       - int_var - variance of intensity during each event.
%       - int_cum - cumulative intensity across each event.
%       - xloc - location of each event in x-dimension of TEMP.
%       - yloc - location of each event in y-dimension of TEMP. 
%
%   mclim - A 3D matrix (m-by-n-by-366) containing climatologies.
%
%   m90 - A 3D matrix (m-by-n-by-366) containing thresholds.
%
%   mhw_ts - A 3D matrix
%   (m-by-n-by-(datenum(MHW_end,1,1)-datenum(MHW_start)+1)) containing 
%   spatial intensity of MHW/MCS in each day.
%
%   category_ts - A 3D matrix
%   (m-by-n-by-(datenum(MHW_end,1,1)-datenum(MHW_start)+1)) containing 
%   category of MHW/MCS in each day.


% vEvent = 'MHW';
% vThreshold = 0.9;
% vWindowHalfWidth = 5;
% vsmoothPercentileWidth = 31;
% vminDuration = 5;
% vmaxGap = 2;
% ClimTemp = temp;
% ClimTime = time;
% 
paramNames = {'Event','Threshold','WindowHalfWidth','smoothPercentileWidth','minDuration',...
    'maxGap','ClimTemp','ClimTime','percentile'};
defaults   = {'MHW',0.9,5,31,5,2,temp,time,'matlab'};
% 输入参数有事件类型MHW，阈值0.9，平均窗口期5，滑动平均窗口期31，温度3d矩阵数据输入，对应的温度3d矩阵数据时间跨度，百分位计算方式matlab环境

% varargin = reshape(varargin,2,length(varargin)/2);
% 
% for i = 1:length(defaults)
%     if any(ismember(varargin(1,:),paramNames{i}))
%        feval(@()assignin('caller',paramNames{i},varargin{2,ismember(varargin(1,:),paramNames{i})}))
%     end      
% end  
%这段代码的核心是解析键值对可变参数，并将参数值动态赋值到调用者工作区；

[vEvent, vThreshold,vWindowHalfWidth,vsmoothPercentileWidth,vminDuration,vmaxGap,ClimTemp,ClimTime,vpercentile]...
    = internal.stats.parseArgs(paramNames, defaults, varargin{:});

EventNames = {'MHW','MCS'};
vEvent = internal.stats.getParamVal(vEvent,EventNames,...
    '''Event''');  %解析事件类型  
PercentileNames={'matlab','python'};
vpercentile=internal.stats.getParamVal(vpercentile,PercentileNames,...
    '''matlab''');  %百分位算法检查

%%  "What if cli_start-window or cli_end+window exceeds the time range of data"
%举例：气候态基期为1982.1.1-2005.12.31   海温数据时间是1982.1.1-2016.12.31
ahead_date=ClimTime(1)-(cli_start-vWindowHalfWidth); %an举例:1982.1.1-(1982.1.1-5)=5  % >0表示在气候态窗在前端的日期越界了，超出了我们输入的海温数据范围，<=0则表示气候态窗口没有越界
after_date=cli_end+vWindowHalfWidth-ClimTime(end);   %2005.12.31+5-2016.12.31<0       %与上面一样的道理，都表示所超出来的天数
temp_clim=ClimTemp(:,:,ClimTime>=cli_start-vWindowHalfWidth & ClimTime<=cli_end+vWindowHalfWidth); %在给出气候态基期后，我们下一步要在气候态基期
% 最左右两边提取气候态数据段，我们的气候态基期温度数据要选取到数据源最左右两边的各加多5天时间。


if ahead_date>0 && after_date>0  %气候态窗口在前后两端都超出了海温数据源的日期界限，对超出部分的天数数据进行赋nan填充，补齐理论窗口长度
    temp_clim=cat(3,NaN(size(temp_clim,1),size(temp_clim,2),ahead_date), ...
    temp_clim,NaN(size(temp_clim,1),size(temp_clim,2),after_date));  
elseif ahead_date>0 && after_date<=0
    temp_clim=cat(3,NaN(size(temp_clim,1),size(temp_clim,2),ahead_date), ...
    temp_clim);  %只在前端超过日期界限，只补充前端窗口  %这里要注意cat函数的用法，在变量temp_clim的前后页拼接是有区别的
elseif ahead_date<=0 && after_date>0
        temp_clim=cat(3, ...
            temp_clim,NaN(size(temp_clim,1),size(temp_clim,2),after_date));  %在后端超过界限，只补充后端的窗口
else
    %这种情况是窗口不越界，不用补齐数据窗口长度
end

temp_mhw=temp(:,:,time>=mhw_start & time<=mhw_end);  %只提取用于选定的MHW检测的时间段的海温数据

%% Calculating climatology and thresholds

date_true=datevec(cli_start-vWindowHalfWidth:cli_end+vWindowHalfWidth);  %气候态基期的时间轴
date_true=date_true(:,1:3);  %只取年月日前三个要素

date_false = date_true;
date_false(:,1) = 2012;  %取2012年作为闰年 有366天

fake_doy = day(datetime(date_false),'dayofyear');  %day函数的用法 'dayofyear'，提取在当年的第几天  拿2012这个闰年给一年中的第几天做doy，任何闰年都可以，目的是索引这些年月在闰年里所处的doy
ind = 1:length(date_false);  %时间索引（1,2,3,…）对于气候基期的时间，我们索引他排在第几位（天）

mclim=NaN(size(temp,1),size(temp,2),366); 
m90=NaN(size(temp,1),size(temp,2),366);    %初始化输出矩阵，对于每一个格点都有一个气候态基期的平均值和百分位阈值

for i=1:366
    if i == 60    %因为是366天，所以这里的i=60是指2月29日 一年中的第60天  这里先不算2月29日 因为样本太少后面单独处理
         
    else
        ind_fake=ind;
        ind_fake(fake_doy==i & ~ismember(datenum(date_true),cli_start:cli_end))=nan; %首先判断第i天在8776天的一年中的天数是否相等，就是判断第i天在8776已经排好的doy是不是一致，然后再得到判断超出气候基期的列向量，当两者条件都为1（ture）时，表示在第i天，但不在气候基期里面，赋值为nan。
    data_thre=num2cell(temp_clim(:,:,any(ind_fake'>=(ind_fake(fake_doy == i)-vWindowHalfWidth) & ind_fake' <= (ind_fake(fake_doy ==i)+vWindowHalfWidth),2)),3);
    switch vpercentile
        case 'matlab'
            m90(:,:,i) = quantile(temp_clim(:,:,any(ind_fake'>=(ind_fake(fake_doy == i)-vWindowHalfWidth) & ind_fake' <= (ind_fake(fake_doy ==i)+vWindowHalfWidth),2)),vThreshold,3);
            mclim(:,:,i) = mean(temp_clim(:,:,any(ind_fake'>=(ind_fake(fake_doy == i)-vWindowHalfWidth) & ind_fake' <= (ind_fake(fake_doy ==i)+vWindowHalfWidth),2)),3,'omitnan');
            %这里就是先判断8776已经排好的doy的第i天是否相等，判断完后第i天对应的doy则有一个逻辑列向量，然后再被索引为8776doy的排位（时间点），然后再被索引为这些时间点在窗口内的doy范围，这里ind_fake'判断隐含了ind_fake'：所有时间点（Ntime
            %× 1），lower/upper（就是左右窗口）：多个年份中第 i 天 ±window 的索引（1 ×
            %Nyear），产生[Ntime × Nyear]
            %的逻辑矩阵，每一行：某一个具体时间点；每一列：某一年的窗口，然后用any（:,:,2）表明每一个时间点落进某一年的窗口里面，这个时间点就要。
        case 'python'
            
            m90(:,:,i) = cellfun(@percentile, data_thre,repmat({vThreshold},size(temp,1),size(temp,2)));
            mclim(:,:,i) = mean(temp_clim(:,:,any(ind_fake'>=(ind_fake(fake_doy == i)-vWindowHalfWidth) & ind_fake' <= (ind_fake(fake_doy ==i)+vWindowHalfWidth),2)),3,'omitnan');
    end
    
    end
end
% Dealing with Feb29
m90(:,:,60) = mean(m90(:,:,[59 61]),3,'omitnan');
mclim(:,:,60) = mean(mclim(:,:,[59 61]),3,'omitnan');

% Does running averages of threshold and clim..做滑动平均处理

m90long=smoothdata(cat(3,m90,m90,m90),3,'movmean',vsmoothPercentileWidth);  %沿第三维拼接3次使得用90百分位阈值数据首尾相连，来做31天的滑动平均。31天滑动平均是指对每一个像元 (x,y)，对每一个像元 (x,y)，来平滑每天的 m90。
m90=m90long(:,:,367:367+365); %m90取m90long的中间的一年366天，因为中间这一年首尾都有数据，滑动平均时数据不会缺少，最合理。
mclimlong=smoothdata(cat(3,mclim,mclim,mclim),3,'movmean',vsmoothPercentileWidth);
mclim=mclimlong(:,:,367:367+365); %和上面同理。总而言之，取前后5天作为每个格点的窗口数据是为了增大样本量，做年内31天滑动平均是为了90%分位阈值和气候态均值能保持年内的连续性、平滑性。

[x_size,y_size]=deal(size(m90,1),size(m90,2));

mbigadd=temp_mhw;  %把要检测的热浪数据复制一份

date_mhw=datevec(mhw_start:mhw_end);  %热浪检测时间范围
date_mhw(:,1)=2000;   %把检测年份全部改为2000年  取闰年的doy排序，和上面取2012年的doy做法一样
indextocal = day(datetime(date_mhw),'dayofyear'); %把检测热浪的时间也按照闰年的doy排序来索引  

ts=str2double(string(datestr(mhw_start:mhw_end,'YYYYmmdd')));  %把热浪检测时间从datenum数值转换为19930101-20161231

mhw_ts=NaN(x_size,y_size,length(ts)); %创建一个空的矩阵，热浪检测时间的每一天，都有热浪检测区域的网格矩阵， [x_size, y_size, 时间长度]，保存每天每个格点的热浪造成的异常温度值，nan值表示陆地、非MHW天、缺测；正值表示MHW 日的强度。
category_ts=NaN(x_size,y_size,length(ts)); %创建空的热浪等级时间序列

MHW=[];  %初始化一个空的 MHW 事件结构体数组

%% Detecting MHW/MCS in each grid

switch vEvent
    case 'MHW'
        
        for i=1:x_size
            for j=1:y_size   %热浪检测区域的网格i*j
                
                mhw_ts(i,j,isnan(squeeze(mbigadd(i,j,:))))=nan; %第一件事：处理缺测点（NaN），如果这个格点某天 SST 是 NaN，那一天 不可能是 MHW，对应的 mhw_ts 也设为 NaN
                
                if sum(isnan(squeeze(mbigadd(i,j,:))))~=size(mbigadd,3) %识别不是陆地的区域 原理是判断这个格点不是“整条时间序列全是 NaN”，那就是海洋了，下面才做计算
                    
                    maysum=zeros(size(mbigadd,3),1);  %对于每个格点，0表示没有超出阈值，1表示超出阈值，对于每个格点我们先创建一个时间维度上全0的矩阵，然后再检测90%阈值，超出阈值则赋值1
                    
                    maysum(squeeze(mbigadd(i,j,:))>=squeeze(m90(i,j,indextocal)))=1;  %检测mhw是否超过阈值，对于每一个格点，对应的所在doy，是否超出了对应m90的doy的90%分位阈值
                    
                    trigger=0;  %trigger = 0：现在不在事件里；trigger = 1：正在一个事件里
                    potential_event=[];
                    
                    for n=1:size(maysum,1)  %情况 1：事件开始   一天一天判断
                        if trigger==0 && maysum(n)==1   %开始判断mhw事件，n=1处判断起始于哪一天
                            start_here=n;
                            trigger=1;  %判断有超过阈值的day，改为1，表示mhw在该事件里
                        elseif trigger==1 && maysum(n)==0   %情况 2：事件结束
                            end_here=n-1;  %这里索引mhw结束的doy  这是正确索引mhw结束的时间在哪里
                            trigger=0;  %表示mhw事件结束，现在是mhw事件后一个格点 无mhw 赋值0  从 1 变成 0 → 事件结束
                            potential_event=[potential_event;[start_here end_here]];  %把mhw所处时间存起来
                        elseif n==size(maysum,1) && trigger==1 && maysum(n)==1  %情况 3：一直到最后一天还在事件里，补上最后一个事件
                            trigger=0;
                            end_here=n;
                            potential_event=[potential_event;[start_here end_here]];
                        end
                    end
                    
                    if ~isempty(potential_event) %首先判断非空的mhw事件集合才能进行下一步计算
                        
                        potential_event=potential_event((potential_event(:,2)-potential_event(:,1)+1)>=vminDuration,:); %这里筛选>=5天的时间的事件，然后符合要求的被列出来
                         
                        if ~isempty(potential_event)
                            
                            gaps=(potential_event(2:end,1) - potential_event(1:(end-1),2) - 1); %把发生mhw事件之间的事件间隔算出来；gap = 下一个事件开始 - 上一个事件结束 - 1得出不同事件之间隔的天数

                            
                            while min(gaps)<=vmaxGap  %判断这些mhw事件的间隔如果小于等于两天，则合并为一个事件
                                %                                  potential_event(find(gaps<=vmaxGap),2)=potential_event(find(gaps<=vmaxGap)+1,2);
                                %                                  potential_event(find(gaps<=vmaxGap)+1,:)=[];
                                %                                  gaps=(potential_event(2:end,1) - potential_event(1:(end-1),2) - 1);
                                
                                potential_event(find(gaps<=vmaxGap),2)=potential_event(find(gaps<=vmaxGap)+1,2); %筛选间隔小于2天的事件，把前一个事件结束的时间替换为后一个事件结束的时间
                                loc_should_del=(find(gaps<=vmaxGap)+1); %标记要删除的行  这些行是已经被合并进前一个事件的“后续事件”
                                loc_should_del=loc_should_del(~ismember(loc_should_del,find(gaps<=vmaxGap)));  %这里是为了避免删除错误的行。如果出现了某些连续的gap<=2的情况下，如果后两个符合gap条件的事件都被合并前一个事件当中去，这行的代码就保证了只删除mhw事件合并后最右边的代表的行的mhw事件。
                                potential_event(loc_should_del,:)=[];  %真正删除的行，把已经“合并进去”的事件删掉
                                gaps=(potential_event(2:end,1) - potential_event(1:(end-1),2) - 1);  %对已经合并完的事件重新计算gaps；因为事件数量、顺序都变了，gaps 必须重新算。直到min(gaps) > vmaxGap，所有事件之间都间隔足够远，不能再合并了。
                            end
                            %为每一个最终热浪事件，准备好“结果存储位置”，每一行 = 一场热浪。
                            mhwstart=NaN(size(potential_event,1),1);
                            mhwend=NaN(size(potential_event,1),1);
                            mduration=NaN(size(potential_event,1),1);
                            mhwint_max=NaN(size(potential_event,1),1);
                            mhwint_mean=NaN(size(potential_event,1),1);
                            mhwint_var=NaN(size(potential_event,1),1);  %强度方差
                            mhwint_cum=NaN(size(potential_event,1),1);
                            mhwcategory=NaN(size(potential_event,1),1); %热浪事件等级 （hobday 2018）
                            
                            for le=1:size(potential_event,1)  %已经确定好的最终所有mhw事件
                                event_here=potential_event(le,:);  %把每个mhw事件提取出来
                                endtime=ts(event_here(2));  %热浪结束时间
                                starttime=ts(event_here(1)); %热浪开始时间  时间格式为：YYYYMMDD
                                mcl=squeeze(mclim(i,j,indextocal(event_here(1):event_here(2))));
                                mth=squeeze(m90(i,j,indextocal(event_here(1):event_here(2))));  %mhw事件期间取出来每个网格的90%分位阈值和气候态均值，哪一天到哪一天。
                                mrow=squeeze(mbigadd(i,j,event_here(1):event_here(2)));  %把热浪期间的数据提取出来
                                manom=mrow-mcl;  %计算温度异常
                                mca=mth-mcl;     %90%分位阈值减去气候态均值温度  阈值距平
                                mhw_ts(i,j,event_here(1):event_here(2))=manom;  %保存每个网格的每日温度异常数值
                                 
                                [maxanom,~]=nanmax(squeeze(manom));  % 计算该海洋热浪事件期间的最大温度异常（忽略 NaN）
                                
                                mhwint_max(le)=...
                                    maxanom;
                                mhwint_mean(le)=...
                                    mean(manom);
                                mhwint_var(le)=...
                                    std(manom);   %热浪异常温度的波动性
                                mhwint_cum(le)=...
                                    sum(manom);   %累计热暴露
                                mhwstart(le)=starttime;
                                mhwend(le)=endtime;
                                mduration(le)=event_here(2)-event_here(1)+1;
                                mhwcategory(le)=nanmax(floor(manom./mca));  %hobday2018将90%分位-气候态均值作为基准差值（Δ），将mhw温度异常值除以基准差值（Δ）然后向负方向取整得到mhw类型。
                                category_ts(i,j,event_here(1):event_here(2))=nanmax(floor(manom./mca));
                            end
                            mhwcategory(mhwcategory>4)=4;
                            MHW=[MHW;[mhwstart mhwend mduration mhwint_max mhwint_mean mhwint_var mhwint_cum repmat(i,size(mhwstart,1),1) repmat(j,size(mhwstart,1),1) mhwcategory]];
                        end   %mhw事件输出，一行表示一个mhw事件。把当前网格 (i,j) 中检测到的所有海洋热浪事件，按照“开始时间、结束时间、持续时间、强度指标、空间位置和等级”的格式，逐行存入总的 MHW 事件表格中。
                    end
                end
            end
        end
        
    case 'MCS'
        
        for i=1:x_size
            for j=1:y_size
                
                mhw_ts(i,j,isnan(squeeze(mbigadd(i,j,:))))=nan;
                
                if sum(isnan(squeeze(mbigadd(i,j,:))))~=size(mbigadd,3)
                    
                    maysum=zeros(size(mbigadd,3),1);
                    
                    maysum(squeeze(mbigadd(i,j,:))<=squeeze(m90(i,j,indextocal)))=1;
                    
                    trigger=0;
                    potential_event=[];
                    
                    for n=1:size(maysum,1)
                        if trigger==0 && maysum(n)==1
                            start_here=n;
                            trigger=1;
                        elseif trigger==1 && maysum(n)==0
                            end_here=n-1;
                            trigger=0;
                            potential_event=[potential_event;[start_here end_here]];
                        elseif n==size(maysum,1) && trigger==1 && maysum(n)==1
                            trigger=0;
                            end_here=n;
                            potential_event=[potential_event;[start_here end_here]];
                        end
                    end
                    
                    if ~isempty(potential_event)
                        
                        potential_event=potential_event((potential_event(:,2)-potential_event(:,1)+1)>=vminDuration,:);
                        
                        if ~isempty(potential_event)
                            
                            gaps=(potential_event(2:end,1) - potential_event(1:(end-1),2) - 1);
                            
                            while min(gaps)<=vmaxGap
                                %                                  potential_event(find(gaps<=vmaxGap),2)=potential_event(find(gaps<=vmaxGap)+1,2);
                                %                                  potential_event(find(gaps<=vmaxGap)+1,:)=[];
                                %                                  gaps=(potential_event(2:end,1) - potential_event(1:(end-1),2) - 1);
                                
                                potential_event(find(gaps<=vmaxGap),2)=potential_event(find(gaps<=vmaxGap)+1,2);
                                loc_should_del=(find(gaps<=vmaxGap)+1);
                                loc_should_del=loc_should_del(~ismember(loc_should_del,find(gaps<=vmaxGap)));
                                potential_event(loc_should_del,:)=[];
                                gaps=(potential_event(2:end,1) - potential_event(1:(end-1),2) - 1);
                            end
                            
                            mhwstart=NaN(size(potential_event,1),1);
                            mhwend=NaN(size(potential_event,1),1);
                            mduration=NaN(size(potential_event,1),1);
                            mhwint_max=NaN(size(potential_event,1),1);
                            mhwint_mean=NaN(size(potential_event,1),1);
                            mhwint_var=NaN(size(potential_event,1),1);
                            mhwint_cum=NaN(size(potential_event,1),1);
                            mhwcategory=NaN(size(potential_event,1),1);
                            
                            for le=1:size(potential_event,1)
                                event_here=potential_event(le,:);
                                endtime=ts(event_here(2));
                                starttime=ts(event_here(1));
                                mcl=squeeze(mclim(i,j,indextocal(event_here(1):event_here(2))));
                                mth=squeeze(m90(i,j,indextocal(event_here(1):event_here(2))));
                                mrow=squeeze(mbigadd(i,j,event_here(1):event_here(2)));
                                manom=mrow-mcl;
                                mca=mth-mcl;
                                mhw_ts(i,j,event_here(1):event_here(2))=manom;
                                
                                [maxanom,~]=nanmin(squeeze(manom));
                                
                                mhwint_max(le)=...
                                    maxanom;
                                mhwint_mean(le)=...
                                    mean(manom);
                                mhwint_var(le)=...
                                    std(manom);
                                mhwint_cum(le)=...
                                    sum(manom);
                                mhwstart(le)=starttime;
                                mhwend(le)=endtime;
                                mduration(le)=event_here(2)-event_here(1)+1;
                                mhwcategory(le)=nanmax(floor(manom./mca));
                                category_ts(i,j,event_here(1):event_here(2))=nanmax(floor(manom./mca));
                            end
                            
                            mhwcategory(mhwcategory>4)=4;
                            MHW=[MHW;[mhwstart mhwend mduration mhwint_max mhwint_mean mhwint_var mhwint_cum repmat(i,size(mhwstart,1),1) repmat(j,size(mhwstart,1),1) mhwcategory]];
                            
                        end
                    end
                end
            end
        end
end

category_ts(category_ts>4)=4;
MHW=table(MHW(:,1),MHW(:,2),MHW(:,3),MHW(:,4),MHW(:,5),MHW(:,6),MHW(:,7),MHW(:,8),MHW(:,9),MHW(:,10),...
    'variablenames',{'mhw_onset','mhw_end','mhw_dur','int_max','int_mean','int_var','int_cum','xloc','yloc','category'});

function p=percentile(data,thre)
if nansum(isnan(data))~=length(data)
    data=data(~isnan(data));
    pos=1+(length(data)-1)*thre;
    p=interp1((1:length(data))',sort(data(:)),pos);
else
    p=nan;
end