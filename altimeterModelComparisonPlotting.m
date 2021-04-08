% altModelComparisonPlotting
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Script to make plots for comparing WW3 output to altimeter measurements %
% of significant wave height and wind speed                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Author:       Clarence Olin Collins III, Ph.D.                          %
% Affiliation:  ERDC - FRF                                                %
% created:      02/21/2020                                                %
% version:      1.0                                                       %
% updates:      wind 04/05/202                                            %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%% filter out bad model data
[j , k] = find(pData.mdHs < 0 | pData.altHs < 0.75);
keepIndex = setdiff(1:length(pData.time),j);
fields = fieldnames(pData);
for jj = 1:length(fields)
    if strcmp(fields{jj},'meta') | strcmp(fields{jj},'options')
        continue
    end
    temp = pData.(fields{jj});
    pData.(fields{jj}) = temp(keepIndex) ;
end

%% covert longitude to -180:180 for Atlantic
pData.lon = deg180(pData.lon);
% or do 360 for Pacific
% pData.lon = deg360(pData.lon);

%% check figure directory
[glyph, ~] = giveGlyph;
direct = dir(figurePath);
if isempty(direct)
    mkdir(figurePath);
end
%% time series plot
lonLatTimeSeries(pData.lon,pData.lat,pData.time,pData.altHs,pData.mdHs)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
plotName = 'timeSeries';
plotPath = [figurePath glyph plotName];
savefig(plotPath)
export_fig(plotPath,'-transparent')
%% calculate global stats
globalStats = Ocstatsp(pData.altHs,pData.mdHs,0.1)
% save('globalStats','globalStats')
%%
% distributions, histograms
%% q-q plots
qqplot(pData.altHs,pData.mdHs)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
plotName = 'q-q';
plotPath = [figurePath glyph plotName];
savefig(plotPath)
export_fig(plotPath,'-transparent')
%% scatter plot
colorScatterQQ(pData.altHs,pData.mdHs,100)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
plotName = 'colorScatter';
plotPath = [figurePath glyph plotName];
savefig(plotPath)
export_fig(plotPath,'-transparent')
%% taylor diagram

%% map of stats
clear statGridCenterLat statGridCenterLon statGridEdgesLat statGridEdgesLon...
    bias rmse nbias nrmse scatIn meanAlt meanMod dHs

minLON = min(pData.lon);
maxLON = max(pData.lon);
minLAT = min(pData.lat);
maxLAT = max(pData.lat);
blockSize = 1; %[degrees]
statGridEdgesLon = minLON:blockSize:maxLON;
statGridEdgesLat = minLAT:blockSize:maxLAT;
statGridCenterLon = statGridEdgesLon(1:end-1) + diff(statGridEdgesLon)/2;
statGridCenterLat = statGridEdgesLat(1:end-1) + diff(statGridEdgesLat)/2;

for i = 1:length(statGridCenterLon)
    for j = 1:length(statGridCenterLat)
        [index, ~] = find(pData.lon >= statGridEdgesLon(i) & pData.lon <= statGridEdgesLon(i+1) &...
            pData.lat >= statGridEdgesLat(j) & pData.lat <= statGridEdgesLat(j+1));
        meanAlt(i,j) = nanmean(pData.altHs(index));
        meanMod(i,j) = nanmean(pData.mdHs(index));
        dWIND = pData.altHs(index)-pData.mdHs(index);
        bias(i,j) = meanMod(i,j) - meanAlt(i,j);
        nbias(i,j) = bias(i,j)./meanAlt(i,j);
        rmse(i,j) = sqrt(mean(dWIND.^2));
        nrmse(i,j) = sqrt(mean(dWIND.^2)/meanAlt(i,j).^2);
        scatIn(i,j) = rmse(i,j)/meanAlt(i,j);
    end
end
%%
% bHS = pData.altHs-pData.mdHs;
% rHS = sqrt((pData.altHs-pData.mdHs).^2);

% figure
% m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
% m_pcolor(deg180(mdCol.lon),mdCol.lat,mean(mdCol.hs,3)')
% shading('interp')
% m_coast('patch',[.8 .8 .8]);
% m_grid('box','fancy','tickdir','out');
% title('mean model wave height (native res.)')
% fontsize(16,16,16,16)
% h = colorbar;
% ylabel(h, '[m]')
% colormap('lansey')
% caxis([0 10])

figure
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,meanMod')
shading('interp')
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
title('mean model wave height')
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('lansey')
caxis([0 10])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
plotName = 'meanModelMap';
plotPath = [figurePath glyph plotName];
savefig(plotPath)
export_fig(plotPath,'-transparent')

figure
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,meanAlt')
shading('interp')
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
title('mean altimeter wave height')
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('lansey')
caxis([0 10])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
plotName = 'meanAltimeterMap';
plotPath = [figurePath glyph plotName];
savefig(plotPath)
export_fig(plotPath,'-transparent')

figure
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,rmse')
shading('interp')
hold on
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
title('rmse')
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('jet')
caxis([0 1])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
plotName = 'rmseMap';
plotPath = [figurePath glyph plotName];
savefig(plotPath)
export_fig(plotPath,'-transparent')

figure
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,bias')
shading('interp')
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
title('bias')
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('lansey')
caxis([-1 1])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
plotName = 'biasMap';
plotPath = [figurePath glyph plotName];
savefig(plotPath)
export_fig(plotPath,'-transparent')

figure
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,nrmse')
shading('interp')
hold on
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
title('nrmse')
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('jet')
caxis([0 1])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
plotName = 'normalizedRmseMap';
plotPath = [figurePath glyph plotName];
savefig(plotPath)
export_fig(plotPath,'-transparent')

figure
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,scatIn')
shading('interp')
hold on
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
title('scatIn')
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('jet')
caxis([0 1])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
plotName = 'scatterIndexMap';
plotPath = [figurePath glyph plotName];
savefig(plotPath)
export_fig(plotPath,'-transparent')

figure
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,nbias')
shading('interp')
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
title('nbias')
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('lansey')
caxis([-1 1])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
plotName = 'normalizedBiasMap';
plotPath = [figurePath glyph plotName];
savefig(plotPath)
export_fig(plotPath,'-transparent')



%% months and seasons
if options.seasonalPlots == 1
    date = datevec(pData.time);
    
    for month = 1:12
        monthIndex = date(:,2) == month;
        eval(['LON' num2str(month) '= pData.lon(monthIndex);'])
        eval(['LAT' num2str(month) '= pData.lat(monthIndex);'])
        eval(['TIME' num2str(month) '= pData.time(monthIndex);'])
        eval(['HSobs' num2str(month) '= pData.altHs(monthIndex);'])
        eval(['HSmd' num2str(month) '= pData.mdHs(monthIndex);'])
    end
    
    for i = 1:length(statGridCenterLon)
        for j = 1:length(statGridCenterLat)
            for month = 1:12
                LON = eval(['LON' num2str(month)]);
                LAT = eval(['LAT' num2str(month)]);
                HSobs = eval(['HSobs' num2str(month)]);
                HSmd = eval(['HSmd' num2str(month)]);
                [index , ~] = find(LON >= statGridEdgesLon(i) & LON <= statGridEdgesLon(i+1) &...
                    LAT >= statGridEdgesLat(j) & LAT <= statGridEdgesLat(j+1));
                dWIND = HSmd(index)-HSobs(index);
                eval(['bias' num2str(month) '(i,j) = mean(dHS);'])
                eval(['rmse' num2str(month) '(i,j) = sqrt(mean(dHS.^2));'])
            end
        end
    end
    
    LONs1 = [LON3; LON4; LON5];
    LONs2 = [LON6; LON7; LON8];
    LONs3 = [LON9; LON10; LON11];
    LONs4 = [LON12; LON1; LON2];
    LATs1 = [LAT3; LAT4; LAT5];
    LATs2 = [LAT6; LAT7; LAT8];
    LATs3 = [LAT9; LAT10; LAT11];
    LATs4 = [LAT12; LAT1; LAT2];
    HSobss1 = [HSobs3; HSobs4; HSobs5];
    HSobss2 = [HSobs6; HSobs7; HSobs8];
    HSobss3 = [HSobs9; HSobs10; HSobs11];
    HSobss4 = [HSobs12; HSobs1; HSobs2];
    HSmds1 = [HSmd3; HSmd4; HSmd5];
    HSmds2 = [HSmd6; HSmd7; HSmd8];
    HSmds3 = [HSmd9; HSmd10; HSmd11];
    HSmds4 = [HSmd12; HSmd1; HSmd2];
    
    for i = 1:length(statGridCenterLon)
        for j = 1:length(statGridCenterLat)
            for s = 1:4
                LON = eval(['LONs' num2str(s)]);
                LAT = eval(['LATs' num2str(s)]);
                HSobs = eval(['HSobss' num2str(s)]);
                HSmd = eval(['HSmds' num2str(s)]);
                [index , ~] = find(LON >= statGridEdgesLon(i) & LON <= statGridEdgesLon(i+1) &...
                    LAT >= statGridEdgesLat(j) & LAT <= statGridEdgesLat(j+1));
                dWIND = HSmd(index)-HSobs(index);
                eval(['biasS' num2str(s) '(i,j) = mean(dHS);'])
                eval(['rmseS' num2str(s) '(i,j) = sqrt(mean(dHS.^2));'])
            end
        end
    end
    %
    figure
    subplot(2,2,1)
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_pcolor(statGridCenterLon,statGridCenterLat,rmseS1')
    shading('interp')
    m_coast('patch',[.8 .8 .8]);
    m_grid('box','fancy','tickdir','out');
    title('RMSE - Spring')
    fontsize(16,16,16,16)
    h = colorbar;
    ylabel(h, '[m]')
    colormap('jet')
    caxis([0 1])
    subplot(2,2,2)
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_pcolor(statGridCenterLon,statGridCenterLat,rmseS2')
    shading('interp')
    m_coast('patch',[.8 .8 .8]);
    m_grid('box','fancy','tickdir','out');
    title('Summer')
    fontsize(16,16,16,16)
    h = colorbar;
    ylabel(h, '[m]')
    colormap('jet')
    caxis([0 1])
    subplot(2,2,3)
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_pcolor(statGridCenterLon,statGridCenterLat,rmseS3')
    shading('interp')
    m_coast('patch',[.8 .8 .8]);
    m_grid('box','fancy','tickdir','out');
    title('Fall')
    fontsize(16,16,16,16)
    h = colorbar;
    ylabel(h, '[m]')
    colormap('jet')
    caxis([0 1])
    subplot(2,2,4)
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_pcolor(statGridCenterLon,statGridCenterLat,rmseS4')
    shading('interp')
    m_coast('patch',[.8 .8 .8]);
    m_grid('box','fancy','tickdir','out');
    title('Winter')
    fontsize(16,16,16,16)
    h = colorbar;
    ylabel(h, '[m]')
    colormap('jet')
    caxis([0 1])
    packfig(2,2)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    plotName = 'seasonalRmseMap';
    plotPath = [figurePath glyph plotName];
    savefig(plotPath)
    export_fig(plotPath,'-transparent')
    
    figure
    subplot(2,2,1)
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_pcolor(statGridCenterLon,statGridCenterLat,biasS1')
    shading('interp')
    m_coast('patch',[.8 .8 .8]);
    m_grid('box','fancy','tickdir','out');
    title('bias - Spring')
    fontsize(16,16,16,16)
    h = colorbar;
    ylabel(h, '[m]')
    colormap('lansey')
    caxis([-1 1])
    subplot(2,2,2)
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_pcolor(statGridCenterLon,statGridCenterLat,biasS2')
    shading('interp')
    m_coast('patch',[.8 .8 .8]);
    m_grid('box','fancy','tickdir','out');
    title('Summer')
    fontsize(16,16,16,16)
    h = colorbar;
    ylabel(h, '[m]')
    colormap('lansey')
    caxis([-1 1])
    subplot(2,2,3)
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_pcolor(statGridCenterLon,statGridCenterLat,biasS3')
    shading('interp')
    m_coast('patch',[.8 .8 .8]);
    m_grid('box','fancy','tickdir','out');
    title('Fall')
    fontsize(16,16,16,16)
    h = colorbar;
    ylabel(h, '[m]')
    colormap('lansey')
    caxis([-1 1])
    subplot(2,2,4)
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_pcolor(statGridCenterLon,statGridCenterLat,biasS4')
    shading('interp')
    m_coast('patch',[.8 .8 .8]);
    m_grid('box','fancy','tickdir','out');
    title('Winter')
    fontsize(16,16,16,16)
    h = colorbar;
    ylabel(h, '[m]')
    colormap('lansey')
    caxis([-1 1])
    packfig(2,2)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    plotName = 'seasonalBiasMap';
    plotPath = [figurePath glyph plotName];
    savefig(plotPath)
    export_fig(plotPath,'-transparent')
    
end
%% WIND
if pData.meta.options.runWindEval == 1    
    lonLatTimeSeries(pData.lon,pData.lat,pData.time,pData.altWind,pData.mdWind)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    plotName = 'timeSeriesWind';
    plotPath = [figurePath glyph plotName];
    savefig(plotPath)
    export_fig(plotPath,'-transparent')
    %% calculate global stats
    globalStats = Ocstatsp(pData.altWind,pData.mdWind,0.1)
    % save('globalStats','globalStats')
    %%
    % distributions, histograms
    %% q-q plots
    qqplot(pData.altWind,pData.mdWind)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    plotName = 'q-qWind';
    plotPath = [figurePath glyph plotName];
    savefig(plotPath)
    export_fig(plotPath,'-transparent')
    %% scatter plot
    colorScatter(pData.altWind,pData.mdWind,100)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    plotName = 'colorScatterWind';
    plotPath = [figurePath glyph plotName];
    savefig(plotPath)
    export_fig(plotPath,'-transparent')
    %% taylor diagram
    
    %% map of stats
    clear statGridCenterLat statGridCenterLon statGridEdgesLat statGridEdgesLon...
        bias rmse nbias nrmse scatIn meanAlt meanMod dHs
    
    minLON = min(pData.lon);
    maxLON = max(pData.lon);
    minLAT = min(pData.lat);
    maxLAT = max(pData.lat);
    blockSize = 1; %[degrees]
    statGridEdgesLon = minLON:blockSize:maxLON;
    statGridEdgesLat = minLAT:blockSize:maxLAT;
    statGridCenterLon = statGridEdgesLon(1:end-1) + diff(statGridEdgesLon)/2;
    statGridCenterLat = statGridEdgesLat(1:end-1) + diff(statGridEdgesLat)/2;
    
    for i = 1:length(statGridCenterLon)
        for j = 1:length(statGridCenterLat)
            [index dumy] = find(pData.lon >= statGridEdgesLon(i) & pData.lon <= statGridEdgesLon(i+1) &...
                pData.lat >= statGridEdgesLat(j) & pData.lat <= statGridEdgesLat(j+1));
            meanAlt(i,j) = nanmean(pData.altWind(index));
            meanMod(i,j) = nanmean(pData.mdWind(index));
            dWIND = pData.altWind(index)-pData.mdWind(index);
            bias(i,j) = meanMod(i,j) - meanAlt(i,j);
            nbias(i,j) = bias(i,j)./meanAlt(i,j);
            rmse(i,j) = sqrt(mean(dWIND.^2));
            nrmse(i,j) = sqrt(mean(dWIND.^2)/meanAlt(i,j).^2);
            scatIn(i,j) = rmse(i,j)/meanAlt(i,j);
        end
    end
    %%
    % bHS = pData.altHs-pData.mdHs;
    % rHS = sqrt((pData.altHs-pData.mdHs).^2);
    
    % figure
    % m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    % m_pcolor(deg180(mdCol.lon),mdCol.lat,mean(mdCol.hs,3)')
    % shading('interp')
    % m_coast('patch',[.8 .8 .8]);
    % m_grid('box','fancy','tickdir','out');
    % title('mean model wave height (native res.)')
    % fontsize(16,16,16,16)
    % h = colorbar;
    % ylabel(h, '[m]')
    % colormap('lansey')
    % caxis([0 10])
    
    figure
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_pcolor(statGridCenterLon,statGridCenterLat,meanMod')
    shading('interp')
    m_coast('patch',[.8 .8 .8]);
    m_grid('box','fancy','tickdir','out');
    title('mean model wind speed')
    fontsize(16,16,16,16)
    h = colorbar;
    ylabel(h, '[m/s]')
    colormap('lansey')
    caxis([0 20])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    plotName = 'meanModelMapWind';
    plotPath = [figurePath glyph plotName];
    savefig(plotPath)
    export_fig(plotPath,'-transparent')
    
    figure
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_pcolor(statGridCenterLon,statGridCenterLat,meanAlt')
    shading('interp')
    m_coast('patch',[.8 .8 .8]);
    m_grid('box','fancy','tickdir','out');
    title('mean altimeter wind speed')
    fontsize(16,16,16,16)
    h = colorbar;
    ylabel(h, '[m/s]')
    colormap('lansey')
    caxis([0 20])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    plotName = 'meanAltimeterMapWind';
    plotPath = [figurePath glyph plotName];
    savefig(plotPath)
    export_fig(plotPath,'-transparent')
    
    figure
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_pcolor(statGridCenterLon,statGridCenterLat,rmse')
    shading('interp')
    hold on
    m_coast('patch',[.8 .8 .8]);
    m_grid('box','fancy','tickdir','out');
    title('rmse')
    fontsize(16,16,16,16)
    h = colorbar;
    ylabel(h, '[m/s]')
    colormap('jet')
    caxis([0 5])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    plotName = 'rmseMapWind';
    plotPath = [figurePath glyph plotName];
    savefig(plotPath)
    export_fig(plotPath,'-transparent')
        
    figure
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_pcolor(statGridCenterLon,statGridCenterLat,100*nrmse')
    shading('interp')
    hold on
    m_coast('patch',[.8 .8 .8]);
    m_grid('box','fancy','tickdir','out');
    title('nrmse')
    fontsize(16,16,16,16)
    h = colorbar;
    ylabel(h, '[%]')
    colormap('jet')
    caxis([0 100])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    plotName = 'normalizedRmseMapWind';
    plotPath = [figurePath glyph plotName];
    savefig(plotPath)
    export_fig(plotPath,'-transparent')
    
    figure
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_pcolor(statGridCenterLon,statGridCenterLat,scatIn')
    shading('interp')
    hold on
    m_coast('patch',[.8 .8 .8]);
    m_grid('box','fancy','tickdir','out');
    title('scatIn')
    fontsize(16,16,16,16)
    h = colorbar;
    ylabel(h, '[]')
    colormap('jet')
    caxis([0 1])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    plotName = 'scatterIndexMapWind';
    plotPath = [figurePath glyph plotName];
    savefig(plotPath)
    export_fig(plotPath,'-transparent')
    
    figure
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_pcolor(statGridCenterLon,statGridCenterLat,100*nbias')
    shading('interp')
    m_coast('patch',[.8 .8 .8]);
    m_grid('box','fancy','tickdir','out');
    title('nbias')
    fontsize(16,16,16,16)
    h = colorbar;
    ylabel(h, '[%]')
    colormap('lansey')
    caxis([-25 25])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    plotName = 'normalizedBiasMapWind';
    plotPath = [figurePath glyph plotName];
    savefig(plotPath)
    export_fig(plotPath,'-transparent')

    figure
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_pcolor(statGridCenterLon,statGridCenterLat,bias')
    shading('interp')
    m_coast('patch',[.8 .8 .8]);
    m_grid('box','fancy','tickdir','out');
    title('nbias')
    fontsize(16,16,16,16)
    h = colorbar;
    ylabel(h, '[m/s]')
    colormap('lansey')
    caxis([-5 5])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    plotName = 'biasMapWind';
    plotPath = [figurePath glyph plotName];
    savefig(plotPath)
    export_fig(plotPath,'-transparent')
    
    
    %% months and seasons
    if options.seasonalPlots == 1
        date = datevec(pData.time);
        
        for month = 1:12
            monthIndex = date(:,2) == month;
            eval(['LON' num2str(month) '= pData.lon(monthIndex);'])
            eval(['LAT' num2str(month) '= pData.lat(monthIndex);'])
            eval(['TIME' num2str(month) '= pData.time(monthIndex);'])
            eval(['WINDobs' num2str(month) '= pData.altWind(monthIndex);'])
            eval(['WINDmd' num2str(month) '= pData.mdWind(monthIndex);'])
        end
        
        for i = 1:length(statGridCenterLon)
            for j = 1:length(statGridCenterLat)
                for month = 1:12
                    LON = eval(['LON' num2str(month)]);
                    LAT = eval(['LAT' num2str(month)]);
                    WINDobs = eval(['WINDobs' num2str(month)]);
                    WINDmd = eval(['WINDmd' num2str(month)]);
                    [index , ~] = find(LON >= statGridEdgesLon(i) & LON <= statGridEdgesLon(i+1) &...
                        LAT >= statGridEdgesLat(j) & LAT <= statGridEdgesLat(j+1));
                    dWIND = WINDmd(index)-WINDobs(index);
                    eval(['bias' num2str(month) '(i,j) = mean(dWIND);'])
                    eval(['rmse' num2str(month) '(i,j) = sqrt(mean(dWIND.^2));'])
                end
            end
        end
        
        LONs1 = [LON3; LON4; LON5];
        LONs2 = [LON6; LON7; LON8];
        LONs3 = [LON9; LON10; LON11];
        LONs4 = [LON12; LON1; LON2];
        LATs1 = [LAT3; LAT4; LAT5];
        LATs2 = [LAT6; LAT7; LAT8];
        LATs3 = [LAT9; LAT10; LAT11];
        LATs4 = [LAT12; LAT1; LAT2];
        WINDobss1 = [WINDobs3; WINDobs4; WINDobs5];
        WINDobss2 = [WINDobs6; WINDobs7; WINDobs8];
        WINDobss3 = [WINDobs9; WINDobs10; WINDobs11];
        WINDobss4 = [WINDobs12; WINDobs1; WINDobs2];
        WINDmds1 = [WINDmd3; WINDmd4; WINDmd5];
        WINDmds2 = [WINDmd6; WINDmd7; WINDmd8];
        WINDmds3 = [WINDmd9; WINDmd10; WINDmd11];
        WINDmds4 = [WINDmd12; WINDmd1; WINDmd2];
        
        for i = 1:length(statGridCenterLon)
            for j = 1:length(statGridCenterLat)
                for s = 1:4
                    LON = eval(['LONs' num2str(s)]);
                    LAT = eval(['LATs' num2str(s)]);
                    WINDobs = eval(['WINDobss' num2str(s)]);
                    WINDmd = eval(['WINDmds' num2str(s)]);
                    [index , ~] = find(LON >= statGridEdgesLon(i) & LON <= statGridEdgesLon(i+1) &...
                        LAT >= statGridEdgesLat(j) & LAT <= statGridEdgesLat(j+1));
                    dWIND = WINDmd(index)-WINDobs(index);
                    eval(['biasS' num2str(s) '(i,j) = mean(dWIND);'])
                    eval(['rmseS' num2str(s) '(i,j) = sqrt(mean(dWIND.^2));'])
                end
            end
        end
        %
        figure
        subplot(2,2,1)
        m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
        m_pcolor(statGridCenterLon,statGridCenterLat,rmseS1')
        shading('interp')
        m_coast('patch',[.8 .8 .8]);
        m_grid('box','fancy','tickdir','out');
        title('RMSE - Spring')
        fontsize(16,16,16,16)
        h = colorbar;
        ylabel(h, '[m/s]')
        colormap('jet')
        caxis([0 1])
        subplot(2,2,2)
        m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
        m_pcolor(statGridCenterLon,statGridCenterLat,rmseS2')
        shading('interp')
        m_coast('patch',[.8 .8 .8]);
        m_grid('box','fancy','tickdir','out');
        title('Summer')
        fontsize(16,16,16,16)
        h = colorbar;
        ylabel(h, '[m/s]')
        colormap('jet')
        caxis([0 1])
        subplot(2,2,3)
        m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
        m_pcolor(statGridCenterLon,statGridCenterLat,rmseS3')
        shading('interp')
        m_coast('patch',[.8 .8 .8]);
        m_grid('box','fancy','tickdir','out');
        title('Fall')
        fontsize(16,16,16,16)
        h = colorbar;
        ylabel(h, '[m/s]')
        colormap('jet')
        caxis([0 1])
        subplot(2,2,4)
        m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
        m_pcolor(statGridCenterLon,statGridCenterLat,rmseS4')
        shading('interp')
        m_coast('patch',[.8 .8 .8]);
        m_grid('box','fancy','tickdir','out');
        title('Winter')
        fontsize(16,16,16,16)
        h = colorbar;
        ylabel(h, '[m/s]')
        colormap('jet')
        caxis([0 1])
        packfig(2,2)
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        plotName = 'seasonalRmseMapWind';
        plotPath = [figurePath glyph plotName];
        savefig(plotPath)
        export_fig(plotPath,'-transparent')
        
        figure
        subplot(2,2,1)
        m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
        m_pcolor(statGridCenterLon,statGridCenterLat,biasS1')
        shading('interp')
        m_coast('patch',[.8 .8 .8]);
        m_grid('box','fancy','tickdir','out');
        title('bias - Spring')
        fontsize(16,16,16,16)
        h = colorbar;
        ylabel(h, '[m/s]')
        colormap('lansey')
        caxis([-1 1])
        subplot(2,2,2)
        m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
        m_pcolor(statGridCenterLon,statGridCenterLat,biasS2')
        shading('interp')
        m_coast('patch',[.8 .8 .8]);
        m_grid('box','fancy','tickdir','out');
        title('Summer')
        fontsize(16,16,16,16)
        h = colorbar;
        ylabel(h, '[m/s]')
        colormap('lansey')
        caxis([-1 1])
        subplot(2,2,3)
        m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
        m_pcolor(statGridCenterLon,statGridCenterLat,biasS3')
        shading('interp')
        m_coast('patch',[.8 .8 .8]);
        m_grid('box','fancy','tickdir','out');
        title('Fall')
        fontsize(16,16,16,16)
        h = colorbar;
        ylabel(h, '[m/s]')
        colormap('lansey')
        caxis([-1 1])
        subplot(2,2,4)
        m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
        m_pcolor(statGridCenterLon,statGridCenterLat,biasS4')
        shading('interp')
        m_coast('patch',[.8 .8 .8]);
        m_grid('box','fancy','tickdir','out');
        title('Winter')
        fontsize(16,16,16,16)
        h = colorbar;
        ylabel(h, '[m/s]')
        colormap('lansey')
        caxis([-1 1])
        packfig(2,2)
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        plotName = 'seasonalBiasMapWind';
        plotPath = [figurePath glyph plotName];
        savefig(plotPath)
        export_fig(plotPath,'-transparent')
        
    end
end