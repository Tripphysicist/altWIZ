%% get all North Atlantic Hurricane Data
clear 
close all
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
stormName = 'all';
stormYear = [ ];
altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/RY19/'; % path to Ribal and Young data
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/altimeterHurricane';
options.altDatabase     = 'RY19'; 
options.maxTimeDiff     =  60/(24*60); % 1 hour
options.maxDistance     =  250; % km
options.QC              =  2; % take coastal data
options.save            =  1;
[stormObs, stormI] = altimeterStormPairing(stormName,stormYear,codePath,altPath,savePath,options)

%%
clear 
close all
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
mdPath = '/Users/tripp/D/Analysis/altimeterComparison/modelData/AtlanticYearRun/';
altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/RY19/'; % path to Ribal and Young data
options.averagingMethod =  1;
options.maxTimeDiff     =  30/(24*60);
options.maxDistance     =  25;
options.minNumberObs    =  7;
options.loopSize        =  3;
options.QC              =  1;
options.save            =  1;
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/Atlantic_2017';
[pData] = altimeterMNodelComparison(codePath,mdPath,altPath,savePath,options)


%%

clear
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
mdPath = '/Users/tripp/D/Analysis/altimeterComparison/modelData/AtlanticYearRun/';
altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/Ribal_Young_2019/'; % path to Ribal and Young data
options.averagingMethod =  2;
options.maxTimeDiff     =  30/(24*60);
options.maxDistance     =  25;
options.minNumberObs    =  7;
options.loopSize        =  10;
options.QC              =  1;
options.save            =  1;
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/Atlantic_2017_Feb_Bub';
[pData] = altimeterWW3Comparison(codePath,mdPath,altPath,savePath,options)
%
clear
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
mdPath = '/Users/tripp/D/Analysis/altimeterComparison/modelData/AtlanticYearRun/';
altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/Ribal_Young_2019/'; % path to Ribal and Young data
options.averagingMethod =  3;
options.maxTimeDiff     =  30/(24*60);
options.maxDistance     =  25;
options.minNumberObs    =  7;
options.loopSize        =  3;
options.QC              =  2;
options.save            =  1;
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/Atlantic_2017_Feb_none_QC2';
[pData] = altimeterWW3Comparison(codePath,mdPath,altPath,savePath,options)
%%
clear
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
mdPath = '/Users/tripp/D/Analysis/altimeterComparison/modelData/AtlanticYearRun/';
altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/Ribal_Young_2019/'; % path to Ribal and Young data
options.averagingMethod =  3;
options.maxTimeDiff     =  30/(24*60);
options.maxDistance     =  25;
options.minNumberObs    =  7;
options.loopSize        =  3;
options.QC              =  1;
options.save            =  0;
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/Atlantic_2017_Feb_none_QC5';
[pData] = altimeterWW3Comparison(codePath,mdPath,altPath,savePath,options)
%% make video
[LO LA] = meshgrid(mdCol.lon,mdCol.lat);
minLON = min(deg180(mdCol.lon));
minLAT = min(mdCol.lat);
maxLON = max(deg180(mdCol.lon));
maxLAT = max(mdCol.lat);

vidObj = VideoWriter('WIS_Atlantic_Feb_2017_color)scale_16_alt.avi');
open(vidObj);
figure(7)
[altInd k] = find(pData.time >= mdCol.time(1) - 1/48 &...
    pData.time < mdCol.time(1) + 1/48);
[aa bb] =  sort(pData.time(altInd));
altLon = pData.lon(altInd(bb));
altLat = pData.lat(altInd(bb));

for i = 1:length(altInd)
h = figure(7);
clf
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_coast('patch',[.8 .8 .8]);
hold on
m_pcolor(deg180(LO),LA,mdCol.hs(:,:,1)')
m_plot(deg180(altLon(i)),altLat(i),'ko')
shg
m_grid('box','fancy','tickdir','out');
h = colorbar;
ylabel(h, '[m]')
colormap('lansey')
fontsize(16,16,16,16)
caxis([0 16])
shg
currFrame = getframe(gcf);
writeVideo(vidObj,currFrame);
end
close(vidObj);    
%% make video
pData.lon = deg180(pData.lon);
minLON = min(pData.lon);
maxLON = max(pData.lon);
minLAT = min(pData.lat);
maxLAT = max(pData.lat);
[LO LA] = meshgrid(mdCol.lonO,mdCol.lat);
vidObj = VideoWriter('WIS_Atlantic_Feb_2017_CS_16_alt_k_keepTrack.avi');
open(vidObj);
figure(7)
date = datevec(mdCol.time);
for i = 1:length(mdCol.time);
[altIndCum k] = find(pData.time >= mdCol.time(1) - 1/48 &...
    pData.time < mdCol.time(i) + 1/48);
[altIndCur k] = find(pData.time >= mdCol.time(i) - 1/48 &...
    pData.time < mdCol.time(i) + 1/48);
h = figure(7);
clf
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_coast('patch',[.8 .8 .8]);
hold on
m_pcolor(deg180(LO),LA,mdCol.hs(:,:,i)')
m_scatter(deg180(pData.lon(altIndCum)),pData.lat(altIndCum),'ko')
m_scatter(deg180(pData.lon(altIndCur)),pData.lat(altIndCur),'go')
shg
m_grid('box','fancy','tickdir','out');
title([num2str(date(i,1)) '-' num2str(date(i,2)) '-' num2str(date(i,3)) ' '...
    num2str(date(i,4)) ':' num2str(date(i,5)) ':' num2str(date(i,6)) ' UTC' ])
h = colorbar;
ylabel(h, '[m]')
colormap('lansey')
fontsize(16,16,16,16)
caxis([0 16])
shg
currFrame = getframe(gcf);
writeVideo(vidObj,currFrame);
end
    close(vidObj);
