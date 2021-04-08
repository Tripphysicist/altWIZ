
% house cleaning
clear
close all

%% MODEL: Test for wind functionality
% this is tutorial for how to pair altimeter measurements with model data,
% calculate statistics and plot the results

% Part I: house cleaning
clear
close all

% Part II: deisgnate paths
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
mdPath   = '/Users/tripp/D/Analysis/altimeterComparison/modelData/test/';
altPath  = '/Volumes/LaCie/RY19/'; % path to Ribal and Young data
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/wind_test';
figurePath = '/Users/tripp/D/Analysis/altimeterComparison/figures/wind_test';

% Part III: set options
options.altDatabase     =  'RY19';
options.averagingMethod =  'none'; % 'none' or 'bubble'
options.maxTimeDiff     =  30/(24*60);
options.maxDistance     =  50;
options.minNumberObs    =  5;
options.loopSize        =  3;
options.QC              =  1;
options.save            =  1;
options.seasonalPlots   =  0;
options.runWindEval     =  1;

% Part IV: run the pairing code
[pData] = altimeterModelPairing(codePath,mdPath,altPath,savePath,options);

% Part V: make the plots
% if the pairing has alread been saved, then just load in the pdata and run
% the following. This section has addtional dependancies. JLab, m_map,
% export_fig

% load(savePath)
altimeterModelComparisonPlotting
%% BUOY: CDIP 430, 26m duck pier 
% deisgnate paths
tic
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
buoyInfo.source = 'CDIP';
buoyInfo.name = '166'; %CDIP ID for Ocean Station Papa
altPath = '/Volumes/LaCie/RY19/'; % path to ESA
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/CDIP_430.mat';

% set options
options.altDatabase     = 'RY19';
options.maxTimeDiff     =  30/(24*60);
options.maxDistance     =  100;
options.minNumberObs    =  1;
options.QC              =  1;
options.save            =  0;
[pData , rawAlt , buoyData] = altimeterBuoyPairing(codePath,buoyInfo,altPath,savePath,options)
toc
%% BUOY: NDBC 41002, Atlantic Ocean, RY19 QC-1 vs QC-2
% This was done for comparing QC-1 and QC-2 in response to comments on MDPI 
% paper on altimeter data and hurricane

% house cleaning
clear
close all

% deisgnate paths
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
buoyInfo.source = 'NDBC';
buoyInfo.name = '46075';% station number
buoyInfo.lat = 53.983;    % N
buoyInfo.lon = -160.817; % E
altPath = '/Volumes/LaCie/RY19/'; % path to Ribal and Young data

% set options
options.altDatabase     = 'RY19';
options.maxTimeDiff     =  60/(24*60);
options.maxDistance     =  50;
options.minNumberObs    =  5;
options.QC              =  1;
options.save            =  1;

savePath = ['/Users/tripp/D/Analysis/altimeterComparison/data/buoys/'...
    buoyInfo.source buoyInfo.name options.altDatabase 'qc' num2str(options.QC)];

[pData rawAlt buoyData] = altimeterBuoyPairing(codePath,buoyInfo,altPath,savePath,options)

%% MODEL: WIS Atlantic Field Files for hurricane months

% Part I: house cleaning
clear
close all
% Part II: deisgnate paths
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
path(codePath, path); %add code path to search path
mdPath   = '/Users/tripp/D/Analysis/altimeterComparison/modelData/atlanticTropicalCyclones/WIS/';
altPath  = '/Volumes/LaCie/RY19/'; % path to Ribal and Young data
savePath = '/Users/tripp/D/Analysis/altimeterTropicalCyclones/figures/WIS/allData/';
figurePath = '/Users/tripp/D/Analysis/altimeterTropicalCyclones/figures/WIS/allData/';

% Part III: set options
options.altDatabase     =  'RY19';
options.averagingMethod =  'none'; % 'none' or 'bubble'
options.maxTimeDiff     =  30/(24*60);
options.maxDistance     =  50;
options.minNumberObs    =  5;
options.loopSize        =  3;
options.QC              =  2;
options.save            =  1;
options.seasonalPlots   =  1;

% Part IV: run the pairing code
[pData] = altimeterModelPairing(codePath,mdPath,altPath,savePath,options);

% Part V: make the plots
% if the pairing has alread been saved, then just load in the pdata and run
% the following. This section has addtional dependancies. JLab, m_map,
% export_fig

% load(savePath)
altimeterModelComparisonPlotting


%% STORM: get all North Atlantic Hurricane Data RY19 - QC-2 
clear
close all
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
stormName = 'all';
stormYear = 1985:2020;
altPath = '/Volumes/LaCie/RY19/'; % path to Ribal and Young data
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/atlanticHurricanes/altimeterHurricaneTCOBS1985v2.1test';
options.altDatabase     = 'RY19';
options.maxTimeDiff     =  30/(24*60); % 1/2 hour
options.maxDistance     =  500; % km
options.QC              =  2; % 
options.save            =  1;
[stormObs, stormI] = altimeterStormPairing(stormName,stormYear,codePath,altPath,savePath,options);


%% individual storms for Peter: Irene (2011), Sandy (2012), Matthew (2016), 
% Irma (2017), Florence (2018), Dorian (2019)
clear
close all
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
stormName = 'FLORENCE';
stormYear = 2018;
% altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/RY19/'; % path to Ribal and Young data
altPath  = '/Volumes/LaCie/ESA_sea_state/'; % path to ESA
savePath = '/Users/tripp/D/Analysis/altimeterTropicalCyclones/forPeter/Florence2018-1000km';
% options.altDatabase     = 'RY19';
options.altDatabase     = 'ESA';
options.maxTimeDiff     =  1*60/(24*60); % 1 hour
options.maxDistance     =  1000; % km
options.QC              =  1; % only offshore data
options.save            =  1;
[stormObs, stormI] = altimeterStormPairing(stormName,stormYear,codePath,altPath,savePath,options);

%% STORM: get all North Atlantic Hurricane Data IBTrACS offshore data only (ESA)
clear
close all
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
stormName = 'all';
stormYear = 1985:2020;
altPath = '/Volumes/LaCie/ESA_sea_state/'; % path to Ribal and Young data
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/altimeterHurricaneTCOBS1985ESA';
options.altDatabase     = 'ESA';
options.maxTimeDiff     =  1*60/(24*60); % 1 hour
options.maxDistance     =  500; % km
options.QC              =  1; % only offshore data
options.save            =  1;
[stormObs, stormI] = altimeterStormPairing(stormName,stormYear,codePath,altPath,savePath,options);
%% STORM: get all North Atlantic Hurricane Data IBTrACS offshore data only
clear
close all
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
stormName = 'all';
stormYear = 1985:2020;
altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/RY19/'; % path to Ribal and Young data
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/altlanticHurricanes/altimeterHurricaneTCOBS1985-1000km';
options.altDatabase     = 'RY19';
options.maxTimeDiff     =  1*60/(24*60); % 1 hour
options.maxDistance     =  1000; % km
options.QC              =  2; % only offshore data
options.save            =  1;
[stormObs, stormI] = altimeterStormPairing(stormName,stormYear,codePath,altPath,savePath,options);


%% MODEL: WIS Atlantic 2018 box method
% this is tutorial for how to pair altimeter measurements with model data,
% calculate statistics and plot the results

% Part I: house cleaning
clear
close all
% Part II: deisgnate paths
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
path(codePath, path); %add code path to search path
mdPath   = '/Users/tripp/D/Analysis/altimeterComparison/modelData/WIS-PAC-2018_fields/basin_l1/';
altPath  = '/Users/tripp/D/Datasets/Satellite/Altimeter/RY19/'; % path to Ribal and Young data
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/WIS/WIS-PAC-2018-basin-l1';
figurePath = '/Users/tripp/D/Analysis/altimeterComparison/figures/WIS-PAC-2018-basin-l1';

% Part III: set options
options.altDatabase     =  'RY19';
options.averagingMethod =  'none'; % 'none' or 'bubble'
options.maxTimeDiff     =  30/(24*60);
options.maxDistance     =  50;
options.minNumberObs    =  5;
options.loopSize        =  3;
options.QC              =  1;
options.save            =  1;
options.seasonalPlots   =  1;

% Part IV: run the pairing code
[pData] = altimeterModelPairing(codePath,mdPath,altPath,savePath,options);

% Part V: make the plots
% if the pairing has alread been saved, then just load in the pdata and run
% the following. This section has addtional dependancies. JLab, m_map,
% export_fig

% load(savePath)
altimeterModelComparisonPlotting

%% BUOY: NDBC 46029 local copy
% The buoy is currently not fully functional. One has to manually put in
% buoy information unless its a CDIP location. Pairing uses the bubble
% method keeping both the average and the nearest (in space) altimeter
% measurement.

% house cleaning
clear
close all
% deisgnate paths
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel/';
buoyInfo.source = '3DMG'; %SIO mini buoys
buoyInfo.name = ''; %SIO mini buoys
altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/RY19/'; % path to ESA
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/NDBC_46029_3DMG_RY19.mat';

% set options
options.altDatabase     = 'RY19';
options.maxTimeDiff     =  30/(24*60); %days
options.maxDistance     =  50; %km
options.minNumberObs    =  1;
options.QC              =  2;
options.save            =  1;

[pData , rawAlt , buoyData] = altimeterBuoyPairing(codePath,buoyInfo,altPath,savePath,options);
%
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel/';
buoyInfo.source = 'HIPPY'; %SIO mini buoys
buoyInfo.name = ''; %SIO mini buoys
altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/RY19/'; % path to ESA
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/NDBC_46029_HIPPY_RY19.mat';

% set options
options.altDatabase     = 'RY19';
options.maxTimeDiff     =  30/(24*60); %days
options.maxDistance     =  50; %km
options.minNumberObs    =  1;
options.QC              =  2;
options.save            =  1;

[pData , rawAlt , buoyData] = altimeterBuoyPairing(codePath,buoyInfo,altPath,savePath,options);

codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel/';
buoyInfo.source = 'WR'; %SIO mini buoys
buoyInfo.name = ''; %SIO mini buoys
altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/RY19/'; % path to ESA
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/NDBC_46029_WR_RY19.mat';

% set options
options.altDatabase     = 'RY19';
options.maxTimeDiff     =  30/(24*60); %days
options.maxDistance     =  50; %km
options.minNumberObs    =  1;
options.QC              =  2;
options.save            =  1;

[pData , rawAlt , buoyData] = altimeterBuoyPairing(codePath,buoyInfo,altPath,savePath,options);

%% BUOY: SIO mini buoys, RY19 & ESA
% The buoy is currently not fully functional. One has to manually put in
% buoy information unless its a CDIP location. Pairing uses the bubble
% method keeping both the average and the nearest (in space) altimeter
% measurement.

% house cleaning
clear
close all
% deisgnate paths
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel/';
buoyInfo.source = 'SIO'; %SIO mini buoys
buoyInfo.name = ''; %SIO mini buoys
altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/RY19/'; % path to ESA
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/SIOminiBuoyRY19.mat';

% set options
options.altDatabase     = 'RY19';
options.maxTimeDiff     =  60/(24*60);
options.maxDistance     =  100;
options.minNumberObs    =  1;
options.QC              =  2;
options.save            =  1;

[pData , rawAlt , buoyData] = altimeterBuoyPairing(codePath,buoyInfo,altPath,savePath,options)

% ESA
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/SIOminiBuoyESA.mat';
options.altDatabase     = 'ESA';
options.QC              =  1;
altPath  = '/Volumes/LaCie/ESA_sea_state/'; % path to ESA
[pData , rawAlt , buoyData] = altimeterBuoyPairing(codePath,buoyInfo,altPath,savePath,options)


%% BUOY: CDIP 166, Ocean Station Papa, RY19
% The buoy is currently not fully functional. One has to manually put in
% buoy information unless its a CDIP location. Pairing uses the bubble
% method keeping both the average and the nearest (in space) altimeter
% measurement.

% house cleaning
clear
close all
% deisgnate paths
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
buoyInfo.source = 'CDIP';
buoyInfo.name = '179'; %CDIP ID for Ocean Station Papa
altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/RY19/'; % path to ESA
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/CDIP_179.mat';

% set options
options.altDatabase     = 'RY19';
options.maxTimeDiff     =  30/(24*60);
options.maxDistance     =  50;
options.minNumberObs    =  5;
options.QC              =  2;
options.save            =  1;
[pData , rawAlt , buoyData] = altimeterBuoyPairing(codePath,buoyInfo,altPath,savePath,options)

%% MODEL: WIS Atlantic 2017 Feburary box
% this is tutorial for how to pair altimeter measurements with model data,
% calculate statistics and plot the results

% Part I: house cleaning
clear
close all

% Part II: deisgnate paths
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
mdPath   = '/Users/tripp/D/Analysis/altimeterComparison/modelData/test/';
altPath  = '/Volumes/LaCie/ESA_sea_state/'; % path to Ribal and Young data
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/Atlantic_2017_Feb_ESA_box';
figurePath = '/Users/tripp/D/Analysis/altimeterComparison/figures/Atlantic_2017_Feb_ESA_box';

% Part III: set options
options.altDatabase     =  'ESA';
options.averagingMethod =  'box'; % 'none' or 'bubble'
options.maxTimeDiff     =  30/(24*60);
options.maxDistance     =  50;
options.minNumberObs    =  5;
options.loopSize        =  3;
options.QC              =  1;
options.save            =  1;
options.seasonalPlots   =  0;

% Part IV: run the pairing code
[pData] = altimeterModelPairing(codePath,mdPath,altPath,savePath,options);

% Part V: make the plots
% if the pairing has alread been saved, then just load in the pdata and run
% the following. This section has addtional dependancies. JLab, m_map,
% export_fig

% load(savePath)
altimeterModelComparisonPlotting

%% MODEL: WIS Atlantic 2017 Feburary bubble
% this is tutorial for how to pair altimeter measurements with model data,
% calculate statistics and plot the results

% Part I: house cleaning
clear
close all

% Part II: deisgnate paths
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
mdPath   = '/Users/tripp/D/Analysis/altimeterComparison/modelData/test/';
altPath  = '/Volumes/LaCie/ESA_sea_state/'; % path to Ribal and Young data
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/Atlantic_2017_Feb_ESA_bubble';
figurePath = '/Users/tripp/D/Analysis/altimeterComparison/figures/Atlantic_2017_Feb_ESA_bubble';

% Part III: set options
options.altDatabase     =  'ESA';
options.averagingMethod =  'bubble'; % 'none' or 'bubble'
options.maxTimeDiff     =  30/(24*60);
options.maxDistance     =  50;
options.minNumberObs    =  5;
options.loopSize        =  8;
options.QC              =  1;
options.save            =  1;
options.seasonalPlots   =  0;

% Part IV: run the pairing code
[pData] = altimeterModelPairing(codePath,mdPath,altPath,savePath,options);

% Part V: make the plots
% if the pairing has alread been saved, then just load in the pdata and run
% the following. This section has addtional dependancies. JLab, m_map,
% export_fig

% load(savePath)
altimeterModelComparisonPlotting
%% BUOY: CDIP 166, Ocean Station Papa, RY19
% The buoy is currently not fully functional. One has to manually put in
% buoy information unless its a CDIP location. Pairing uses the bubble
% method keeping both the average and the nearest (in space) altimeter
% measurement.

% house cleaning
clear
close all

% deisgnate paths
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
buoyInfo = '166'; %CDIP ID for Ocean Station Papa
altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/RY19/'; % path to Ribal and Young data
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/papaTestv2.0';

% set options
options.altDatabase     = 'RY19';
options.maxTimeDiff     =  30/(24*60);
options.maxDistance     =  50;
options.minNumberObs    =  5;
options.QC              =  1;
options.save            =  1;
[pData rawAlt buoyData] = altimeterBuoyPairing(codePath,buoyInfo,altPath,savePath,options)

%% MODEL: get WIS Atlantic 2017 basin level 1 (basin_l1)
% this is tutorial for how to pair altimeter measurements with model data,
% calculate statistics and plot the results

% house cleaning
clear
close all

% deisgnate paths
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
mdPath   = '/Users/tripp/D/Analysis/altimeterComparison/modelData/WIS-ATL-2017_fields/eastc_l2/';
altPath  = '/Volumes/LaCie/ESA_sea_state/'; % path to ESA sea state
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/WIS-ATL-2017_eastc_l2_ESA';
figurePath = '/Users/tripp/D/Analysis/altimeterComparison/figures/WIS-ATL-2017_eastc_l2_ESA';

% set options
options.altDatabase     =   'ESA';
options.averagingMethod =  3;
options.maxTimeDiff     =  30/(24*60);
options.maxDistance     =  50;
options.minNumberObs    =  5;
options.loopSize        =  3;
options.QC              =  1;
options.save            =  1;
options.seasonalPlots   =  1;

% run the pairing code
[pData] = altimeterModelPairing(codePath,mdPath,altPath,savePath,options);

% calculate stats and make the plots
% if the pairing has alread been saved, then just load in the pdata and run
% the following. This section has addtional dependancies. JLab, m_map,
% export_fig

% load(savePath)
altimeterModelComparisonPlotting
%% MODEL: get WIS Atlantic 2017 basin level 1 (basin_l1)
% this is tutorial for how to pair altimeter measurements with model data,
% calculate statistics and plot the results

% Part I: house cleaning
clear
close all

% Part II: deisgnate paths
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
mdPath   = '/Users/tripp/D/Analysis/altimeterComparison/modelData/WIS-ATL-2017_fields/basin_l1/';
altPath  = '/Users/tripp/D/Datasets/Satellite/Altimeter/RY19/'; % path to Ribal and Young data
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/WIS-ATL-2017_basin_l1';
figurePath = '/Users/tripp/D/Analysis/altimeterComparison/figures/WIS-ATL-2017_basin_l1';

% Part III: set options
options.altDatabase     =   'RY19';
options.averagingMethod =  1;
options.maxTimeDiff     =  30/(24*60);
options.maxDistance     =  50;
options.minNumberObs    =  5;
options.loopSize        =  3;
options.QC              =  1;
options.save            =  1;
options.seasonalPlots   =  1;

% Part IV: run the pairing code
[pData] = altimeterModelPairing(codePath,mdPath,altPath,savePath,options);

% Part V: make the plots
% if the pairing has alread been saved, then just load in the pdata and run
% the following. This section has addtional dependancies. JLab, m_map,
% export_fig

% load(savePath)
% altimeterModelComparisonPlotting

%% STORM: get all North Atlantic Hurricane Data ESA
clear
close all
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
stormName = 'all';
stormYear = 1985:2020;
altPath  = '/Volumes/LaCie/ESA_sea_state/'; % path to ESA sea state
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/altimeterHurricaneESA';
options.altDatabase     = 'ESA';
options.maxTimeDiff     =  3*60/(24*60); % 3 hour
options.maxDistance     =  250; % km
options.QC              =  1; % take coastal data
options.save            =  1;
[stormObs, stormI] = altimeterStormPairing(stormName,stormYear,codePath,altPath,savePath,options);

%% STORM: get all North Atlantic Hurricane Data
clear
close all
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
stormName = 'all';
stormYear = 1985:2020;
altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/RY19/'; % path to Ribal and Young data
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/altimeterHurricane1985';
options.altDatabase     = 'RY19';
options.maxTimeDiff     =  3*60/(24*60); % 3 hour
options.maxDistance     =  250; % km
options.QC              =  2; % take coastal data
options.save            =  1;
[stormObs, stormI] = altimeterStormPairing(stormName,stormYear,codePath,altPath,savePath,options);


%% MODEL: get WIS Atlantic 2017 basin level 2: east coast (eastc_l2)
% this is a model run of the altimeter pairing code, the final output is a
% paired data structure. This is also saved at the loacation indicated in
% the save path, if save option is used.

clear
close all
codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
mdPath   = '/Users/tripp/D/Analysis/altimeterComparison/modelData/WIS-ATL-2017_fields/eastc_l2/';
altPath  = '/Users/tripp/D/Datasets/Satellite/Altimeter/RY19/'; % path to Ribal and Young data
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/WIS-ATL-2017_eastc_l2';
options.altDatabase     =   'RY19';
options.averagingMethod =  3;
options.maxTimeDiff     =  30/(24*60);
options.maxDistance     =  50;
options.minNumberObs    =  5;
options.loopSize        =  3;
options.QC              =  2;
options.save            =  1;
[pData] = altimeterModelPairing(codePath,mdPath,altPath,savePath,options);

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
