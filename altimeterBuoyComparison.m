function [pData rawAlt buoyData] = altimeterBuoyComparison(codePath,buoyPath,altPath,savePath,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Script to compare wave buoy measurements to contemperaneous altimeter   %
% measurements.                                                           %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%
%  OUTPUT
%  pData - structure of averaged wave height data from altimeter and paired
%  with buoy dara on a shared lon, lat, time
%  pData.lon - longitude in degrees 
%  pData.lat - latitude in degrees 
%  pData.time - time stamp
%  pData.altHs average altimeter significant wave height 
%  pData.wind average altimeter wind speed
%  pData.altHsStd standard deviation from  averaging method 
%  pData.altHsNoSam - number of samples
%  pData.buoyHs - buoy significant wave height
%  pData.options - documentation of options used
%  rawAlt - structure of 1 Hz altimeter data
%  rawAlt.lon 
%  rawAlt.lat 
%  rawAlt.time
%  rawAlt.Hs 
%  rawAlt.wind 
%
%  INPUT
%  Paths - in order path to the directory for code, directory with model
%  files, directory with altimeter altimeter data, and full path and file
%  name for saving results
%
%  OPTIONS structure (not optional)
%  options.averagingMethod = 1 = 'box'; 2 = 'bubble'; 3 = 'none'; Box
%  and bubble method are more fully described below. Default for buoy is
%  bubble method. Results are similar.
%
%  options for pairing data
%  options.maxTimeDiff = maxTimeDiffMinutes/(24*60); %maximum time
%                difference [days]
%  options.maxDistance = 25; %maximum space distance [km] (radius for
%                bubble method.
%  options.minNumberObs = 7; %minimum number of observations for
%                averaging
%  generic options
%  options.QC = 1; strictest quality control or 2 to include coastal data
%  options.save -> save output? logical 1 | 0
%
% e.g.
% codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
% buoyPath = '/Users/tripp/D/Analysis/altimeterComparison/mini-buoys/AO_MWB.mat';
% altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/Ribal_Young_2019/'; % path to Ribal and Young data
% savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/atlantic2016May-2017Nov';
% options.maxTimeDiff     =  30/(24*60);
% options.maxDistance     =  25;
% options.minNumberObs    =  7;
% options.QC              =  1;
% options.save            =  0;
% [pData rawAlt buoyData] = altimeterBuoyComparison(codePath,buoyPath,altPath,savePath,options)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Author:       Clarence Olin Collins III, Ph.D.                          %
% Affiliation:  ERDC - FRF                                                %
% created:      10/28/2019                                                %
% version:      1.0                                                       %
% updates:                                                                %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% 
% Purpose:
% This code design to take data output from a buoy or buoys  and
% compare it with nearby, contemporaneous satelitte altimeter estimations
% of significant wave height and wind speed

tic
%% Paths
path(codePath, path); %add code path to search path

%% Global Parameters
% There are two main methods for reducing altimeter data, which generally
% outputs data at 1 Hz. One way to do it is to use "find" and set a maximum
% distance and time. I vizualize this as a bubble and hence it is dubbed
% the bubble method. Another way to do it is to discritize the the grid and
% mean everything in each bin or cobe. This I call the box method. The box
% method, is currently about 10 times faster on gridded data, but doesn't
% make much sense for point measurements. Thus, only the bubble method is
% used here.

% averagingMethod = 'bubble';
% averagingMethod = 'none';

%% parameters for bubble method od altimeter data averaging

maxTimeDiff = options.maxTimeDiff;
maxDistance = options.maxDistance;
minNumberObs = options.minNumberObs;

%% PAPA TEST
%{
papaFile = 'D:\Datasets\Papa\166p1_historic.nc';
papaInfo = ncinfo(papaFile);
papaTimeOffset = datenum([1970 01 01 00 00 00]);
papaTime = double(ncread(papaFile,'waveTime'))/(60*60*24) + papaTimeOffset;
papaHs = ncread(papaFile,'waveHs');
papaGPStime = double(ncread(papaFile,'gpsTime'))/(60*60*24) + papaTimeOffset;
papaGPSlat = ncread(papaFile,'gpsLatitude');
papaGPSlon = ncread(papaFile,'gpsLongitude');
papaLat = interp1(papaGPStime, papaGPSlat, papaTime);
papaLon = interp1(papaGPStime, papaGPSlon, papaTime) + 360;
papaLat(1)= papaLat(2); %interpolation made this a NaN
papaLon(1)= papaLon(2); %interpolation made this a NaN
%}

% buoyData.time = papaTime;
% buoyData.lon  = papaLon;
% buoyData.lat  = papaLat;
% buoyData.hs   = papaHs;

%% CDIP 
% cdipURL = buoyPath;
% cdipTimeOffset = datenum([1970 01 01 00 00 00]);
% cdipTime = double(ncread(cdipURL,'waveTime'))/(60*60*24) + cdipTimeOffset;
% cdipHs = ncread(cdipURL,'waveHs');
% cdipGPStime = double(ncread(cdipURL,'gpsTime'))/(60*60*24) + cdipTimeOffset;
% cdipGPSlat = ncread(cdipURL,'gpsLatitude');
% cdipGPSlon = ncread(cdipURL,'gpsLongitude');
% cdipLat = interp1(cdipGPStime, cdipGPSlat, cdipTime);
% cdipLon = interp1(cdipGPStime, cdipGPSlon, cdipTime) + 360;
% cdipLat(1)= cdipLat(2); %interpolation made this a NaN
% cdipLon(1)= cdipLon(2); %interpolation made this a NaN
%
% 
% buoyData.time = cdipTime;
% buoyData.lon  = cdipLon;
% buoyData.lat  = cdipLat;
% buoyData.hs   = cdipHs;


%% NDBC 
% buoy = '46029';
% buoy = '41048';
% buoy = '46001';
% buoy = '44011';
% type = 'h';
% load([baseDir 'Analysis' glyph 'Jensen Buoy Work' glyph 'test' glyph buoy glyph buoy type '.mat'])
% buoyData.time = timeCat;
% buoyData.lon  = lonCat;
% buoyData.lat  = latCat;
% buoyData.hs   = hsCat;
% clear timeCat lonCat latCat hsCat

%% SIO MINI BUOYS

% load(buoyPath)
% 
% time = [];
% lon  = []; %from negative W to E poisitive 0 - 360 
% lat  = [];
% Hs   = [];
% 
% %concatenate data into one vector
% fields = fieldnames(MWB);
% for i = 1:length(fields)
%     time = [time; MWB.(fields{i}).time];
%     lon = [lon; MWB.(fields{i}).lon];
%     lat = [lat; MWB.(fields{i}).lat];
%     Hs = [Hs; MWB.(fields{i}).Hs];
% end
% 
% buoyData.time = time;
% buoyData.lon  = 360 + lon; %from negative W to E poisitive 0 - 360 
% buoyData.lat  = lat;
% buoyData.hs   = Hs;

%% ITOP
%EASI-N
% load([buoyPath 'EASI_Nv2.mat'])
% buoyData.time = EN.var.yday + datenum([2010 01 01 00 00 00]);
% buoyData.lon  =  126.968.*ones(length(buoyData.time),1);
% buoyData.lat  =  21.238.*ones(length(buoyData.time),1);
% buoyData.hs   = EN.par.int.Hm0;
% clear EN

%EASI-S
load([buoyPath 'EASI_Sv5.mat'])
buoyData.time = ES.var.yday + datenum([2010 01 01 00 00 00]);
buoyData.lon  =  127.258.*ones(length(buoyData.time),1);
buoyData.lat  =  19.683.*ones(length(buoyData.time),1);
buoyData.hs   = ES.par.int.Hm0;
clear ES

%% load altimeter data, version 1
% This code uses the Ribal and Young (2019) dataset. You can find it here:
% https://www.nature.com/articles/s41597-019-0083-9
% The dataset is large, ~100Gb.

loadSatList = defineSatList(buoyData.time, altPath);

% use model lat - lon to narrow down to files loaded, and load data from 
% local netCDF files 

obs = getAltimeterObsBuoy(loadSatList,buoyData.lat, buoyData.lon, altPath);

%% quality control
% flags: In the present database, a series of data flags defined as 1, 2,
% 3, 4, and 9 represent Good_data, Probably_good_data, SAR-mode data or
% possible hardware error (only used for CRYOSAT-2), Bad_data and
% Missing_data, respectively, have been used. We will retain only good data
% for now

for i = 1:length(obs)
    switch options.QC
        case 2
            qcPassInd = find(obs(i).hsKqc == 1 | obs(i).hsKqc == 2);
        otherwise
            qcPassInd = find(obs(i).hsKqc == 1);
    end
    obs(i).time = obs(i).time(qcPassInd);
    obs(i).lat = obs(i).lat(qcPassInd);
    obs(i).lon = obs(i).lon(qcPassInd);
    obs(i).hsKcal = obs(i).hsKcal(qcPassInd);
    obs(i).hsKqc = obs(i).hsKqc(qcPassInd);
    % obs(i).hsK = obs(i).hsK(qcPassInd );
    % obs(i).hsKno = obs(i).hsKno(qcPassInd );
    % obs(i).hsKstd = obs(i).hsKstd(qcPassInd);
    %WIND
    % obs(i).wind = obs(i).wind(qcPassInd );
    obs(i).windCal = obs(i).windCal(qcPassInd );
    % obsLength(i) = length(obs(1).time);
end
%% reduce data based on time and grid status

% concatenate altimeter data
LONobsNA  = vertcat(obs(:).lon);
LATobsNA  = vertcat(obs(:).lat);
TIMEobsNA = vertcat(obs(:).time);
HSobsNA   = vertcat(obs(:).hsKcal);
WINDobsNA = vertcat(obs(:).windCal);
DISTobsNA = vertcat(obs(:).dist2coast);

[obsIndx dum2] = find(TIMEobsNA >= min(buoyData.time) - maxTimeDiff...
    & TIMEobsNA <= max(buoyData.time) + maxTimeDiff);

rawAlt.lon = LONobsNA(obsIndx);
rawAlt.lat = LATobsNA(obsIndx);
rawAlt.time = TIMEobsNA(obsIndx);
rawAlt.Hs = HSobsNA(obsIndx);
rawAlt.wind = WINDobsNA(obsIndx);
rawAlt.dist = DISTobsNA(obsIndx);

%% Reducing data by averaging over space and time
% METHOD I
% this method takes the x-y-time grid from the model, then finds all the
% data within a maxDistance radius and within maxTimeDiff. This should form
% a bubble around the data and we get the average in that bubble.

% There may be a bug somewhere in the code because it never compares as well
% as the gridded method, which takes data from much further away.

% to do: add gaussian weighted average
[pData.lon pData.lat pData.time pData.altHsMean pData.altHsStd pData.altHsNearest...
    pData.altWindMean pData.altWindStd pData.altWindNearest pData.buoyIndNaN]...
    = bubbleMethodBuoy(LONobsNA,LATobsNA, TIMEobsNA, HSobsNA, WINDobsNA,...
    buoyData.lon, buoyData.lat,buoyData.time, maxDistance, maxTimeDiff, minNumberObs);
% average altimetered data paired with buoy data
pData.buoyHs = buoyData.hs(pData.buoyIndNaN);
pData.options = options;
pData.options.paths.buoyPath =buoyPath;
pData.options.paths.savePath = savePath;
pData.options.paths.codePath = codePath;
pData.options.paths.altPath = altPath;
% toc
%% Done
% The unaveraged altimeter data lives in *obsNA, and the averaged stuff in
% alt*, the original buoy data is in the structure buoyData. The paried
% stuff is paired*, alt*, buoy*.

% E.g., lon, lat, time, wave height are in pairedLon, pairedLat,
% pairedTime, buoyHs, altHs

%% save the matched pairs
if options.save
    save(savePath,'pData','rawAlt','buoyData')
end
