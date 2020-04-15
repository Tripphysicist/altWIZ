% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Script to compare buoy output to altimeter measurements of significant  %
% wave height and wind speed                                              %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Author:       Clarence Olin Collins III, Ph.D.                          %
% Affiliation:  ERDC - FRF                                                %
% created:      02/25/2020                                                %
% version:      1.0                                                       %
% updates:                                                                %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

%
% Purpose
% This code was designed to take data of wave height and wind speed from 
% a buoy and compare it with nearby, contemporaneous satelitte altimeter 
% estimations of wave height and wind speed

%% house cleaning
tic % start a timer
close all
clear

%% Paths

codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
path(codePath, path); %add code path to search path
buoyPath = '/Users/tripp/D/Analysis/altimeterComparison/mini-buoys/';

% This code uses the Ribal and Young (2019) dataset. You can find it here:
% https://www.nature.com/articles/s41597-019-0083-9
% The dataset is large, ~100Gb. Please indicate the full path to the
% directory here:

altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/Ribal_Young_2019/'; % path to Ribal and Young data
% choose patha and file name to save your results
savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/atlantic2016May-2017Nov'; 

% The following is to ensure compatibility across platforms, only glyph
% matters
[glyph baseDir] = giveGlyph;

%% Global Parameters
% There are two main methods for reducing altimeter data, which generally
% outputs data at 1 Hz. One way to do it is to use "find" and set a maximum
% distance and time. I vizualize this as a bubble and hence it is dubbed
% the bubble method. Another way to do it is to discritize the the grid and
% mean everything in each bin or cobe. This I call the box method. The box
% method, is currently about 10 times faster on gridded data, but doesn't
% make much sense for point measurements. Thus, only the bubble method is
% used here.

averagingMethod = 'bubble';
% averagingMethod = 'none';

%% parameters for bubble method od altimeter data averaging

maxTimeDiffMinutes =30; %minutes
maxTimeDiff = maxTimeDiffMinutes/(24*60); %days
maxDistance = 100; %km radius
minNumberObs = 5;

%% BUOY data
% PAPA TEST
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

% buoyTest.time = papaTime;
% buoyTest.lon  = papaLon;
% buoyTest.lat  = papaLat;
% buoyTest.hs   = papaHs;

% NDBC 
% buoy = '46029';
% buoy = '41048';
% buoy = '46001';
% buoy = '44011';
% type = 'h';
% load([baseDir 'Analysis' glyph 'Jensen Buoy Work' glyph 'test' glyph buoy glyph buoy type '.mat'])
% buoyTest.time = timeCat;
% buoyTest.lon  = lonCat;
% buoyTest.lat  = latCat;
% buoyTest.hs   = hsCat;
% clear timeCat lonCat latCat hsCat

% CDIP MINI BUOYS
% enter the full path to your buoy data here
buoy = 'AO_MWB';
buoyDataFile = [buoyPath buoy '.mat'];
load(buoyDataFile)

time = [];
lon  = []; %from negative W to E poisitive 0 - 360 
lat  = [];
Hs   = [];

%concatenate data into one vector
fields = fieldnames(MWB);
for i = 1:length(fields)
    time = [time; MWB.(fields{i}).time];
    lon = [lon; MWB.(fields{i}).lon];
    lat = [lat; MWB.(fields{i}).lat];
    Hs = [Hs; MWB.(fields{i}).Hs];
end

buoyTest.time = time;
buoyTest.lon  = 360 + lon; %from negative W to E poisitive 0 - 360 
buoyTest.lat  = lat;
buoyTest.hs   = Hs;

%% PART ?: load altimeter data, version 1
% using the Ribal and Young database

loadSatList = defineSatList(buoyTest.time, altPath);

% use model lat - lon to narrow down to files loaded, and load data from 
% local netCDF files 

obs = getAltimeterObsBuoy(loadSatList,buoyTest.lat, buoyTest.lon, altPath);

%% quality control
% flags: In the present database, a series of data flags defined as 1, 2,
% 3, 4, and 9 represent Good_data, Probably_good_data, SAR-mode data or
% possible hardware error (only used for CRYOSAT-2), Bad_data and
% Missing_data, respectively, have been used. We will retain only good data
% for now

for i = 1:length(obs)
    qcPassInd = find(obs(i).hsKqc ==1);
    
    obs(i).time = obs(i).time(qcPassInd );
    obs(i).lat = obs(i).lat(qcPassInd );
    obs(i).lon = obs(i).lon(qcPassInd );
    obs(i).hsKcal = obs(i).hsKcal(qcPassInd );
    obs(i).hsKqc = obs(i).hsKqc(qcPassInd );
    %     obs(i).hsK = obs(i).hsK(qcPassInd );
    %     obs(i).hsKno = obs(i).hsKno(qcPassInd );
    %     obs(i).hsKstd = obs(i).hsKstd(qcPassInd );|
    %WIND
    %     obs(i).wind = obs(i).wind(qcPassInd );
    obs(i).windCal = obs(i).windCal(qcPassInd );
    obsLength(i) = length(obs(1).time);
end
%% reduce data based on time and grid status

% concatenate altimeter data
LONobsNA  = vertcat(obs(:).lon);
LATobsNA  = vertcat(obs(:).lat);
TIMEobsNA = vertcat(obs(:).time);
HSobsNA   = vertcat(obs(:).hsKcal);
WINDobsNA = vertcat(obs(:).windCal);

[obsIndx dum2] = find(TIMEobsNA >= min(buoyTest.time) - maxTimeDiff...
    & TIMEobsNA <= max(buoyTest.time) + maxTimeDiff);

LONobsNA  = LONobsNA(obsIndx);
LATobsNA  = LATobsNA(obsIndx);
TIMEobsNA = TIMEobsNA(obsIndx);
HSobsNA   = HSobsNA(obsIndx);
WINDobsNA = WINDobsNA(obsIndx);

%% Reducing data by averaging over space and time
% METHOD I
% this method takes the x-y-time grid from the model, then finds all the
% data within a maxDistance radius and within maxTimeDiff. This should form
% a bubble around the data and we get the average in that bubble.

% There may be a bug somewhere in the code because it never compares as well
% as the gridded method, which takes data from much further away.

% to do: add gaussian weighted average
tic
switch averagingMethod
    case 'bubble'
        [pairedLon pairedLat pairedTime altHs altHsStd altWind altWindStd buoyIndNaN] = bubbleMethodBuoy(LONobsNA, LATobsNA, TIMEobsNA, HSobsNA, WINDobsNA, buoyTest.lon, buoyTest.lat, buoyTest.time, maxDistance, maxTimeDiff, minNumberObs);
    %{
    case 'box'
        
        %% Method 2
        % This creates a gridded cube with demsions lat, lon, & time and
        % the average in each of those cubes. It takes too long (added 30 minutes
        % to a 2 x 2 degree grid (from like 30 seconds).
        
        [LONobs LATobs TIMEobs HSobs] = boxMethodBuoy(LONobsNA, LATobsNA, TIMEobsNA, HSobsNA, buoyTest.lon, buoyTest.lat, maxTimeDiff, minNumberObs);
        %}
    case 'none'
        
        LONobs  = LONobsNA;
        LATobs  = LATobsNA;
        TIMEobs = TIMEobsNA;
        HSobs   = HSobsNA;
        
        clear LONobsNA LATobsNA TIMEobsNA HSobsNA
end

% paired buoy data
buoyHs = buoyTest.hs(buoyIndNaN);




toc
%% Done
% The unaveraged altimeter data lives in *obsNA, and the averaged stuff in
% alt*, the original buoy data is in the structure buoyTest. The paried
% stuff is paired*, alt*, buoy*.

% E.g., lon, lat, time, wave height are in pairedLon, pairedLat,
% pairedTime, buoyHs, altHs

%% save the matched pairs
% save(savePath)
