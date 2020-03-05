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

%% The Code
% Purpose
% This code design to take data output from the model WAVEWATCHIII and
% compare it with nearby, contemporaneous satelitte altimeter estimations
% of significant wave height and wind speed

%% Part I: house cleaning
tic
close all
clear
[glyph baseDir] = giveGlyph;
%% Global Parameters
% There are two main methods for reducing altimeter data, which generally
% outputs data at 1 Hz. One way to do it is to use "find" and set a maximum
% distance and time. I vizualize this as a bubble and hence it is dubbed
% the bubble method. Another way to do it is to discritize the the grid and
% mean everything in each bin or cobe. This I call the box method. The box
% method, is currently about 10 times faster.

% averagingMethod = 'bubble';
averagingMethod = 'bubble';
% averagingMethod = 'none';

%% parameters for data averaging

maxTimeDiffMinutes =30; %minutes
maxTimeDiff = maxTimeDiffMinutes/(24*60); %days
maxDistance = 50; %km radius
minNumberObs = 5;

%% PART ?: load coastline data
% It will be useful to know the proximity of the comparison point to the
% coast where the general quality of altimeter data degrades

%% BUOY test


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

buoy = 'AO_MWB';
load([baseDir 'Analysis' glyph 'AltimerComparison' glyph 'mini-buoys' glyph buoy '.mat'])

time = [];
lon  = []; %from negative W to E poisitive 0 - 360 
lat  = [];
Hs   = [];

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

loadSatList = defineSatList(buoyTest.time);

% use model lat - lon to narrow down to files loaded, and load data from 
% local netCDF files 

obs = getAltimeterObs(loadSatList,buoyTest.lat, buoyTest.lon);

%% reduce data based on time and grid status

% concatenate data
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
        [LONobs LATobs TIMEobs HSobs HSobsSTD WINDobs WINDobsSTD indNaN] = bubbleMethodBuoy(LONobsNA, LATobsNA, TIMEobsNA, HSobsNA, WINDobsNA, buoyTest.lon, buoyTest.lat, buoyTest.time, maxDistance, maxTimeDiff, minNumberObs);
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
toc
%}

%% %% PART ?: collocation 
Ocstatsp(buoyTest.hs(indNaN),HSobs,.12)

% F = scatteredInterpolant(LONobs,LATobs,TIMEobs,HSobs,'linear','none');
% HSalt = F(papaLon,papaLat,papaTime);
% indNaN = ~isnan(HSalt);
% HSalt = HSalt(indNaN);

% phase 1: interpolation

%interp3d
%{
        LONobs  = vertcat(obs(:).lon);
        LATobs  = vertcat(obs(:).lat);
        TIMEobs = vertcat(obs(:).time);
        HSobs   = vertcat(obs(:).hsKcal);
        HSmd = interp1(papaTime, mdTest.hs, TIMEobs);
        
        indNaN = ~isnan(HSmd);
        HSmd = HSmd(indNaN);
        
        HSobs  = HSobs(indNaN);
        LONobs = LONobs(indNaN);
        LATobs = LATobs(indNaN);
        TIMEobs = TIMEobs(indNaN);
%}

% phase 2: kd tree
figure
plot(buoyTest.time, buoyTest.hs,'.')
hold on
errorbar(TIMEobs, HSobs, HSobsSTD)
grid on
datetick('x','yy')




%% save the matched pairs, clear data and loop through the whole domain

% save('yearResults2010')
