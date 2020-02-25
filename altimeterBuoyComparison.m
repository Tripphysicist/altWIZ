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

% Credits
% This code is very much based on the 2019 Wave Summer School, WAVEWATCHIII
% tutorial (Day 5), written by Stylianos 'Stelios' Flampouris
% (stylianos.flampouris@gmail.com) and taught by Stelios and Ricardo M.
% Campos (riwave@gmail.com). The code structure is derrivative of their
% original work.

%% Part I: house cleaning
tic
close all
clear

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

%% PART ?: load coastline data
% It will be useful to know the proximity of the comparison point to the
% coast where the general quality of altimeter data degrades

%% PAPA test
% this was a test designed to interpolate altimeter data to a buoy
% position, the buoy chosen was ocean station papa. The interpolation works
% but is not meaningful because the altimeter data is interpolated over
% long distances. I think we need a way to reduce altimeter data, or to
% inforce some gap length. To run the test again, simply uncomment the PAPA
% test sections.

% PAPA TEST

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


%% PART ?: load WW3 data, version 1
% phase 1: one .nc file ACTIVE
% phase 2: loop through a directory

%WIND
% UmdCat = [];
% VmdCat = [];

%WIND
% calculate wind from u and v components
% uwnd = ncread([mdPath 'ww3.201010.nc'], 'uwnd');
% vwnd = ncread([mdPath 'ww3.201010.nc'], 'vwnd');
% mdCol.wind = sqrt(uwnd.^2 + vwnd.^2);

% PAPA TEST

%PAPA TEST
buoyTest.time = papaTime;
buoyTest.lon  = papaLon;
buoyTest.lat  = papaLat;
buoyTest.hs   = papaHs;

%% PART ?: load altimeter data, version 1
% using the Ribal and Young database

% Altimeter	  Freq.-Band Latitude-coverage Initial-Date Final-Date
% GEOSAT	  Ku	     -73 to 72	       31/03/1985	31/12/1989
% ERS-1	      Ku	     -81.5 to 81.5	   01/08/1991	02/06/1996
% TOPEX	      Ku C	     -66 to 66	       25/09/1992	08/10/2005
% ERS-2	      Ku	     -81.5 to 81.5	   29/04/1995	11/05/2009
% GFO	      Ku	     -73 to 72	       07/06/2000	07/09/2008
% JASON-1	  Ku C	     -66.15 to 66.15   15/01/2002	21/06/2013
% ENVISAT	  Ku S	     -82 to 82	       14/05/2002	08/04/2012
% JASON-2	  Ku C	     -66.15 to 66.15   04/07/2008	Ongoing
% CRYOSAT-2	  Ku	     -88 to 88	       14/07/2010	Ongoing
% HY-2A	      Ku C	     -81 to 80	       01/10/2011	06/06/2018
% SARAL	      Ka	     -81.49 to 81.49   14/03/2013	Ongoing
% JASON-3	  Ku C	     -66.15 to 66.15   12/02/2016	Ongoing
% SENTINEL-3A Ku C	     -78 to 81	       01/03/2016	Ongoing

%what to do with altimeters with 2 bands?

altPath = 'd:\Datasets\Satellite\Altimeter\Ribal_Young_2019\';
satList  = dir(altPath);
satList = rmfield (satList,{'date','bytes','isdir','datenum'});
satList = satList(3:end);

for i = 1:length(satList)
    satFileList = dir([altPath satList(i).name]);
    satFileList = rmfield (satFileList,{'date','bytes','isdir','datenum'});
    satFileList = satFileList(3:end);
    info = ncinfo([altPath satList(i).name '\' satFileList(1).name]);
    %    timeOffsetMeta(i) = info.Variables(1).Attributes(3);
    altTime = ncread([altPath satList(i).name '\' satFileList(1).name],'TIME');
    satList(i).timeStart = altTime(1)+datenum([1985 01 01 00 00 00]);
    satList(i).timeEnd = altTime(end)+datenum([1985 01 01 00 00 00]);
end

% use model time to neglect some of the altimeter data and build a list
count1 = 0;
for i = 1:length(satList)
    satSpan = floor(satList(i).timeStart):ceil(satList(i).timeEnd);
    mdSpan  = floor(buoyTest.time(1)):ceil(buoyTest.time(end));
    
    if  intersect(satSpan,mdSpan)
        count1 = count1 +1;
        loadSatList(count1) = satList(i);
    end
end

%% use model lat - lon to narrow down to files loaded
count2 = 0;
for i = 1:length(loadSatList)
    altFilePath = [altPath loadSatList(i).name '\'];
    altFileList = dir(altFilePath);
    altFileList = rmfield (altFileList,{'date','bytes','isdir','datenum'});
    altFileList = altFileList(3:end);
    for j = 1:length(altFileList)
        altNorth  = str2num(altFileList(j).name(end-16:end-14));
        if altFileList(j).name(end-13) == 'S'
            altNorth = -altNorth;
        else
        end
        altEast   = str2num(altFileList(j).name(end-11:end-9));
        
        % search for data that is within 1 degree of the min and max buoy coordinate
        if altNorth <= max(buoyTest.lat) + 1 & altNorth >= min(buoyTest.lat) - 1 & ...
                altEast <= max(buoyTest.lon) + 1 & altEast >= min(buoyTest.lon) -1
            %         if altNorth >= 40 & altNorth <= 41 & ...
            %                 altEast >= 200 & altEast <= 201
            fileName = [altFilePath altFileList(j).name];
            %            altInfo = ncinfo(fileName);
            currentAltInfo = ncinfo(fileName);
            count2 = count2 + 1;
            
            % load data
            satTimeOffset = datenum([1985 01 01 00 00 00]); %same for all
            
            currentTime   = ncread(fileName, 'TIME') + satTimeOffset;
            
            obs(count2).time = ncread(fileName, 'TIME') + satTimeOffset;
            obs(count2).lat  = ncread(fileName, 'LATITUDE');
            obs(count2).lat  = double(obs(count2).lat);
            obs(count2).lon  = ncread(fileName, 'LONGITUDE');
            obs(count2).lon  = double(obs(count2).lon);
            
            % try C band
            %             obs(count2).hsC   = ncread(fileName, 'SWH_C');
            %             obs(count2).hsCqc = ncread(fileName, 'SWH_C_quality_control');
            %             obs(count2).hsCno = ncread(fileName, 'SWH_C_num_obs');
            %             obs(count2).hsCstd = ncread(fileName, 'SWH_C_std_dev');
            
            if  strcmp(loadSatList(i).name,'SARAL')
                %                 obs(count2).hsK   = ncread(fileName, 'SWH_KA');
                obs(count2).hsKcal   = ncread(fileName, 'SWH_KA_CAL');
                obs(count2).hsKqc = ncread(fileName, 'SWH_KA_quality_control');
                %                 obs(count2).hsKno = ncread(fileName, 'SWH_KA_num_obs');
                %                 obs(count2).hsKstd = ncread(fileName, 'SWH_KA_std_dev');
            else
                %                 obs(count2).hsK   = ncread(fileName, 'SWH_KU');
                obs(count2).hsKcal   = ncread(fileName, 'SWH_KU_CAL');
                obs(count2).hsKqc = ncread(fileName, 'SWH_KU_quality_control');
                %                 obs(count2).hsKno = ncread(fileName, 'SWH_KU_num_obs');
                %                 obs(count2).hsKstd = ncread(fileName, 'SWH_KU_std_dev');
            end
            
            %WIND
            %             obs(count2).wind    = ncread(fileName, 'WSPD');
            %             obs(count2).windCal = ncread(fileName, 'WSPD_CAL');
            
            
        end
    end
end


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
    %     obs(i).windCal = obs(i).windCal(qcPassInd );
    obsLength(i) = length(obs(1).time);
end


%% reduce data based on time and grid status

% concatenate data
LONobsNA  = vertcat(obs(:).lon);
LATobsNA  = vertcat(obs(:).lat);
TIMEobsNA = vertcat(obs(:).time);
HSobsNA   = vertcat(obs(:).hsKcal);
%WIND
% WINDobsNA = vertcat(obs(:).windCal);

halfHour = 0.0208;

[obsIndx dum2] = find(TIMEobsNA >= min(buoyTest.time) - halfHour...
    & TIMEobsNA <= max(buoyTest.time) + halfHour);

LONobsNA  = LONobsNA(obsIndx);
LATobsNA  = LATobsNA(obsIndx);
TIMEobsNA = TIMEobsNA(obsIndx);
HSobsNA   = HSobsNA(obsIndx);
%WIND
%WINDobsNA = WINDobsNA*obsIndx);

%
%% parameters for data averaging

maxTimeDiffMinutes =30; %minutes
maxTimeDiff = maxTimeDiffMinutes/(24*60); %days
maxDistance = 25; %km radius
minNumberObs = 7;

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
        [LONobs LATobs TIMEobs HSobs indNaN] = bubbleMethodBuoy(LONobsNA, LATobsNA, TIMEobsNA, HSobsNA, buoyTest.lon, buoyTest.lat, buoyTest.time, maxDistance, maxTimeDiff, minNumberObs);
    case 'box'
        
        %% Method 2
        % This creates a gridded cube with demsions lat, lon, & time and
        % the average in each of those cubes. It takes too long (added 30 minutes
        % to a 2 x 2 degree grid (from like 30 seconds).
        
        [LONobs LATobs TIMEobs HSobs] = boxMethodBuoy(LONobsNA, LATobsNA, TIMEobsNA, HSobsNA, buoyTest.lon, buoyTest.lat, maxTimeDiff, minNumberObs);
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




%% save the matched pairs, clear data and loop through the whole domain

clearvars -except TIMEobs LONobs LATobs HSobs HSmd TIMEobsCat...
    LONobsCat LATobsCat HSobsCat HSmdCat mdCol bigLoopLat...
    bigLoopLon averagingMethod mdPath mdFileList loopSize...
    lonLength latLength lonIndexLength latIndexLength fileNum

% save('yearResults2010')
