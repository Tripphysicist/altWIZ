function [pData , rawAlt , buoyData] = altimeterBuoyPairing(codePath,buoyInfo,altPath,savePath,options)
% function [pData rawAlt buoyData] = altimeterBuoyPairing(codePath,buoyPath,altPath,savePath,options)
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% % Script to compare wave buoy measurements to contemperaneous altimeter %
% % measurements.                                                         %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%
%  INPUT
%  Paths - in order path to the directory for code, directory with model
%  files, directory with altimeter altimeter data, and full path and file
%  name for saving results
%
%  OPTIONS structure (not optional)
%  options.altDatabase -> Ribal and Young (2019) or ESA Sea State
%  options.maxTimeDiff -> maxTimeDiffMinutes/(24*60); %maximum time
%                difference [days]
%  options.maxDistance -> maximum space distance [km] (radius for
%                bubble method)
%  options.minNumberObs -> minimum number of observations for averaging
%  options.QC -> quality control level (1 is strictest, 2 will include
%                coastal data if using RY19)
%  options.save -> save output? logical 1 | 0
%
%  OUTPUT
%  pData - structure of averaged wave height data from altimeter and paired
%  with buoy dara on a shared lon, lat, time
%  pData.lon  - longitude in degrees
%  pData.lat  - latitude in degrees
%  pData.time - time stamp
%  pData.altHsMean    - average altimeter significant wave height
%  pData.altHsNearest - nearest altimeter significant wave height
%  pData.altHsStd     - standard deviation from  averaging method
%  pData.altWindMean    - average altimeter wind speed
%  pData.altWindNearest - nearest altimeter wind speed
%  pData.altWindStd     - standard deviation of altimeter wind speed
%  pData.altHsNoSam     - number of samples for each average
%  pData.buoyHs         - buoy significant wave height
%  pData.meta - documentation of options used
%
%  rawAlt - structure of 1 Hz altimeter data
%  rawAlt.lon
%  rawAlt.lat
%  rawAlt.time
%  rawAlt.Hs
%  rawAlt.wind
%
% e.g.
% codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
% buoyInfo = '/Users/tripp/D/Analysis/altimeterComparison/mini-buoys/AO_MWB.mat';
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
%
% To do: given a buoy name and date, look for buoy data locally, if not
% found, find on-line and download into local directories (if wanted) or
% just run from workspace

%tic
%% Paths
path(codePath, path); %add code path to search path

%% parameters for bubble method od altimeter data averaging
% rename into variables for ease of reading code

maxTimeDiff = options.maxTimeDiff;
maxDistance = options.maxDistance;
minNumberObs = options.minNumberObs;
%%
switch buoyInfo.source
    case 'CDIP'
        
        % CDIP
        cdipURL = ['http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/'...
            buoyInfo.name 'p1/' buoyInfo.name 'p1_historic.nc'];
        cdipTimeOffset = datenum([1970 01 01 00 00 00]);
        cdipTime = double(ncread(cdipURL,'waveTime'))/(60*60*24) + cdipTimeOffset;
        cdipHs = ncread(cdipURL,'waveHs');
        cdipGPStime = double(ncread(cdipURL,'gpsTime'))/(60*60*24) + cdipTimeOffset;
        cdipGPSlat = ncread(cdipURL,'gpsLatitude');
        cdipGPSlon = ncread(cdipURL,'gpsLongitude');
        cdipLat = interp1(cdipGPStime, cdipGPSlat, cdipTime);
        cdipLon = interp1(cdipGPStime, cdipGPSlon, cdipTime) + 360; %0 - 360
        cdipLat(1)= cdipLat(2); %interpolation made this a NaN
        cdipLon(1)= cdipLon(2); %interpolation made this a NaN
        
        % remove fill values
        [goodIndex, ~] = find(cdipHs>0);
        buoyData.time = cdipTime(goodIndex);
        buoyData.lon  = cdipLon(goodIndex);
        buoyData.lat  = cdipLat(goodIndex);
        buoyData.hs   = cdipHs(goodIndex);
        %}
                
    case 'NDBC'
        % grab hfile, I have processed hfiles here:
        load('sData_v20200721.mat')
        % this is from the NDBC website, which is considered the "official
        % archive", however there is not meta-data, i.e. lon, lat, buoy
        % hull, etc. One needs to supply this seperately for now, will make
        % this functional in the future.
        
        for i = 1:length(sData)
            currentFile = sData(i).hFile;
            if isempty(currentFile)
                continue
            end
            bIndex(i) = contains(currentFile,buoyInfo.name);
        end
        
        Hs = vertcat(sData(bIndex).hHs);
        windSpeed = vertcat(sData(bIndex).hWindSpeed);
        time = vertcat(sData(bIndex).htime);
        time = time(~isnan(Hs));
        Hs = Hs(~isnan(Hs));
        windSpeed = windSpeed(~isnan(Hs));
        lon = buoyInfo.lon.*ones(length(Hs),1);
        lat = buoyInfo.lat.*ones(length(Hs),1);
        
        buoyData.time = time;
        buoyData.lon  = 360 + lon; %from negative W to E poisitive 0 - 360
        buoyData.lat  = lat;
        buoyData.hs   = Hs;
        buoyData.windSpeed = windSpeed;
                
    case 'SIO'
        
        % SIO MINI BUOYS
        load('/Users/tripp/D/Analysis/altimeterComparison/mini-buoys/AO_MWB.mat')
        
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
        
        buoyData.time = time;
        buoyData.lon  = 360 + lon; %from negative W to E poisitive 0 - 360
        buoyData.lat  = lat;
        buoyData.hs   = Hs;
        
    case '3DMG'
        load('/Users/tripp/Google Drive/Work/Jensen Buoy Work/46029_analysis.mat',...
            'MG','parMG');

        buoyData.time = MG.time; 
        buoyData.lon = MG.lon + 360;
        buoyData.lat = MG.lat;
        buoyData.hs  = parMG.Hm0';
        
    case 'HIPPY'
        
        load('/Users/tripp/Google Drive/Work/Jensen Buoy Work/46029_analysis.mat',...
            'HP','parHP');

        buoyData.time = HP.time; 
        buoyData.lon = HP.lon + 360;
        buoyData.lat = HP.lat;
        buoyData.hs  = parHP.Hm0';
        
    case 'WR'    
        
        load('/Users/tripp/Google Drive/Work/Jensen Buoy Work/46029_analysis.mat',...
            'WR','parWR');

        buoyData.time = WR.time; 
        buoyData.lon = WR.lon + 360;
        buoyData.lat = WR.lat;
        buoyData.hs  = parWR.Hm0';
        
    case 'ITOP'
        % ITOP
        % EASI-N
        % load([buoyPath 'EASI_Nv2.mat'])
        % buoyData.time = EN.var.yday + datenum([2010 01 01 00 00 00]);
        % buoyData.lon  =  126.968.*ones(length(buoyData.time),1);
        % buoyData.lat  =  21.238.*ones(length(buoyData.time),1);
        % buoyData.hs   = EN.par.int.Hm0;
        % clear EN
        
        % EASI-S
        % load([buoyPath 'EASI_Sv5.mat'])
        % buoyData.time = ES.var.yday + datenum([2010 01 01 00 00 00]);
        % buoyData.lon  =  127.258.*ones(length(buoyData.time),1);
        % buoyData.lat  =  19.683.*ones(length(buoyData.time),1);
        % buoyData.hs   = ES.par.int.Hm0;
        % clear ES
end

%% test plot
figure
plot(buoyData.time,buoyData.hs,'.')
datetick('x')
title('Hs Test Plot')
ylabel('Hs [m]')
datetick('x','mm-YY')

%% load altimeter data

switch options.altDatabase
    case 'RY19'
        loadSatList = defineSatListRY19(buoyData.time,altPath);
    case 'ESA'
        loadSatList = defineSatListESA(buoyData.time,altPath);
end

if isempty(loadSatList)
    disp(['no data in database for ' buoyInfo ' from ' num2str(buoyData.time(1))...
        ' to ' num2str(buoyData.time(end))]);
else
    for i =1:length(loadSatList)
        disp([loadSatList(i).name ' was active with buoy ' buoyInfo.name]);
    end
    disp('Searching for data...');
end

% use model lat - lon to narrow down to files loaded, and load data from
% local netCDF files
switch options.altDatabase
    case 'RY19'
        obs = getRY19AltimeterObsBuoy(loadSatList, buoyData.lat, buoyData.lon, altPath, options.QC);
        
        % if there are no obs, skip this loop
        if isempty(obs)
            disp(['no altimeter observations for buoy ' buoyInfo]);
        else
            disp('raw altimeter observations exist')
        end
        
        %% reduce data based on time and grid status
        
        % concatenate altimeter data
        LONobsNA  = vertcat(obs(:).lon);
        LATobsNA  = vertcat(obs(:).lat);
        TIMEobsNA = vertcat(obs(:).time);
        HSobsNA   = vertcat(obs(:).hs);
        WINDobsNA = vertcat(obs(:).wind);
        
        %set ranges for time, lat, and lon, lat and lon within
        
        [obsIndx , ~] = find(TIMEobsNA >= min(buoyData.time) - maxTimeDiff...
            & TIMEobsNA <= max(buoyData.time) + maxTimeDiff...
            & LONobsNA  >= min(buoyData.lon) - 4 ...
            & LONobsNA  <= max(buoyData.lon) + 4 ...
            & LATobsNA  >= min(buoyData.lat) - 4 ...
            & LATobsNA  <= max(buoyData.lat) + 4);
        
        rawAlt.lon  = LONobsNA(obsIndx);
        rawAlt.lat  = LATobsNA(obsIndx);
        rawAlt.time = TIMEobsNA(obsIndx);
        rawAlt.hs   = HSobsNA(obsIndx);
        rawAlt.wind = WINDobsNA(obsIndx);
        
    case 'ESA'
        
        rawAlt.lon = [];
        rawAlt.lat = [];
        rawAlt.time = [];
        rawAlt.hs   = [];
        rawAlt.wind = [];
        
        for monthI=1:floor(length(buoyData.time)/(48*30))
            if monthI == floor(length(buoyData.time)/(48*30))
                monthIndex = (monthI-1)*1440 + 1: length(buoyData.time);
            else
                monthIndex = (monthI-1)*1440 + 1:monthI*1440;
            end
            
            obs = getESAAltimeterObs(loadSatList, buoyData.time(monthIndex), altPath, options.QC);
            
            % if there are no obs, skip this loop
            if isempty(obs)
                disp(['no altimeter observations for buoy ' buoyInfo]);
                continue
            else
                disp(['raw altimeter observations exist for month ' num2str(monthI)...
                    ' of ' num2str(floor(length(buoyData.time)/(48*30)))])
            end
            
            %% reduce data based on time and grid status
            
            % concatenate altimeter data
            LONobsNA  = vertcat(obs(:).lon);
            LONobsNA(LONobsNA<0) = LONobsNA(LONobsNA<0)+360; %0 - 360
            LATobsNA  = vertcat(obs(:).lat);
            TIMEobsNA = vertcat(obs(:).time);
            HSobsNA   = vertcat(obs(:).hs);
            WINDobsNA = vertcat(obs(:).wind);
            
            % set ranges for time, lat, and lon, lat and lon within
            [obsIndx , ~] = find(TIMEobsNA >= min(buoyData.time) - maxTimeDiff...
                & TIMEobsNA <= max(buoyData.time) + maxTimeDiff...
                & LONobsNA  >= min(buoyData.lon) - 4 ...
                & LONobsNA  <= max(buoyData.lon) + 4 ...
                & LATobsNA  >= min(buoyData.lat) - 4 ...
                & LATobsNA  <= max(buoyData.lat) + 4);
            
            rawAlt.lon  = [rawAlt.lon; LONobsNA(obsIndx)];
            rawAlt.lat  = [rawAlt.lat; LATobsNA(obsIndx)];
            rawAlt.time = [rawAlt.time; TIMEobsNA(obsIndx)];
            rawAlt.hs   = [rawAlt.hs; HSobsNA(obsIndx)];
            rawAlt.wind = [rawAlt.wind; WINDobsNA(obsIndx)];
        end
end

%% Reducing data by averaging over space and time
% to do: add gaussian weighted average
[pData.lon , pData.lat , pData.time , pData.altHsMean , pData.altHsStd ,...
    pData.altHsNearest , pData.altWindMean , pData.altWindStd ,...
    pData.altWindNearest , pData.buoyIndNaN] = bubbleMethodBuoy(rawAlt.lon ,...
    rawAlt.lat , rawAlt.time , rawAlt.hs , rawAlt.wind , buoyData.lon, buoyData.lat...
    , buoyData.time , maxDistance, maxTimeDiff , minNumberObs);

pData.altHsNearest = pData.altHsNearest';
pData.altWindNearest = pData.altWindNearest';

% average altimetered data paired with buoy data
pData.buoyHs = buoyData.hs(pData.buoyIndNaN);
pData.buoyWindSpeed = buoyData.windSpeed(pData.buoyIndNaN);
pData.meta.options = options;
pData.meta.buoyInfo = buoyInfo;
pData.meta.paths.savePath = savePath;
pData.meta.paths.codePath = codePath;
pData.meta.paths.altPath = altPath;
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
