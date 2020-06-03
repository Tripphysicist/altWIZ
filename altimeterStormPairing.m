function [stormObs, stormI] = altimeterStormPairing(stormName,stormYear,codePath,altPath,savePath,options)
% function [stormObs, stormI] = altimeterStormPairing(stormName,stormYear,codePath,altPath,savePath,options)
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% % Script to get altimeter data for tropical storms                      %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%
%  OUTPUT
%  pData - structure of wave height data from altimeter and models with
%  shared lon, lat, time
%  .lon - longitude in degrees, .lat - latitude in degrees, .time - time
%  stamp, .altHs output from averaging method, .altHsStd standard deviation
%  from  averaging method, .altHsNoSam - number of samples
%
%  INPUT
%  Paths - in order path to the directory for code, directory with model
%  files, directory with altimeter altimeter data, and full path and file
%  name for saving results
%
%  OPTIONS structure (not optional)
%  options.averagingMethod = 1 = 'box'; 2 = 'bubble'; 3 = 'none'; Box
%  and bubble method are more fully described below. Default for model is
%  box because it is much quicker. Results are similar.
%  options for pairing data
%  options.maxTimeDiff = maxTimeDiffMinutes/(24*60); %maximum time
%                difference [days]
%  options.maxDistance = 25; %maximum space distance [km] (radius for
%                bubble method.
%  options.minNumberObs = 7; %minimum number of observations for
%                averaging
%  generic options
%  options.loopSize = 5; number_of_loops = number_of_files*loopSize^2
%                consider increasing this if you have memory issues
%  options.QC = 1; strictest quality control
%  options.save -> save output? logical 1 | 0
%
% e.g.
% codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
% stormName = 'BONNIE';
% stormYear = 2016;
% altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/Ribal_Young_2019/'; % path to Ribal and Young data
% savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/testing';
% options.altDatabase     = 'RY19'; % or 'ESA'
% options.maxTimeDiff     =  60/(24*60);
% options.maxDistance     =  250;
% options.QC              =  1;
% options.save            =  0;
% [stormObs, stormI] = altimeterStormPairing(stormName,stormYear,codePath,altPath,savePath,options)
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% % Author:       Clarence Olin Collins III, Ph.D.                        %
% % Affiliation:  ERDC - FRF                                              %
% % created:      06/01/2020                                              %
% % version:      1.0                                                     %
% % updates:                                                              %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

%% Paths
path(codePath, path); %add code path to search path



%% parameters for data averaging

maxTimeDiff = options.maxTimeDiff; % [days]
maxDistance = options.maxDistance; % [km] radius

%% load Storm data
try
    stormData = getStormTrack(stormYear,stormName);
catch
    disp(['couldn''t find ' stormName ', try all capital letters']);
end
%% subsample
stormData.lon(stormData.lon<0) = stormData.lon(stormData.lon<0) + 360; %0 - 360
%interpolate to maxTimeDifference
fields = fieldnames(stormData);
oneIndex = find(stormData.index == 1);
for i = 1:length(fields)
    stormI.(fields{i}) = [];
end

if length(oneIndex) == 1
    currentIndex = (1:length(stormData.time));
    originalTime = stormData.time(currentIndex);
    subsampleTime = originalTime(1):maxTimeDiff:originalTime(end);
    for j = 1:length(fields)
        temp = interp1(originalTime, stormData.(fields{j})(currentIndex), subsampleTime);
        stormI.(fields{j}) = [stormI.(fields{j}) temp];
    end
    stormI.stormCount = ones(1,length(temp));
    
else
    stormCount = 0;
    stormI.stormCount = [];
    for i = 1:length(oneIndex)-1
        if i < length(oneIndex)-1
            currentIndex = (oneIndex(i):oneIndex(i+1)-1);
        else
            currentIndex = (oneIndex(i):length(stormData.time));
        end
        if length(currentIndex) == 1
            continue
        end
        originalTime = stormData.time(currentIndex);
        subsampleTime = originalTime(1):maxTimeDiff:originalTime(end);
        for j = 1:length(fields)
            temp = interp1(originalTime, stormData.(fields{j})(currentIndex), subsampleTime);
            stormI.(fields{j}) = [stormI.(fields{j}) temp];
        end
        stormCount = stormCount + 1;
        stormI.stormCount = [stormI.stormCount ones(1,length(temp)).*stormCount];
    end
end

%% load altimeter data, using the Ribal and Young database

for i=1:max(stormI.stormCount)
    index = find(stormI.stormCount == i);
    clear loadSatList
    %use model time range to exclude some satellite missions
    if i == 1
        switch options.altDatabase
            case 'RY19'
                [loadSatList, satList] = defineSatListRY19(stormI.time(index),altPath);
            case 'ESA'
                loadSatList = defineSatListESA(stormI.time(indeX),altPath);
        end
    else
        satListLength = 0;
        for j = 1:length(satList)
            satSpan = floor(satList(j).timeStart):ceil(satList(j).timeEnd);
            mdSpan  = floor(min(stormI.time(index))):ceil(max(stormI.time(index)));
            if  intersect(satSpan,mdSpan)
                satListLength = satListLength +1;
                loadSatList(satListLength) = satList(j);
            end
        end
        if satListLength == 0
            loadSatList = [];
        end
    end
    
    if isempty(loadSatList)
        disp(['no data in database for ' stormName ' from ' num2str(stormI.time(index(1)))...
            ' to ' num2str(stormI.time(index(end)))]);
        continue
    end
    
    
    % use model lat - lon to narrow down to files loaded, and load data from
    % local netCDF files
    
    switch options.altDatabase
        case 'RY19'
            obs = getRY19AltimeterObsBuoy(loadSatList, stormI.lat(index), stormI.lon(index), altPath, options.QC);
        case 'ESA'
            obs = getESAAltimeterObs(loadSatList, stormI.time(index), altPath, options.QC);
    end
    
    % if there are no obs, skip this loop
    if isempty(obs)
        disp(['no altimeter observations for storm ' num2str(i)]);
        continue
    else
        disp('raw altimeter observations exist')
    end
    
    
    % reduce data based on lon, lat, time
    
    % concatenate data
    LONobs  = vertcat(obs(:).lon);
    % make sure this is 0 - 360
    LONobs(LONobs<0) = LONobs(LONobs<0)+360; %0 - 360
    LATobs  = vertcat(obs(:).lat);
    TIMEobs = vertcat(obs(:).time);
    HSobs   = vertcat(obs(:).hs);
    % HSERobs   = vertcat(obs(:).hsEr); % error estimate
    % HSQCobs   = vertcat(obs(:).hsQC); % QC flag
    % SATIDobsNA = vertcat(obs(:).satID)); % Satellite Mission ID
    %WIND
    % WINDobs = vertcat(obs(:).wind);
    
    for idx = 1:length(index)
        distance = latlon2dist(LATobs,LONobs,stormI.lat(index(idx)),stormI.lon(index(idx)));
        [nearStormIndex , ~] = find(distance <= maxDistance &...
            abs(stormI.time(index(idx)) - TIMEobs) <= maxTimeDiff/2);
        stormObs(index(idx)).lat   = LATobs(nearStormIndex);
        stormObs(index(idx)).lon   = LONobs(nearStormIndex);
        stormObs(index(idx)).time  = TIMEobs(nearStormIndex);
        stormObs(index(idx)).hs    = HSobs(nearStormIndex);
        
        %         plot(stormObs(count).time,stormObs(count).altHs,'.')
        %         hold on
        %         pause
        %         shg
        %         datetick('x')
    end
end
%%
%WIND
%WINDobsNA = WINDobsNA(obsIndx);


% get rid of NaNs

% indNaNobs = ~isnan(stormObs.hs);
% stormObs.lat = stormObs.lat(indNaNobs);
% stormObs.lon = stormObs.lon(indNaNobs);
% stormObs.time = stormObs.time(indNaNobs);
% stormObs.hs = stormObs.hs(indNaNobs);

% stormObs.options = options;

if options.save
    save(savePath,'stormObs','stormI','options')
end