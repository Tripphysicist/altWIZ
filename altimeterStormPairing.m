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

%% load Single Storm data or single year data
if numel(stormYear) == 1
    try
        stormData = getStormTrackTCOBS(stormYear,stormName);
        disp(['in TC-OBS database']);
    catch
        disp(['couldn''t find ' stormName ', in TC-OBS database trying IBTrACS']);
        try
            stormData = getStormTrackIBTrACS(stormYear,stormName);
            disp(['in IBTrACS database']);
        catch
            disp(['couldn''t find ' stormName ', in IBTrACS either'...
                ' maybe try all upper case']);
        end
    end
    
    
    %% special case of all altimeter data
    % I want to take TCOBS data when possible, but fill in the rest with
    % IBTrACS because it covers a longer time period. Of the area of
    % crossover, there is one storm in IBTrACS that is not in TCOBS -
    % Barbara (the first storm in 2013), most likely because it was a Pacific Basin Storm. Becuase we
    % are not currently considering Pacific storms, we simply remove
    % Barbara (entry 439 when considering 1985+ data), this is semi-hard coded
    % at the moment, should improve this implementation in the future
elseif numel(stormYear) == 36 & stormYear(1) == 1985 & strcmp(stormName,'all')
    disp('!! special case looking at all atlimeter data !!')
    stormDataIB = getStormTrackIBTrACS(stormYear,stormName);
    stormDataTC = getStormTrackTCOBS(stormYear,stormName);
    % fields = fieldnames(stormDataIB);
    [index13, ~] = find(vertcat(stormDataIB(:).year)==2013);
    stormData=[stormDataIB(1:index13(1)-1) stormDataIB(index13(1)+1:end)];
    [index89, ~] = find(vertcat(stormDataIB(:).year)==1989);
    stormData(index89(1):index89(1)+numel(stormDataTC)-1)=stormDataTC;
end

%% reshape HURDAT2 data
% no longer using HURDAT2
% count = 0;
% fields = fieldnames(stormData);
% for i =1:length(stormData.time)
%     if stormData.index(i) == 1
%         count = count + 1;
%     end
%     for jj = 1:numel(fields)
%     stormData2(count).(fields{jj})(stormData.index(i)) = stormData.(fields{jj})(i);
%     end
% end


%% subsample
%
%interpolate to maxTimeDifference into a new structure stormI
fields = fieldnames(stormData);
for ii=1:numel(stormData)
    stormData(ii).lon(stormData(ii).lon<0) = stormData(ii).lon(stormData(ii).lon<0) + 360; %0 - 360
    originalTime = stormData(ii).time;
    subsampleTime = originalTime(1):maxTimeDiff*2:originalTime(end);
    if numel(originalTime) == numel(subsampleTime)
        for j = 1:length(fields)
            stormI(ii).(fields{j}) = stormData(ii).(fields{j});
        end
        continue
    end
    for j = 1:length(fields)
        if strcmp((fields{j}),'year') | strcmp((fields{j}),'name')
            stormI(ii).(fields{j}) = stormData(ii).(fields{j});
            continue
        end
        if sum(isnan(stormData(ii).(fields{j}))) == length(stormData(ii).(fields{j})) %no information to interp
            stormI(ii).(fields{j}) = NaN(1,length(subsampleTime));
        else
            % switched this from interp1 to interp1gap (file exchange), to handle NaNs, max gap 6 hours
            if strcmp((fields{j}),'dir') %handling directional data interp
                cosTemp  = cosd(stormData(ii).(fields{j}));
                sinTemp  = sind(stormData(ii).(fields{j}));
                cosTempI = interp1gap(originalTime, cosTemp, subsampleTime,.5);
                sinTempI = interp1gap(originalTime, sinTemp, subsampleTime,.5);
                stormI(ii).(fields{j}) = mod((180/pi)*atan2(sinTempI,cosTempI),360);
                clear cosTemp sinTemp cosTempi sinTempi
            else
                stormI(ii).(fields{j}) = interp1gap(originalTime, stormData(ii).(fields{j}), subsampleTime,.5);
            end
        end
    end
end

% make sure all are vectors are oriented the same way
for i =1:numel(stormI)
    [~, m] = size(stormI(i).time);
    if m ~= 1
        for j = 1:length(fields)
            if strcmp((fields{j}),'year') | strcmp((fields{j}),'name')
                continue
            end
            stormI(i).(fields{j}) = transpose(stormI(i).(fields{j}));
        end
    end
end

%%
count = 0;
for i=1:length(stormI)
    clear loadSatList
    % use storm time range to exclude some satellite missions
    if i == 1
        switch options.altDatabase
            case 'RY19'
                [loadSatList, satList] = defineSatListRY19(stormI(i).time,altPath);
            case 'ESA'
                [loadSatList, satList] = defineSatListESA(stormI(i).time,altPath);
        end
    else
        satListLength = 0;
        for j = 1:length(satList)
            satSpan = floor(satList(j).timeStart):ceil(satList(j).timeEnd);
            mdSpan  = floor(min(stormI(i).time)):ceil(max(stormI(i).time));
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
        disp(['no satellites during ' stormI(i).name ' ' num2str(stormI(i).year)]);
        continue
    end
    
    % use model lat - lon to narrow down to files loaded, and load data from
    % local netCDF files
    
    switch options.altDatabase
        case 'RY19'
            obs = getRY19AltimeterObsBuoy(loadSatList, stormI(i).lat, stormI(i).lon, altPath, options.QC);
        case 'ESA'
            obs = getESAAltimeterObs(loadSatList, stormI(i).time, altPath, options.QC);
    end
    
    % if there are no obs, skip this loop
    if isempty(obs)
        disp(['no altimeter observations for storm ' stormI(i).name ' ' num2str(stormI(i).year)]);
        continue
    else
        disp(['altimeter observations exist for ' stormI(i).name ' ' num2str(stormI(i).year)])
    end
    
    
    % reduce data based on lon, lat, time
    
    % concatenate data
    LONobs  = vertcat(obs(:).lon);
    % make sure this is 0 - 360
    LONobs(LONobs<0) = LONobs(LONobs<0)+360; %0 - 360
    LATobs  = vertcat(obs(:).lat);
    TIMEobs = vertcat(obs(:).time);
    HSobs   = vertcat(obs(:).hs);
    WINDobs   = vertcat(obs(:).wind);
    % HSERobs   = vertcat(obs(:).hsEr); % error estimate
    % HSQCobs   = vertcat(obs(:).hsQC); % QC flag
    % SATIDobsNA = vertcat(obs(:).satID)); % Satellite Mission ID
    for idx = 1:length(stormI(i).time)
        distance = latlon2dist(LATobs,LONobs,stormI(i).lat(idx),stormI(i).lon(idx));
        [nearStormIndex , ~] = find(distance <= maxDistance &...
            abs(stormI(i).time(idx) - TIMEobs) <= maxTimeDiff/2);
        count = count + 1;
        stormObs(count).lat   = LATobs(nearStormIndex);
        stormObs(count).lon   = LONobs(nearStormIndex);
        stormObs(count).time  = TIMEobs(nearStormIndex);
        stormObs(count).hs    = HSobs(nearStormIndex);
        stormObs(count).wind  = WINDobs(nearStormIndex);
        
        %         plot(stormObs(count).time,stormObs(count).altHs,'.')
        %         hold on
        %         pause
        %         shg
        %         datetick('x')
    end
    %interpolat storm information to each pass
end

%pair data / interpolate

% get rid of NaNs

% indNaNobs = ~isnan(stormObs.hs);
% stormObs.lat = stormObs.lat(indNaNobs);
% stormObs.lon = stormObs.lon(indNaNobs);
% stormObs.time = stormObs.time(indNaNobs);
% stormObs.hs = stormObs.hs(indNaNobs);
% stormObs.wind = stormObs.wind(indNaNobs);

% stormObs.options = options;
disp('great success!')
if options.save
    try
        save(savePath,'stormObs','stormI','options')
    catch
        save(stormData,'stormObs','stormI','options')
    end
end