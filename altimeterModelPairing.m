function [pData] = altimeterModelPairing(codePath,mdPath,altPath,savePath,options)
% [pData] = altimeterModelPairing(codePath,buoyPath,altPath,savePath,options)
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% % Script to compare wave model output to altimeter measurements of .    % 
% % significant wave height .                                             %
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
% mdPath = '/Users/tripp/D/Analysis/altimeterComparison/modelData/AtlanticYearRun/';
% altPath = '/Users/tripp/D/Datasets/Satellite/Altimeter/Ribal_Young_2019/'; % path to Ribal and Young data
% savePath = '/Users/tripp/D/Analysis/altimeterComparison/data/atlantic2016May-2017Nov';
% options.altDatabase     = 'RY19'; % or 'ESA'
% options.averagingMethod =  1; 
% options.maxTimeDiff     =  30/(24*60);
% options.maxDistance     =  25;
% options.minNumberObs    =  7;
% options.loopSize        =  5;
% options.QC              =  1;
% options.save            =  0;
% [pData] = altimeterModelPairing(codePath,mdPath,altPath,savePath,options)
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% % Author:       Clarence Olin Collins III, Ph.D.                        %
% % Affiliation:  ERDC - FRF                                              %
% % created:      06/01/2020                                              %
% % version:      1.1                                                     %
% % updates:                                                              %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% TO DO:
% (1) make it more robust to multiple different grid resolutions, right
%     now it works great if the grid resolution is 0.5 degrees.
% (2) implement kd tree for interpolation
% (3) add a Gaussian weighted average
% Purpose:
% This code design to take data output from the model WAVEWATCHIII and
% compare it with nearby, contemporaneous satelitte altimeter estimations
% of significant wave height and wind speed

tic
%% Paths
path(codePath, path); %add code path to search path

% This code uses the Ribal and Young (2019) dataset. You can find it here:
% https://www.nature.com/articles/s41597-019-0083-9
% The dataset is large, ~100Gb.

%% Global Parameters
% There are two main methods for reducing altimeter data, which generally
% outputs data at 1 Hz. One way to do it is to use "find" and set a maximum
% distance and time. I vizualize this as a bubble and hence it is dubbed
% the bubble method. Another way to do it is to discritize the the grid and
% mean everything in each bin or cube. This I call the box method. The box
% method, is currently about 10 times faster for gridded model data.

switch options.averagingMethod
    case 1
        averagingMethod = 'box';
    case 2
        averagingMethod = 'bubble';
    case 3
        averagingMethod = 'none';
end

%% parameters for data averaging

maxTimeDiff = options.maxTimeDiff; % [days]
maxDistance = options.maxDistance; % [km] radius
minNumberObs = options.minNumberObs;%

%% load WW3 data

%Build a list of model files
mdFileList  = dir(mdPath);
mdFileList = rmfield (mdFileList,{'date','bytes','isdir','datenum'});
mdFileList = mdFileList(3:end);

%we are going to loop through each file and concatenate the results
TIMEobsCat    = [];
LONobsCat     = [];
LATobsCat     = [];
HSobsCat      = [];
HSobsStdCat   = [];
HSobsNoSamCat = [];
HSmdCat       = [];

%WIND
% UmdCat = [];
% VmdCat = [];

loopCount = 0;
for fileNum=1:length(mdFileList)
    %load model data, try a couple of variable naming conventions
    
    %update this section to look at the info file for variable names..
    % modelInfo = ncinfo([mdPath mdFileList(fileNum).name]);
    
    try
        mdCol.lat = ncread([mdPath mdFileList(fileNum).name], 'latitude');
    catch 
        mdCol.lat = ncread([mdPath mdFileList(fileNum).name], 'lat');
        SACSfile
    end
    mdCol.lat = double(mdCol.lat);
    try
        mdCol.lonO = ncread([mdPath mdFileList(fileNum).name], 'longitude');
    catch
        mdCol.lonO = ncread([mdPath mdFileList(fileNum).name], 'lon');
    end
    mdCol.lon = mod(double(mdCol.lonO),360); % 0 - 360 degrees
    mdCol.lonO = mdCol.lon;
    mdCol.lonO(mdCol.lon>180) = mdCol.lon(mdCol.lon>180)-360; %-180 - 180
    
    
    try
        mdCol.hs = ncread([mdPath mdFileList(fileNum).name], 'hs');
    catch
        mdCol.hs = ncread([mdPath mdFileList(fileNum).name], 'waveHs');
    end
    try
        mdCol.gridStatus = ncread([mdPath mdFileList(fileNum).name], 'MAPSTA');
    catch
        mdCol.gridStatus = ncread([mdPath mdFileList(fileNum).name], 'mask');
    end
    %becareful here, the code is not robust to different time stamps, may
    %need to make some changes yourself
    if exist('SACSfile')
        timeInSeconds = ncread([mdPath mdFileList(fileNum).name], 'time');
        mdCol.time = double(timeInSeconds)/(60*60*24) + datenum([1970 01 01 00 00 00]);
    else
        mdCol.time = ncread([mdPath mdFileList(fileNum).name], 'time') + datenum([1990 01 01 00 00 00]);
    end
    
    %WIND
    %calculate wind from u and v components
    % uwnd = ncread([mdPath 'ww3.201010.nc'], 'uwnd');
    % vwnd = ncread([mdPath 'ww3.201010.nc'], 'vwnd');
    % mdCol.wind = sqrt(uwnd.^2 + vwnd.^2);
    
    
    %% reduce the model domain so it can run on my desktop
    
    loopSize = options.loopSize; %loopSize = N, split grid into a NxN chunks
    
    lonLength = length(mdCol.lon);
    latLength = length(mdCol.lat);
    
    lonIndexLength = round(lonLength/loopSize);
    latIndexLength = round(latLength/loopSize);
    
    
    %loop over each grid chunk
    for bigLoopLon = 1:loopSize
        for bigLoopLat = 1:loopSize
            loopCount = loopCount + 1;
            %deal with last loop if not evenly divided
            if bigLoopLon < loopSize
                lonInd =1+lonIndexLength*(bigLoopLon-1):lonIndexLength*(bigLoopLon);
                %needs to be overlapping by 1, dont know why
                lonInd = [lonInd lonInd(end)+1];
            elseif bigLoopLon == loopSize
                lonInd =1+lonIndexLength*(bigLoopLon-1):length(mdCol.lon);
            else
                disp('something wierd has happened')
            end
            %deal with last loop if not evenly divided
            if bigLoopLat < loopSize
                latInd =1+latIndexLength*(bigLoopLat-1):latIndexLength*(bigLoopLat);
                %needs to be overlapping by 1, dont know why
                latInd = [latInd latInd(end)+1];
            elseif bigLoopLat == loopSize
                latInd =1+latIndexLength*(bigLoopLat-1):length(mdCol.lat);
            else
                disp('something wierd has happened')
            end
            %build a new structure for the chunk
            mdTest.lonO  = round(mdCol.lonO(lonInd),1); %keep original lon for plotting
            mdTest.lon  = round(mdCol.lon(lonInd),1);   %convert to deg. E for comparing
            
            % does the grid cross the prime merdian?
            % this needs to be robust to different resolutions
            if sum([ismember(mdTest.lon,359.5); ismember(mdTest.lon,0)]) == 2
                crossesPrime = 1;
            else
                crossesPrime = 0;
            end
            mdTest.lat  = round(mdCol.lat(latInd),1);
            timeInd = 1:length(mdCol.time);
            mdTest.time = mdCol.time(timeInd);
            mdTest.hs   = mdCol.hs(lonInd,latInd,timeInd);
            mdTest.gridStatus = mdCol.gridStatus(lonInd,latInd);
            
            %if the grid status indicates no active points, skip this loop
            if sum(sum(mdTest.gridStatus)) == 0
                disp(['no active grid points for loop ' num2str(loopCount)...
                    ' of ' num2str(loopSize^2*length(mdFileList))]);
                continue
            end
            
            %test plot
            %                         figure
            %                         pcolor(mdTest.lonO,mdTest.lat,mdTest.gridStatus')
            %                         shading('interp')
            %                         colorbar
            
            %% load altimeter data, using the Ribal and Young database
            
            %use model time range to exclude some satellite missions
            switch options.altDatabase
                case 'RY19'
                loadSatList = defineSatListRY19(mdTest.time,altPath);
                case 'ESA'
                loadSatList = defineSatListESA(mdTest.time,altPath);
            end
            
            if isempty(loadSatList)
                disp(['no data in database for the time period ' num2str(mdTest.time(1))...
                    ' to ' num2str(mdTest.time(end)) ' during ' num2str(loopCount)...
                    ' of ' num2str(loopSize^2*length(mdFileList))]);
                continue
            else
%                disp('there were satellites during this time')
            end

            
            % use model lat - lon to narrow down to files loaded, and load data from
            % local netCDF files
            switch options.altDatabase
                case 'RY19'
                obs = getRY19AltimeterObsModel(loadSatList, mdTest.lat, mdTest.lon, mdTest.lonO, altPath, crossesPrime, options.QC);
                case 'ESA'
                obs = getESAAltimeterObsModel(loadSatList, mdTest.time, altPath, options.QC);
            end
            
            % if there are no obs, skip this loop
            if isempty(obs)
                disp(['no altimeter observations for loop ' num2str(loopCount)...
                    ' of ' num2str(loopSize^2*length(mdFileList))]);
                continue
            else
                disp('raw altimeter observations exist')
            end
            
            
            %% reduce data based on lon, lat, time
            
            % concatenate data
            LONobsNA  = vertcat(obs(:).lon);
            % make sure this is 0 - 360
            LONobsNA(LONobsNA<0) = LONobsNA(LONobsNA<0)+360; %0 - 360
            LATobsNA  = vertcat(obs(:).lat);
            TIMEobsNA = vertcat(obs(:).time);
            HSobsNA   = vertcat(obs(:).hs);
            % HSERobs   = vertcat(obs(:).hsEr); % error estimate
            % HSQCobs   = vertcat(obs(:).hsQC); % QC flag
            % SATIDobsNA = vertcat(obs(:).satID)); % Satellite Mission ID
            %WIND
            % WINDobsNA = vertcat(obs(:).wind);
            
            % NEED TO CHECK LON for consistency (180 or 360)
            [obsIndx , ~] = find(TIMEobsNA >= min(mdTest.time) - maxTimeDiff...
                                & TIMEobsNA <= max(mdTest.time) + maxTimeDiff...
                                & LONobsNA  >= min(mdTest.lon) - 1 ...
                                & LONobsNA  <= max(mdTest.lon) + 1 ...
                                & LATobsNA  >= min(mdTest.lat) - 1 ...
                                & LATobsNA  <= max(mdTest.lat) + 1);
            
            LONobsNA  = LONobsNA(obsIndx);
            LATobsNA  = LATobsNA(obsIndx);
            TIMEobsNA = TIMEobsNA(obsIndx);
            HSobsNA   = HSobsNA(obsIndx);
            
            if isempty(HSobsNA)
                disp(['no altimeter observations during querry time, for loop '...
                    num2str(loopCount) ' of ' num2str(loopSize^2*length(mdFileList))]);
                continue
            else
                %                disp('raw altimeter observations exist')
            end
            
            
            %WIND
            %WINDobsNA = WINDobsNA(obsIndx);
            
            %% Reduce data by averaging over space and time
            switch averagingMethod
                case 'bubble'
                    [LONobs , LATobs , TIMEobs , HSobs , HSobsStd , HSobsNoSam] = ...
                        bubbleMethod(LONobsNA, LATobsNA, TIMEobsNA,...
                        HSobsNA, mdTest.lon, mdTest.lat, mdTest.time,...
                        maxDistance, maxTimeDiff, mdTest.gridStatus,...
                        minNumberObs);
                case 'box'
                    % This creates a gridded cube with demsions lat, lon, & time and
                    % the average in each of those cubes. It takes too long (added 30 minutes
                    % to a 2 x 2 degree grid (from like 30 seconds).
                    [LONobs , LATobs , TIMEobs , HSobs , HSobsStd , HSobsNoSam] = ...
                        boxMethod(LONobsNA, LATobsNA, TIMEobsNA, HSobsNA,...
                        mdTest.lon, mdTest.lat, maxTimeDiff, mdTest.gridStatus,...
                        minNumberObs);
                    
                case 'none'
                    
                    LONobs  = LONobsNA;
                    LATobs  = LATobsNA;
                    TIMEobs = TIMEobsNA;
                    HSobs   = HSobsNA;
                    
                    clear LONobsNA LATobsNA TIMEobsNA HSobsNA
                    
            end
            
            
            %% Collocation - interpolation
            % model -> altimeter observations
            % I had tons of issues with this and tried several different
            % methods. Currently, the longitudes, on a -180 - 180 deg.,
            % grid are shifted up so that lon >= 0, then interpolated. 
            
            if crossesPrime
                lonOffset = abs(min(mdTest.lonO));
                mdNewLon = double(mdTest.lonO + lonOffset);
                [LONmd , LATmd , TIMEmd] = ndgrid(mdNewLon,mdTest.lat,mdTest.time);
                Fmodel = scatteredInterpolant(double(LONmd(:)), double(LATmd(:))...
                    , TIMEmd(:), mdTest.hs(:),'linear','none'); % can't run on full domain
                obsNewLon = double(deg180(LONobs) + lonOffset);
                HSmd = Fmodel(obsNewLon, LATobs, TIMEobs);
                indNaN = ~isnan(HSmd);
                HSmd = HSmd(indNaN);
                
            else
                [LONmd , LATmd , TIMEmd] = ndgrid(mdTest.lon,mdTest.lat,mdTest.time);
                Fmodel = scatteredInterpolant(double(LONmd(:)), double(LATmd(:))...
                    , TIMEmd(:), mdTest.hs(:),'linear','none'); % can't run on full domain
                HSmd = Fmodel(LONobs, LATobs, TIMEobs);
                indNaN = ~isnan(HSmd);
                HSmd = HSmd(indNaN);
                
            end
            %%
            HSobs  = HSobs(indNaN);
            
            if length(HSmd) ~= length(HSobs)
                disp('Error: model length does not match observations');
            end
            
            LONobs = LONobs(indNaN);
            LATobs = LATobs(indNaN);
            TIMEobs = TIMEobs(indNaN);
            if ~strcmp(averagingMethod,'none')
                HSobsStd = HSobsStd(indNaN);
                HSobsNoSam = HSobsNoSam(indNaN);
                disp([num2str(length(HSobsNA)) ' raw altimeter measurments have'...
                    ' been reduced to ' num2str(length(HSobs)) ' averages ']);
            else
                indNaNobs = ~isnan(HSobs);
                HSmd = HSmd(indNaNobs);
                HSobs = HSobs(indNaNobs);
                LONobs = LONobs(indNaNobs);
                LATobs = LATobs(indNaNobs);
                TIMEobs = TIMEobs(indNaNobs);
            end
            
            %keep track of loops
            disp(['loop ' num2str(loopCount) ' of ' num2str(loopSize^2*length(mdFileList))])
            if loopCount == 1
                disp(['the first loop took ~' num2str(toc/60) ' minutes'])
            end
            %% save the matched pairs, clear data and loop through the whole domain
            
            clearvars -except TIMEobs LONobs LATobs HSobs HSmd TIMEobsCat...
                LONobsCat LATobsCat HSobsCat HSmdCat mdCol bigLoopLat...
                bigLoopLon averagingMethod mdPath mdFileList loopSize...
                lonLength latLength lonIndexLength latIndexLength fileNum...
                mdPath altPath codePath glyph loopCount maxTimeDiff...
                maxDistance minNumberObs loopCount savePath HSobsStd HSobsStdCat...
                HSobsNoSam HSobsNoSamCat options
            
            TIMEobsCat = [TIMEobsCat; TIMEobs];
            LONobsCat  = [LONobsCat; LONobs];
            LATobsCat  = [LATobsCat; LATobs];
            HSobsCat   = [HSobsCat; HSobs];
            HSmdCat    = [HSmdCat; HSmd];
            if ~strcmp(averagingMethod,'none')
                HSobsStdCat   = [HSobsStdCat; HSobsStd];
                HSobsNoSamCat   = [HSobsNoSamCat; HSobsNoSam];
            end
        end
    end
end

pData.lon = LONobsCat;
pData.lat = LATobsCat;
pData.time = TIMEobsCat;
pData.altHs = HSobsCat;
pData.mdHs = HSmdCat;
pData.options = options;
pData.options.paths.mdPath = mdPath;
pData.options.paths.savePath = savePath;
pData.options.paths.codePath = codePath;
pData.options.paths.altPath = altPath;

if ~strcmp(averagingMethod,'none')
    pData.altHsStd = HSobsStdCat;
    pData.altNoSam = HSobsNoSamCat;
end
if options.save
    save(savePath,'pData','mdCol')
end