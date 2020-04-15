% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Script to compare WW3 output to altimeter measurements of significant   %
% wave height and wind speed                                              %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Author:       Clarence Olin Collins III, Ph.D.                          %
% Affiliation:  ERDC - FRF                                                %
% created:      10/28/2019                                                %
% version:      1.0                                                       %
% updates:                                                                %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

%
% Purpose:
% This code design to take data output from the model WAVEWATCHIII and
% compare it with nearby, contemporaneous satelitte altimeter estimations
% of significant wave height and wind speed

%% house cleaning
tic;
close all
clear

%% Paths

codePath = '/Users/tripp/D/Analysis/altimeterComparison/altVmodel';
path(codePath, path); %add code path to search path
mdPath = '/Users/tripp/D/Analysis/altimeterComparison/modelData/AtlanticYearRun/';

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
% method, is currently about 10 times faster for gridded model data.

% averagingMethod = 'bubble';
averagingMethod = 'box';
% averagingMethod = 'none';


%% parameters for data averaging

maxTimeDiffMinutes =30; %minutes
maxTimeDiff = maxTimeDiffMinutes/(24*60); %days
maxDistance = 25; %km
minNumberObs = 7;

%% TO DO: coastline data
% It will be useful to know the proximity of the comparison point to the
% coast where the general quality of altimeter data degrades

%% load WW3 data

% mdPath = 'D:\Analysis\AltimerComparison\pacificYearRun\';
% mdPath = '/Users/tripp/D/Analysis/AltimerComparison/AtlanticYearRun/';
mdFileList  = dir(mdPath);
mdFileList = rmfield (mdFileList,{'date','bytes','isdir','datenum'});
mdFileList = mdFileList(3:end);

% we are going to loop through each file and concatenate the results
TIMEobsCat = [];
LONobsCat  = [];
LATobsCat  = [];
HSobsCat   = [];
HSmdCat    = [];

% WIND
% UmdCat = [];
% VmdCat = [];

loopCount = 0;
for fileNum=1:length(mdFileList)
    
    
    mdCol.time = ncread([mdPath mdFileList(fileNum).name], 'time') + datenum([1990 01 01 00 00 00]);
    mdCol.lat = ncread([mdPath mdFileList(fileNum).name], 'latitude');
    mdCol.lat = double(mdCol.lat);
    mdCol.lonO = ncread([mdPath mdFileList(fileNum).name], 'longitude');
    mdCol.lon = mod(double(mdCol.lonO),360); %convert to degrees east
    mdCol.hs = ncread([mdPath mdFileList(fileNum).name], 'hs');
    mdCol.gridStatus = ncread([mdPath mdFileList(fileNum).name], 'MAPSTA');
    
    % WIND
    % calculate wind from u and v components
    % uwnd = ncread([mdPath 'ww3.201010.nc'], 'uwnd');
    % vwnd = ncread([mdPath 'ww3.201010.nc'], 'vwnd');
    % mdCol.wind = sqrt(uwnd.^2 + vwnd.^2);
    
    
    %% reduce the model domain so it can run on my desktop
    % the whole grid runs for Pacific: lon 110:0.5:300, lat -64:0.5:64; 190 x 128
    
    loopSize = 5; %split grid into a 5x5 chunks
    
    lonLength = length(mdCol.lon);
    latLength = length(mdCol.lat);
    
    lonIndexLength = round(lonLength/loopSize);
    latIndexLength = round(latLength/loopSize);
    
    
    %loop over each grid chunk
    for bigLoopLon = 1:loopSize;
        for bigLoopLat = 1:loopSize;
            loopCount = loopCount + 1;
            
            if bigLoopLon < loopSize
                lonInd =1+lonIndexLength*(bigLoopLon-1):lonIndexLength*(bigLoopLon);
            elseif bigLoopLon == loopSize
                lonInd =1+lonIndexLength*(bigLoopLon-1):length(mdCol.lon);
            else
                disp('something wierd has happened')
            end
            
            if bigLoopLat < loopSize
                latInd =1+latIndexLength*(bigLoopLat-1):latIndexLength*(bigLoopLat);
            elseif bigLoopLat == loopSize
                latInd =1+latIndexLength*(bigLoopLat-1):length(mdCol.lat);
            else
                disp('something wierd has happened')
            end
            
            %build a new structure for the chunk
            mdTest.lonO  = mdCol.lonO(lonInd); %keep original lon for plotting
            mdTest.lon  = mdCol.lon(lonInd);   %convert to deg. E for comparing
            
            % does the grid cross the prime merdian?
            if sum([ismember(mdTest.lon,359.5); ismember(mdTest.lon,0)]) == 2
                crossesPrime = 1;
            else
                crossesPrime = 0;
            end
            
            
            mdTest.lat  = mdCol.lat(latInd);
            timeInd = 1:length(mdCol.time);
            mdTest.time = mdCol.time(timeInd);
            mdTest.hs   = mdCol.hs(lonInd,latInd,timeInd);
            mdTest.gridStatus = mdCol.gridStatus(lonInd,latInd);
            
            %if the grid status indicates no active points, skip this loop
            if sum(sum(mdTest.gridStatus)) == 0
                continue
            end
            
            % test plot
            %             figure
            %             pcolor(mdTest.lonO,mdTest.lat,mdTest.gridStatus')
            %             shading('interp')
            %             colorbar
            
            %% load altimeter data, version 1
            % using the Ribal and Young database
            
            loadSatList = defineSatList(mdTest.time,altPath);
            
            % use model lat - lon to narrow down to files loaded, and load data from
            % local netCDF files
            
            obs = getAltimeterObsWW3(loadSatList, mdTest.lat, mdTest.lon, mdTest.lonO, altPath, crossesPrime);
            
            % if there are no obs, skip this loop
            
            if ~exist('obs','var')
                continue
                disp('obs don''t exist')
            else
                disp('obs exist')
            end
            %% quality control
            % flags: In the present database, a series of data flags defined as 1, 2,
            % 3, 4, and 9 represent Good_data, Probably_good_data, SAR-mode data or
            % possible hardware error (only used for CRYOSAT-2), Bad_data and
            % Missing_data, respectively, have been used. We will retain only good data
            % for now
            
            % observations <50km from a coast are flagged as probably good.
            % So to access these, one would need to also allow probably
            % good data
            
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
            % now you have all the data in one place
            
            % concatenate data
            LONobsNA  = vertcat(obs(:).lon);
            LATobsNA  = vertcat(obs(:).lat);
            TIMEobsNA = vertcat(obs(:).time);
            HSobsNA   = vertcat(obs(:).hsKcal);
            %WIND
            % WINDobsNA = vertcat(obs(:).windCal);
                        
            [obsIndx dum2] = find(TIMEobsNA >= min(mdTest.time) - maxTimeDiff...
                & TIMEobsNA <= max(mdTest.time) + maxTimeDiff);
            
            LONobsNA  = LONobsNA(obsIndx);
            LATobsNA  = LATobsNA(obsIndx);
            TIMEobsNA = TIMEobsNA(obsIndx);
            HSobsNA   = HSobsNA(obsIndx);
            %WIND
            %WINDobsNA = WINDobsNA*obsIndx);
            
            %
            
            
            %% Reducing data by averaging over space and time
            % METHOD I
            % this method takes the x-y-time grid from the model, then finds all the
            % data within a maxDistance radius and within maxTimeDiff. This should form
            % a bubble around the data and we get the average in that bubble.
            
            % There may be a bug somewhere in the code because it never compares as well
            % as the gridded method, which takes data from much further away.
            
            % to do: add gaussian weighted average
            switch averagingMethod
                case 'bubble'
                    [LONobs LATobs TIMEobs HSobs] = bubbleMethod(LONobsNA,...
                        LATobsNA, TIMEobsNA, HSobsNA, mdTest.lon, mdTest.lat,...
                        mdTest.time, maxDistance, maxTimeDiff, mdTest.gridStatus,...
                        minNumberObs);
                case 'box'
                    
                    %% Method 2
                    % This creates a gridded cube with demsions lat, lon, & time and
                    % the average in each of those cubes. It takes too long (added 30 minutes
                    % to a 2 x 2 degree grid (from like 30 seconds).
                    
                    [LONobs LATobs TIMEobs HSobs] = boxMethod(LONobsNA,...
                        LATobsNA, TIMEobsNA, HSobsNA, mdTest.lon, mdTest.lat,...
                        maxTimeDiff, mdTest.gridStatus, minNumberObs);
                case 'none'
                    
                    LONobs  = LONobsNA;
                    LATobs  = LATobsNA;
                    TIMEobs = TIMEobsNA;
                    HSobs   = HSobsNA;
                    clear LONobsNA LATobsNA TIMEobsNA HSobsNA
            end
            %}
            
            %% PART ?: collocation
            
            % phase 1: interpolation
            
            [LONmd LATmd TIMEmd] = ndgrid(mdTest.lon,mdTest.lat,mdTest.time);
            
            % LONobs  = vertcat(obs(:).lon);
            % LATobs  = vertcat(obs(:).lat);
            % TIMEobs = vertcat(obs(:).time);
            % HSobs   = vertcat(obs(:).hsKcal);
            
            F = scatteredInterpolant(double(LONmd(:)), double(LATmd(:)), TIMEmd(:), mdTest.hs(:),'linear','none'); % can't run on full domain
            HSmd = F(LONobs, LATobs, TIMEobs);
            
            indNaN = ~isnan(HSmd);
            HSmd = HSmd(indNaN);
            
            HSobs  = HSobs(indNaN);
            LONobs = LONobs(indNaN);
            LATobs = LATobs(indNaN);
            TIMEobs = TIMEobs(indNaN);
            
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
            disp(['loop ' num2str(loopCount) ' of ' num2str(loopSize^2*length(mdFileList))])
            if loopCount == 1
            disp(['estimated run time: ~'...
                num2str(toc/loopCount*loopSize^2*length(mdFileList)/3600) ' hours'])
            end
            % toc
            
            
            %% save the matched pairs, clear data and loop through the whole domain
            
            clearvars -except TIMEobs LONobs LATobs HSobs HSmd TIMEobsCat...
                LONobsCat LATobsCat HSobsCat HSmdCat mdCol bigLoopLat...
                bigLoopLon averagingMethod mdPath mdFileList loopSize...
                lonLength latLength lonIndexLength latIndexLength fileNum...
                mdPath altPath codePath glyph loopCount maxTimeDiff...
                maxDistance minNumberObs loopCount savePath

            
            TIMEobsCat = [TIMEobsCat; TIMEobs];
            LONobsCat  = [LONobsCat; LONobs];
            LATobsCat  = [LATobsCat; LATobs];
            HSobsCat   = [HSobsCat; HSobs];
            HSmdCat    = [HSmdCat; HSmd];
            
        end
    end
end

save(savePath)
