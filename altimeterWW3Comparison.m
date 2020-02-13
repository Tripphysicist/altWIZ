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
% original work, and some of their functions are used with little or no
% change.

%% Part I: house cleaning
tic
close all
clear

%%

% averagingMethod = 'bubble';
averagingMethod = 'box';
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

%% PART ?: load WW3 data, version 1
% phase 1: one .nc file ACTIVE
% phase 2: loop through a directory

mdPath = 'D:\Analysis\AltimerComparison\pacificYearRun\';
mdFileList  = dir(mdPath);
mdFileList = rmfield (mdFileList,{'date','bytes','isdir','datenum'});
mdFileList = mdFileList(3:end);


% for fileNum=1:length(mdFileList)


%     mdCol.time = ncread([mdPath mdFileList(fileNum).name], 'time') + datenum([1990 01 01 00 00 00]);
%     mdCol.lat = ncread([mdPath mdFileList(fileNum).name], 'latitude');
%     mdCol.lat = double(mdCol.lat);
%     mdCol.lon = ncread([mdPath mdFileList(fileNum).name], 'longitude');
%     mdCol.lon = double(mdCol.lon);
%     mdCol.hs = ncread([mdPath mdFileList(fileNum).name], 'hs');
%     mdCol.gridStatus = ncread([mdPath mdFileList(fileNum).name], 'MAPSTA');


mdCol.time = ncread([mdPath 'ww3.201010.nc'], 'time') + datenum([1990 01 01 00 00 00]);
mdCol.lat = ncread([mdPath 'ww3.201010.nc'], 'latitude');
mdCol.lat = double(mdCol.lat);
mdCol.lon = ncread([mdPath 'ww3.201010.nc'], 'longitude');
mdCol.lon = double(mdCol.lon);
mdCol.hs = ncread([mdPath 'ww3.201010.nc'], 'hs');
mdCol.gridStatus = ncread([mdPath 'ww3.201010.nc'], 'MAPSTA');

%WIND
% calculate wind from u and v components
% uwnd = ncread([mdPath 'ww3.201010.nc'], 'uwnd');
% vwnd = ncread([mdPath 'ww3.201010.nc'], 'vwnd');
% mdCol.wind = sqrt(uwnd.^2 + vwnd.^2);

% PAPA TEST
%{
%PAPA TEST
mdTest.time = papaTime;
mdTest.lon  = (min(papaLon) - 0.5 : max(papaLon) + 0.5);
mdTest.lat  = min(papaLat)  - 0.5 : max(papaLat) + 0.5;
mdTest.hs   = papaHs;
%}

%% reduce the model domain so it can run on my desktop
% the whole grid runs lon 110:0.5:300, lat -64:0.5:64; 190 x 128


% reduce to 1 degree 40 - 41 N and 200 - 201 East

%quarter of the domain
% around 200 deg lon
% around 40 deg lat

%centerLatInd = 181; %ocean
% centerLonInd =  229; %ocean

% centerLatInd = 6; %ocean
% centerLonInd =  6; %ocean
TIMEobsCat = [];
LONobsCat  = [];
LATobsCat  = [];
HSobsCat   = [];
HSmdCat    = [];

loopSize = 5;

lonLength = length(mdCol.lon);
latLength = length(mdCol.lat);

lonIndexLength = round(lonLength/loopSize);
latIndexLength = round(latLength/loopSize);

for bigLoopLon = 1:loopSize;
    for bigLoopLat = 1:loopSize;
        
        % latInd = centerLatInd - 15 : centerLatInd + 15;
        % lonInd = centerLonInd - 15 : centerLonInd + 15;
        
        % mdTest.lat  = mdCol.lat(latInd);
        % mdTest.lon  = mdCol.lon(lonInd);
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

        
        mdTest.lon  = mdCol.lon(lonInd);    
        mdTest.lat  = mdCol.lat(latInd);
        
%        timeInd = 1:10;
        timeInd = 1:length(mdCol.time);
        mdTest.time = mdCol.time(timeInd);
        
        mdTest.hs   = mdCol.hs(lonInd,latInd,timeInd);
        %whole domain
        % mdTest.lat  = mdCol.lat(:);
        % mdTest.lon  = mdCol.lon(:);
        % mdTest.hs   = mdCol.hs(:,:,:);
        mdTest.gridStatus = mdCol.gridStatus(lonInd,latInd);
        
        % clear mdCol
        
        
        % test plot
        % figure
        % pcolor(mdTest.lon,mdTest.lat,mdTest.gridStatus')
        % shading('interp')
        % colorbar
        
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
            mdSpan  = floor(mdTest.time(1)):ceil(mdTest.time(end));
            
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
                
                % for rectilinear grid, add one to the greater than test because
                % each alt file runs from altNorth to altNorth + 1
                if altNorth + 1 >= mdTest.lat(1) & altNorth <= mdTest.lat(end) & ...
                        altEast + 1 >= mdTest.lon(1) & altEast <= mdTest.lon(end)
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
            %     obs(i).hsKstd = obs(i).hsKstd(qcPassInd );
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
        
        halfHour = 0.0208;
        
        [obsIndx dum2] = find(TIMEobsNA >= min(mdTest.time) - halfHour...
            & TIMEobsNA <= max(mdTest.time) + halfHour);
        
        LONobsNA  = LONobsNA(obsIndx);
        LATobsNA  = LATobsNA(obsIndx);
        TIMEobsNA = TIMEobsNA(obsIndx);
        HSobsNA   = HSobsNA(obsIndx);
        
        %
        %% parameters for data averaging
        
        maxTimeDiffMinutes =30; %minutes
        maxTimeDiff = maxTimeDiffMinutes/(24*60); %days
        maxDistance = 25; %km
        minNumberObs = 7;
        
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
                [LONobs LATobs TIMEobs HSobs] = bubbleMethod(LONobsNA, LATobsNA, TIMEobsNA, HSobsNA, mdTest.lon, mdTest.lat, mdTest.time, maxDistance, maxTimeDiff, mdTest.gridStatus, minNumberObs);
            case 'box'
                
                %% Method 2
                % This creates a gridded cube with demsions lat, lon, & time and
                % the average in each of those cubes. It takes too long (added 30 minutes
                % to a 2 x 2 degree grid (from like 30 seconds).
                
                [LONobs LATobs TIMEobs HSobs] = boxMethod(LONobsNA, LATobsNA, TIMEobsNA, HSobsNA, mdTest.lon, mdTest.lat, maxTimeDiff, mdTest.gridStatus, minNumberObs);
            case 'none'
                
                LONobs  = LONobsNA;
                LATobs  = LATobsNA;
                TIMEobs = TIMEobsNA;
                HSobs   = HSobsNA;
                clear LONobsNA LATobsNA TIMEobsNA HSobsNA
        end
        %}
        
        %% PAPA test
        %{
% LONobs  = vertcat(obs(:).lon);
% LATobs  = vertcat(obs(:).lat);
% TIMEobs = vertcat(obs(:).time);
% HSobs   = vertcat(obs(:).hsKcal);

F = scatteredInterpolant(LONobs,LATobs,TIMEobs,HSobs,'linear','none');
HSalt = F(papaLon,papaLat,papaTime);
indNaN = ~isnan(HSalt);
HSalt = HSalt(indNaN);
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
        
        toc
        
        
        %% save the matched pairs, clear data and loop through the whole domain
        
        clearvars -except TIMEobs LONobs LATobs HSobs HSmd TIMEobsCat...
            LONobsCat LATobsCat HSobsCat HSmdCat mdCol bigLoopLat...
            bigLoopLon averagingMethod mdPath mdFileList loopSize...
            lonLength latLength lonIndexLength latIndexLength
        
        TIMEobsCat = [TIMEobsCat; TIMEobs];
        LONobsCat  = [LONobsCat; LONobs];
        LATobsCat  = [LATobsCat; LATobs];
        HSobsCat   = [HSobsCat; HSobs];
        HSmdCat    = [HSmdCat; HSmd];
    end
end
% end
%% PART ?: calculate statistics

%% PART ?: make simple plots
figure
subplot(1,3,1)
plot(TIMEobsCat,HSobsCat,'.')
hold on
plot(TIMEobsCat,HSmdCat,'.')
grid on
datetick
subplot(1,3,2)
plot(LONobsCat,HSobsCat,'.')
hold on
plot(LONobsCat,HSmdCat,'.')
grid on
subplot(1,3,3)
plot(LATobsCat,HSobsCat,'.')
hold on
plot(LATobsCat,HSmdCat,'.')
grid on

packfig(1,3)

%subplot(1,2,2)
L = Ocstatsp(HSobsCat,HSmdCat,0.1)

% distributions, histograms
% q-q plots
% scatter plot
% taylor diagram

%% map of stats

minLON = min(LONobsCat);
maxLON = max(LONobsCat);
minLAT = min(LATobsCat);
maxLAT = max(LATobsCat);
blockSize = 1; %[degrees]
statGridEdgesLon = minLON:blockSize:maxLON;
statGridEdgesLat = minLAT:blockSize:maxLAT;
statGridCenterLon = statGridEdgesLon(1:end-1) + diff(statGridEdgesLon)/2;
statGridCenterLat = statGridEdgesLat(1:end-1) + diff(statGridEdgesLat)/2;


for i = 1:length(statGridCenterLon)
    for j = 1:length(statGridCenterLat)
        [index dumy] = find(LONobsCat >= statGridEdgesLon(i) & LONobsCat <= statGridEdgesLon(i+1) &...
            LATobsCat >= statGridEdgesLat(j) & LATobsCat <= statGridEdgesLat(j+1));
        dHS = HSobsCat(index)-HSmdCat(index);
        bias(i,j) = mean(dHS);
        rmse(i,j) = sqrt(mean(dHS.^2));
    end
end

figure
subplot(2,1,1)
pcolor(statGridCenterLon,statGridCenterLat,bias')
shading('interp')
subplot(2,1,2)
pcolor(statGridCenterLon,statGridCenterLat,rmse')
shading('interp')
