function obs = getESAAltimeterObs(loadSatList, inputTime, altPath, QC)
% obs = getESAAltimeterObsModel(loadSatList, inputTime, altPath, QC)
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% % Script to load ESA data in to workspace                               %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%
%  OUTPUT
%  obs - structure of wave height data from altimeter and models with
%  shared lon, lat, time
%  obs.lon - longitude in degrees
%  obs.lat - latitude in degrees
%  obs.time - time (MATLAB datenum)
%  obs.hs - denoised data as detailed in Q&C (2019)
%  obs.hsEr - uncertainty as estimated in Ash (2012)
%  obs.hsQC - Quality control of individual altimeter measurements is
%           undertaken with checks on instrument flags and ancillary
%           variables. As a result, the SWH comes with a quality level
%           provided in the swh_quality variable. Its meaning is defined as
%           follows:
% 0 - undefined the measurement value is not defined or relevant (missing
%       value, etc?), no quality check was applied.
% 1 - bad the measurement was qualified as not usable after quality check.
% 2 - acceptable the measurement may be usable for specific applications
%       only or the quality check could not fully assess if it is a bad or
%       good value (suspect).
% 3 - good the measurement is usable.
%
%  INPUT
%  loadSatList - list of satellites to load
%  inputTime  - model, buoy, or storm time
%  altPath - path to ESA data
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% % Author:       Clarence Olin Collins III, Ph.D.                        %
% % Affiliation:  ERDC - FRF                                              %
% % created:      06/01/2020                                              %
% % version:      1.0                                                     %
% % updates:                                                              %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% references:
% Ellis Ash, 2012. GlobWave Annual Quality Control Report, Phase 2, pp. 43.
% Quilfen, Y., and Chapron, B., 2019. Ocean Surface Wave-Current Signatures
%   From Satellite Altimeter Measurements. Geophysical Research Letters 46,
%   253?261. https://doi.org/10.1029/2018GL081029
% ESA Sea State Product User Guide
% function that gets loads altimeter data into the workspace depending on
% the model time provided. This searches the ESA database.




dataVersion = 'v1.1';
[glyph , ~] = giveGlyph;
count2 = 0;

startDay = floor(inputTime(1));
endDay   = ceil(inputTime(end));
days = startDay:1:endDay;
dates = datevec(days);
dirYear = num2str(dates(:,1));
dirMon0  = dates(:,2);
%add zero to month
addZeroIndex = dirMon0<10;
for i = 1:length(dirMon0)
    if addZeroIndex(i)
        dirMon(i,:) = ['0' num2str(dirMon0(i))];
    else
        dirMon(i,:) = num2str(dirMon0(i));
    end
end
dirDay0 = dates(:,3);
%add zero to month
addZeroIndex = dirDay0<10;
for i = 1:length(dirDay0)
    if addZeroIndex(i)
        dirDay(i,:) = ['0' num2str(dirDay0(i))];
    else
        dirDay(i,:) = num2str(dirDay0(i));
    end
end

for i = 1:length(loadSatList)
    
    for j = 1: length(days);
        
        altFilePath = [altPath loadSatList(i).name glyph dataVersion glyph...
            glyph dirYear(j,:) glyph dirMon(j,:) glyph dirDay(j,:) glyph];
        altFileList = dir(altFilePath);
        altFileList = rmfield (altFileList,{'date','bytes','isdir','datenum'});
        altFileList = altFileList(3:end);
        
        for k = 1:length(altFileList)
            count2 = count2 + 1;
            fileName = [altFilePath altFileList(k).name];
            
            % load data
            %ping the file to make sure its working
            try 
                info   = ncread(fileName);
            catch 
                disp(['error opening ' fileName])
                continue
            end
            
            satTimeOffset = datenum([1981 01 01 00 00 00]); %same for all
            timeTemp = ncread(fileName, 'time')./(24*60*60); %convert from seconds to days
            obs(count2).time = timeTemp + satTimeOffset;
            obs(count2).lat  = ncread(fileName, 'lat');
            obs(count2).lon  = ncread(fileName, 'lon');
%            obs(count2).satID = loadSatList(i).name;
%            obs(count2).hsAd   = ncread(fileName, 'swh_adjusted');
            obs(count2).hs   = ncread(fileName, 'swh_denoised');
%            obs(count2).hsEr = ncread(fileName, 'swh_uncertainty');
            obs(count2).hsQC   = ncread(fileName, 'swh_quality');
            %WIND
            %         obs(count2).wind   = ncread(fileName, 'wind_speed_alt');
        end
    end
end
%% quality control
% Quality control of individual altimeter measurements is undertaken with
% checks on instrument flags and ancillary variables. As a result, the SWH
% comes with a quality level provided in the swh_quality variable. Its
% meaning is defined as follows:
% value meaning description
% 0 - undefined the measurement value is not defined or relevant (missing
%       value, etc?), no quality check was applied.
% 1 - bad the measurement was qualified as not usable after quality check.
% 2 - acceptable the measurement may be usable for specific applications
%       only or the quality check could not fully assess if it is a bad or
%       good value (suspect).
% 3 - good the measurement is usable.

for i = 1:length(obs)
    if QC == 5
        continue
    elseif QC == 2
        qcPassInd = find(obs(i).hsQC == 3 | obs(i).hsQC == 2);
    else %default is strict QC
        qcPassInd = find(obs(i).hsQC == 3);
    end
    obs(i).time = obs(i).time(qcPassInd);
    obs(i).lat = obs(i).lat(qcPassInd);
    obs(i).lon = obs(i).lon(qcPassInd);
    obs(i).hs = obs(i).hs(qcPassInd);
%    obs(i).hsEr = obs(i).hsEr(qcPassInd);
    obs(i).hsQC = obs(i).hsQC(qcPassInd);
%    obs(i).satID = obs(i).satID(qcPassInd);
    % WIND
    % obs(i).wind = obs(i).wind(qcPassInd );
end

