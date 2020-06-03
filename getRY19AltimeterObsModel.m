function obs = getRY19AltimeterObsModel(loadSatList, mdLat, mdLon, mdLonO, altPath, crossesPrime, QC)
% obs = getRY19AltimeterObsModel(loadSatList, mdLat, mdLon, mdLonO, altPath, crossesPrime)
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% % Script to load RY19 data in to workspace                              %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%
%  OUTPUT
%  obs - structure of wave height data from altimeter and models with
%  shared lon, lat, time
%  obs.lon - longitude in degrees
%  obs.lat - latitude in degrees
%  obs.time - time (MATLAB datenum)
%  obs.hs - calibrated as detailed in RY19
%  obs.hsEr - the standard deviation of the hs raw data
%  obs.hsQC - flags: In the present database, a series of data flags defined as 1, 2,
%             3, 4, and 9 represent Good_data, Probably_good_data, SAR-mode data or
%             possible hardware error (only used for CRYOSAT-2), Bad_data and
%             Missing_data, respectively, have been used.
%
% QC NOTE: observations <50km from a coast are flagged as probably good.
% So to access these, one would need to also allow probably
% good data
%
%
%  INPUT
%  loadSatList - list of satellites to load
%  mdLon  - model longitude (0 - 360)
%  mdLonO - model longitude (-180 - 180)
%  mdLat  - model latitude
%  altPath - path to RY19 data
%  crossesPrime - logical to indicate whether the longitude grid crosses
%                 the prime merdian
%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% % Author:       Clarence Olin Collins III, Ph.D.                        %
% % Affiliation:  ERDC - FRF                                              %
% % created:      06/01/2020                                              %
% % version:      1.0                                                     %
% % updates:                                                              %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

[glyph , ~] = giveGlyph;
count2 = 0;
for i = 1:length(loadSatList)
    altFilePath = [altPath loadSatList(i).name glyph];
    altFileList = dir(altFilePath);
    altFileList = rmfield (altFileList,{'date','bytes','isdir','datenum'});
    altFileList = altFileList(3:end);
    for j = 1:length(altFileList)
        altNorth  = str2double(altFileList(j).name(end-16:end-14));
        if altFileList(j).name(end-13) == 'S'
            altNorth = -altNorth;
        else
        end
        
        altEast   = str2double(altFileList(j).name(end-11:end-9));
        
        % for rectilinear grid, add one to the greater than test because
        % each alt file runs from altNorth to altNorth + 1
        if crossesPrime
            if altEast > 180
                altEast = altEast - 360;
            end
            if altNorth + 1 >= mdLat(1) & altNorth <= mdLat(end) & ...
                    altEast + 1 >= mdLonO(1) & altEast   <= mdLonO(end)
                
                %                 disp(['model lon - lat ' num2str(mdLon(1)) '-' num2str(mdLon(end))...
                %                     ':' num2str(mdLat(1)) '-' num2str(mdLat(end))])
                %                 disp(['altNorth = ' num2str(altNorth)])
                %                 disp(['prime normal altEast = ' num2str(altEast)])
                
                fileName = [altFilePath altFileList(j).name];
                %            altInfo = ncinfo(fileName);
                count2 = count2 + 1;
                
                % load data
                
                %ping the file to make sure its working
                try
                    info   = ncread(fileName);
                catch
                    disp(['error opening ' fileName])
                    continue
                end
                
                
                satTimeOffset = datenum([1985 01 01 00 00 00]); %same for all
                obs(count2).time = ncread(fileName, 'TIME') + satTimeOffset;
                obs(count2).lat  = ncread(fileName, 'LATITUDE');
                obs(count2).lat  = double(obs(count2).lat);
                obs(count2).lon  = ncread(fileName, 'LONGITUDE');
                obs(count2).lon  = double(obs(count2).lon);
                obs(count2).satID = loadSatList(i).name;
                % try C band
                %             obs(count2).hsC   = ncread(fileName, 'SWH_C');
                %             obs(count2).hsCqc = ncread(fileName, 'SWH_C_quality_control');
                %             obs(count2).hsCno = ncread(fileName, 'SWH_C_num_obs');
                %             obs(count2).hsCstd = ncread(fileName, 'SWH_C_std_dev');
                
                if  strcmp(loadSatList(i).name,'SARAL')
                    %                 obs(count2).hsK   = ncread(fileName, 'SWH_KA');
                    obs(count2).hs   = ncread(fileName, 'SWH_KA_CAL');
                    %                    obs(count2).hsEr   = ncread(fileName, 'SWH_KA_std_dev');
                    obs(count2).hsQC = ncread(fileName, 'SWH_KA_quality_control');
                    %                 obs(count2).hsKno = ncread(fileName, 'SWH_KA_num_obs');
                else
                    %                 obs(count2).hsK   = ncread(fileName, 'SWH_KU');
                    obs(count2).hs   = ncread(fileName, 'SWH_KU_CAL');
                    %                    obs(count2).hsEr   = ncread(fileName, 'SWH_KU_std_dev');
                    obs(count2).hsQC = ncread(fileName, 'SWH_KU_quality_control');
                    %                 obs(count2).hsKno = ncread(fileName, 'SWH_KU_num_obs');
                end
                
                %WIND
                %             obs(count2).wind    = ncread(fileName, 'WSPD');
                %             obs(count2).windCal = ncread(fileName, 'WSPD_CAL');
                
                
            end
        else
            if altNorth + 1 >= mdLat(1) & altNorth <= mdLat(end) & ...
                    altEast + 1 >= mdLon(1) & altEast  <= mdLon(end)
                
                fileName = [altFilePath altFileList(j).name];

                count2 = count2 + 1;
                
                
                %ping the file to make sure its working
                try
                    info   = ncread(fileName);
                catch
                    disp(['error opening ' fileName])
                    continue
                end
                
                % load data
                satTimeOffset = datenum([1985 01 01 00 00 00]); %same for all
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
                    obs(count2).hs   = ncread(fileName, 'SWH_KA_CAL');
                    obs(count2).hsQC = ncread(fileName, 'SWH_KA_quality_control');
                    %                 obs(count2).hsKno = ncread(fileName, 'SWH_KA_num_obs');
                    %                    obs(count2).hsEr = ncread(fileName, 'SWH_KA_std_dev');
                else
                    %                 obs(count2).hsK   = ncread(fileName, 'SWH_KU');
                    obs(count2).hs   = ncread(fileName, 'SWH_KU_CAL');
                    obs(count2).hsQC = ncread(fileName, 'SWH_KU_quality_control');
                    %                 obs(count2).hsKno = ncread(fileName, 'SWH_KU_num_obs');
                    %                 obs(count2).hsEr = ncread(fileName, 'SWH_KU_std_dev');
                end
                %WIND
                %             obs(count2).wind = ncread(fileName, 'WSPD_CAL');
                
                
            end
        end
    end
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
    if QC == 5
        continue
    elseif QC == 2
        qcPassInd = find(obs(i).hsQC == 1 | obs(i).hsQC == 2);
    else %default is strict QC
        qcPassInd = find(obs(i).hsQC == 1);
    end
    obs(i).time = obs(i).time(qcPassInd);
    obs(i).lat = obs(i).lat(qcPassInd);
    obs(i).lon = obs(i).lon(qcPassInd);
    obs(i).hs = obs(i).hs(qcPassInd);
    obs(i).hsQC = obs(i).hsQC(qcPassInd);
    %    obs(i).hsEr = obs(i).hsEr(qcPassInd);
    %    obs(i).satID = obs(i).satID(qcPassInd);
    
    %WIND
    % obs(i).wind = obs(i).wind(qcPassInd );
    
end