function obs = getRY19AltimeterObsBuoy(loadSatList, buoyLat, buoyLon, altPath, QC)

count2 = 0;
[glyph, ~] = giveGlyph;

for i = 1:length(loadSatList)
    if isempty(altPath)
        altPath = [baseDir 'Datasets' glyph 'Satellite' glyph 'Altimeter' glyph 'Ribal_Young_2019' glyph];
    end
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
        % search for data that is within 1 degree of the min and max buoy coordinate
        if altNorth <= max(buoyLat) + 1 & altNorth >= min(buoyLat) - 1 & ...
                altEast <= max(buoyLon) + 1 & altEast >= min(buoyLon) -1
            %         if altNorth >= 40 & altNorth <= 41 & ...
            %                 altEast >= 200 & altEast <= 201
            fileName = [altFilePath altFileList(j).name];
            %            altInfo = ncinfo(fileName);
%            currentAltInfo = ncinfo(fileName);
            count2 = count2 + 1;
            
            % load data
            satTimeOffset = datenum([1985 01 01 00 00 00]); %same for all
            %ping the file to make sure its working
            try 
                info   = ncread(fileName);
            catch
                disp(['error opening ' fileName])
                continue
            end
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
 %               obs(count2).hsKstd = ncread(fileName, 'SWH_KA_std_dev');
            else
                %                 obs(count2).hsK   = ncread(fileName, 'SWH_KU');
                obs(count2).hs   = ncread(fileName, 'SWH_KU_CAL');
                obs(count2).hsQC = ncread(fileName, 'SWH_KU_quality_control');
                %                 obs(count2).hsKno = ncread(fileName, 'SWH_KU_num_obs');
%                obs(count2).hsKstd = ncread(fileName, 'SWH_KU_std_dev');
            end            
            %WIND
            %             obs(count2).wind    = ncread(fileName, 'WSPD');
                obs(count2).wind = ncread(fileName, 'WSPD_CAL');
%                obs(count2).dist2coast = ncread(fileName, 'DIST2COAST');      
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
    switch QC
        case 2
            qcPassInd = find(obs(i).hsQC == 1 | obs(i).hsQC == 2);
        otherwise
            qcPassInd = find(obs(i).hsQC == 1);
    end
    obs(i).time = obs(i).time(qcPassInd);
    obs(i).lat = obs(i).lat(qcPassInd);
    obs(i).lon = obs(i).lon(qcPassInd);
    obs(i).hs = obs(i).hs(qcPassInd);
    obs(i).hsQC = obs(i).hsQC(qcPassInd);
    %WIND
    obs(i).wind = obs(i).wind(qcPassInd );
end