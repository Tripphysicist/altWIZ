function obs = getAltimeterObsBuoy(loadSatList, buoyLat, buoyLon, altPath)

count2 = 0;
[glyph baseDir] = giveGlyph;

for i = 1:length(loadSatList)
    if isempty(altPath)
        altPath = [baseDir 'Datasets' glyph 'Satellite' glyph 'Altimeter' glyph 'Ribal_Young_2019' glyph];
    end
    altFilePath = [altPath loadSatList(i).name glyph];
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
        if altNorth <= max(buoyLat) + 1 & altNorth >= min(buoyLat) - 1 & ...
                altEast <= max(buoyLon) + 1 & altEast >= min(buoyLon) -1
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
 %               obs(count2).hsKstd = ncread(fileName, 'SWH_KA_std_dev');
            else
                %                 obs(count2).hsK   = ncread(fileName, 'SWH_KU');
                obs(count2).hsKcal   = ncread(fileName, 'SWH_KU_CAL');
                obs(count2).hsKqc = ncread(fileName, 'SWH_KU_quality_control');
                %                 obs(count2).hsKno = ncread(fileName, 'SWH_KU_num_obs');
%                obs(count2).hsKstd = ncread(fileName, 'SWH_KU_std_dev');
            end
            
            %WIND
            %             obs(count2).wind    = ncread(fileName, 'WSPD');
                obs(count2).windCal = ncread(fileName, 'WSPD_CAL');
                obs(count2).dist2coast = ncread(fileName, 'DIST2COAST');
            
            
        end
    end
end
