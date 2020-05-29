function obs = getAltimeterObsRY19(loadSatList, mdLat, mdLon, mdLonO, altPath, crossesPrime)
% obs = getAltimeterObsWW3(loadSatList, mdLat, mdLon, mdLonO, altPath, crossesPrime)
% function that gets loads altimeter data into the workspace depending on
% the model grid provided.

[glyph baseDir] = giveGlyph;
count2 = 0;
for i = 1:length(loadSatList)
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
                
                %WIND
                %             obs(count2).wind    = ncread(fileName, 'WSPD');
                %             obs(count2).windCal = ncread(fileName, 'WSPD_CAL');
                
                
            end
        else
            if altNorth + 1 >= mdLat(1) & altNorth <= mdLat(end) & ...
                    altEast + 1 >= mdLon(1) & altEast  <= mdLon(end)
                
%                 disp(['model lon - lat ' num2str(mdLon(1)) '-' num2str(mdLon(end))...
%                     ':' num2str(mdLat(1)) '-' num2str(mdLat(end))])
%                 disp(['altNorth = ' num2str(altNorth)])
%                 disp(['normal altEast = ' num2str(altEast)])
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
%                    obs(count2).hsKstd = ncread(fileName, 'SWH_KA_std_dev');
                else
                    %                 obs(count2).hsK   = ncread(fileName, 'SWH_KU');
                    obs(count2).hsKcal   = ncread(fileName, 'SWH_KU_CAL');
                    obs(count2).hsKqc = ncread(fileName, 'SWH_KU_quality_control');
                    %                 obs(count2).hsKno = ncread(fileName, 'SWH_KU_num_obs');
%                    obs(count2).hsKstd = ncread(fileName, 'SWH_KU_std_dev');
                end
%               obs(count2).dist2coast = ncread(fileName, 'DIST2COAST');
                %WIND
                %             obs(count2).wind    = ncread(fileName, 'WSPD');
                %             obs(count2).windCal = ncread(fileName, 'WSPD_CAL');
                
                
            end
        end
    end
end




