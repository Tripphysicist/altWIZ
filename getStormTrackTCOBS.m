function storm = getStormTrackTCOBS(stormYear,stormName)
% function storm = getStormTrackTCOBS(stormYear,stormName)
% Grab storm data from TCOBS v0.40 netCDF file
%
% NOTES
% I am having a lot of netCDF issues, it seems. vel, dir load but are
% nonsense, pressure doesn't load, name only loads with h5read
% there are 87 vairables, of these we are only interested in a few. See
% their documentation for more information
%
%% TESTING
%stormYear = 1985:2020;
%stormName = 'all';
%%
fileName = '/Users/tripp/D/Analysis/altimeterComparison/storm-tracks/TC-OBS_Database_0.40_20160119.nc';
% fileName = '/Users/tripp/D/Analysis/altimeterComparison/storm-tracks/TC-OBS_Database_0.42_20180209.nc';
time = (ncread(fileName,'TCOBS_time')./(60*60*24)) + datenum([1970 01 01 00 00 00]);
lon = ncread(fileName,'TCOBS_wind_center_position_longitude');
lat = ncread(fileName,'TCOBS_wind_center_position_latitude');
name = h5read(fileName,'/CYCLONE_ATCF_stormname');
listOfBadNames = {'INVEST','ONE','TWO','THREE','FOUR','FIVE',...
    'SIX','SEVEN','EIGHT','NINE','TEN','ELEVEN','TWELVE','THIRTEEN'...
    'FOURTEEN','FIFTEEN','SIXTEEN','SEVENTEEN','EIGHTTEEN','NINETEEN'...
    'TWENTY','TWENTY-ONE','TWENTY-TWO','UNNAMED','SUBTROP','NO_NAME'};

numObs = ncread(fileName, 'CYCLONE_TCOBS_npoints');
% year = ncread(fileName,'CYCLONE_ATCF_stormyear	');
% vel = ncread(fileName,'TCOBS_cyclone_translation_speed').*0.514444;
% dir = ncread(fileName,'TCOBS_cyclone_translation_direction');
maxWind = ncread(fileName, 'TCOBS_maximum_sustained_surface_wind').*0.514444; %m/s
% minPres = ncread(fileName, 'TCOBS_central_pressure'); %[milibar] something
% is wrong
rmw = ncread(fileName, 'TCOBS_maximum_sustained_surface_wind_radius').*1.852; % km

count = 0;
[~, m] = size(time);
for i = 1:m
    date = datevec(time(1:numObs,i));
    currentYear = date(1,1);
    currentName = cell2mat(name(i));
    if strcmp(stormName,'all')
        for j = 1:numel(stormYear)
            if currentYear == stormYear(j)
                disp(['found ' currentName ' ' num2str(stormYear(j))])
                count = count + 1;
                %                currentTime    = time(1:numObs(i),i);
                storm(count).time    = time(1:numObs(i),i);
                storm(count).year    = currentYear;
                storm(count).lat     = lat(1:numObs(i),i);
                storm(count).lon     = lon(1:numObs(i),i);
                storm(count).maxWind = maxWind(1:numObs(i),i);
                % storm(count).minPres = minPres(1:numObs(i),i);
                for ii = 1:numObs(i)
                    if ii ==1
                        storm(count).dir(ii) = NaN;
                        storm(count).vel(ii) = NaN;
                    else
                        [D,dirTemp] = latlon2dist(storm(count).lat(ii),storm(count).lon(ii),...
                            storm(count).lat(ii-1),storm(count).lon(ii-1));
                        timeDiff = (storm(count).time(ii)-storm(count).time(ii-1))*24*60*60; %convert to seconds
                        dirTemp(dirTemp<0) = dirTemp(dirTemp<0)+360;
                        storm(count).dir(ii) = dirTemp;
                        storm(count).vel(ii) = 1000*D./timeDiff; %m/s
                    end
                end
                storm(count).dir     = storm(count).dir';
                storm(count).vel     = storm(count).vel';
                storm(count).rmw     = rmw(1:numObs(i),i);
                for ll =1:length(listOfBadNames)
                    if strcmp(currentName,listOfBadNames{ll})
                        storm(count).name = 'NOT_NAMED';
                        break
                    else
                        storm(count).name    = currentName;
                        
                    end
                end
            end
        end
    else
        if currentYear == stormYear & strcmp(currentName,stormName)
            disp(['found ' stormName ' ' num2str(stormYear)])
            count = count + 1;
            %                currentTime    = time(1:numObs(i),i);
            storm(count).time = time(1:numObs(i),i);
            storm(count).time    = time(1:numObs(i),i);
            storm(count).year    = currentYear;
            storm(count).lat     = lat(1:numObs(i),i);
            storm(count).lon     = lon(1:numObs(i),i);
            storm(count).maxWind = maxWind(1:numObs(i),i);
            % storm(count).minPres = minPres(1:numObs(i),i);
            for ii = 1:numObs(i)
                if ii ==1
                    storm(count).dir(ii) = NaN;
                    storm(count).vel(ii) = NaN;
                else
                    [D,dirTemp] = latlon2dist(storm(count).lat(ii),storm(count).lon(ii),...
                        storm(count).lat(ii-1),storm(count).lon(ii-1));
                    timeDiff = (storm(count).time(ii)-storm(count).time(ii-1))*24*60*60; %convert to seconds
                    dirTemp(dirTemp<0) = dirTemp(dirTemp<0)+360;
                    storm(count).dir(ii) = dirTemp;
                    storm(count).vel(ii) = 1000*D./timeDiff; %m/s
                end
            end
            storm(count).dir     = storm(count).dir';
            storm(count).vel     = storm(count).vel';
            storm(count).rmw     = rmw(1:numObs(i),i);
            for ll =1:length(listOfBadNames)
                if strcmp(currentName,listOfBadNames{ll})
                    storm(count).name = 'NOT_NAMED';
                    break
                else
                    storm(count).name    = currentName;
                    
                end
            end
        end
    end
end
% dealing with NaNs
fields = fieldnames(storm);
for i = 1:length(storm)
    for idx = 1:length(fields)
        if strcmp((fields{idx}),'lon') | strcmp((fields{idx}),'name')
            continue
        end
        nanIndex = storm(i).((fields{idx}))<0;
        storm(i).(fields{idx})(nanIndex) = NaN;
    end
end