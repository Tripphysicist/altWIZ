function storm = getStormTrackIBTrACS(stormYear,stormName)
% function storm = getStormTrackIBTrACS(stormYear,stormName)
% Grab storm data from IBTRaCS netCDF file
% NOTES
% there are 150 vairables, most of them are rudundant information issued by
% 1 of the dozen or so agencies which tracks hurricanes
% The files are split by basin, and there is an global file.
%
%% TESTING
%stormYear = 1985:2020;
%stormName = 'all';
%stormYear = 2012;
%stormName = 'SANDY';
%%
fileName = '/Users/tripp/D/Analysis/altimeterComparison/storm-tracks/IBTrACS.NA.v04r00.nc';
time = ncread(fileName,'time') + datenum([1858 11 17 00 00 00]);
lon = ncread(fileName,'lon');
lat = ncread(fileName,'lat');
name =  ncread(fileName,'name');
numObs = ncread(fileName,'numobs');
year = ncread(fileName,'season');
vel = ncread(fileName,'storm_speed').*0.514444; %[m/s]
dir = ncread(fileName,'storm_dir');
maxWind = ncread(fileName, 'wmo_wind').*0.514444; %[m/s]
% minPres = ncread(fileName, 'wmo_pres'); %[mb]
rmw = ncread(fileName, 'usa_rmw').*1.852; %[km]
% status = ncread(fileName, 'nature');
count = 0;
[~, m] = size(time);
for i = 1:m
    currentYear = year(i);
    currentName = deblank(name(:,i)');
    if strcmp(stormName,'all')
        for j = 1:numel(stormYear)
            if currentYear == stormYear(j)
                disp(['found ' currentName ' ' num2str(stormYear(j))])
                count = count + 1;
                %                currentTime    = time(1:numObs(i),i);
                storm(count).time = time(1:numObs(i),i);
                storm(count).time    = time(1:numObs(i),i);
                storm(count).year    = currentYear;
                storm(count).lat     = lat(1:numObs(i),i);
                storm(count).lon     = lon(1:numObs(i),i);
                storm(count).maxWind = maxWind(1:numObs(i),i);
                % storm(count).minPres = minPres(1:numObs(i),i);
                %                 for ii = 1:numObs(i);
                %                     if ii ==1
                %                         storm(count).dir(ii) = NaN;
                %                         storm(count).vel(ii) = NaN;
                %                     else
                %                         [X,Y,D] = latlon2xy(storm(count).lat(ii),storm(count).lon(ii),...
                %                             storm(count).lat(ii-1),storm(count).lon(ii-1));
                %                         timeDiff = (storm(count).time(ii)-storm(count).time(ii-1))*24*60*60; %convert to seconds
                %                         storm(count).dir(ii) = deg360(90 - (180/pi).*atan2(Y,X));
                %                         storm(count).vel(ii) = 1000*D./timeDiff; %m/s
                %                     end
                %                 end
                %                storm(count).dir     = storm(count).dir';
                %                storm(count).vel     = storm(count).vel';
                storm(count).dir     = dir(1:numObs(i),i);
                storm(count).vel      = vel(1:numObs(i),i);
                storm(count).rmw     = rmw(1:numObs(i),i);
                storm(count).name    = currentName;
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
            %             for ii = 1:numObs(i);
            %                 if ii ==1
            %                     storm(count).dir(ii) = NaN;
            %                     storm(count).vel(ii) = NaN;
            %                 else
            %                     [X,Y,D] = latlon2xy(storm(count).lat(ii),storm(count).lon(ii),...
            %                         storm(count).lat(ii-1),storm(count).lon(ii-1));
            %                     timeDiff = (storm(count).time(ii)-storm(count).time(ii-1))*24*60*60; %convert to seconds
            %                     storm(count).dir(ii) = deg360(90 - (180/pi).*atan2(Y,X));
            %                     storm(count).vel(ii) = 1000*D./timeDiff; %m/s
            %                 end
            %             end
            %             storm(count).dir     = storm(count).dir';
            %             storm(count).vel     = storm(count).vel';
            storm(count).dir     = dir(1:numObs(i),i);
            storm(count).vel      = vel(1:numObs(i),i);
            storm(count).rmw     = rmw(1:numObs(i),i);
            storm(count).name    = currentName;
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