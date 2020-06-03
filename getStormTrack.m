function storm = getStormTrack(stormYear,stormName,subsample)
% function storm = getStormTrack(stormYear,stormName)
% Grab storm data from HURDAT2 file
% NOTES
% (Spaces 1-4) ? Year
% (Spaces 5-6) ? Month
% (Spaces 7-8, before 1st comma) ? Day
% (Spaces 11-12) ? Hours in UTC (Universal Time Coordinate)
% (Spaces 13-14, before 2nd comma) ? Minutes
% (Space 17, before 3rd comma) ? Record identifier (see notes below)
% (Spaces 20-21, before 4th comma) ? Status of system. Options are:
% (Spaces 24-27) ? Latitude
% (Space 28, before 5th comma) ? Hemisphere ? North or South
% (Spaces 31-35) ? Longitude
% (Space 36, before 6th comma) ? Hemisphere ? West or East
% (Spaces 39-41, before 7th comma) ? Maximum sustained wind (in knots)
% (Spaces 44-47, before 8th comma) ? Minimum Pressure (in millibars)
% (Spaces 50-53, before 9th comma) ? 34 kt wind radii maximum extent in northeastern quadrant (in nautical miles)
% (Spaces 56-59, before 10th comma) ? 34 kt wind radii maximum extent in southeastern quadrant (in nautical miles)
% (Spaces 62-65, before 11th comma) ? 34 kt wind radii maximum extent in southwestern quadrant (in nautical miles)
% (Spaces 68-71, before 12th comma) ? 34 kt wind radii maximum extent in northwestern quadrant (in nautical miles)
% (Spaces 74-77, before 13th comma) ? 50 kt wind radii maximum extent in northeastern quadrant (in nautical miles)
% (Spaces 80-83, before 14th comma) ? 50 kt wind radii maximum extent in southeastern quadrant (in nautical miles)
% (Spaces 86-89, before 15th comma) ? 50 kt wind radii maximum extent in southwestern quadrant (in nautical miles)
% (Spaces 92-95, before 16th comma) ? 50 kt wind radii maximum extent in northwestern quadrant (in nautical miles)
% (Spaces 98-101, before 17th comma) ? 64 kt wind radii maximum extent in northeastern quadrant (in nautical miles)
% (Spaces 104-107, before 18th comma) ? 64 kt wind radii maximum extent in southeastern quadrant (in nautical miles)
% (Spaces 110-113, before 19th comma) ? 64 kt wind radii maximum extent in southwestern quadrant (in nautical miles)
% (Spaces 116-119, before 20th comma) ? 64 kt wind radii maximum extent in northwestern quadrant (in nautical miles)
fid = fopen('/Users/tripp/D/Analysis/altimeterComparison/storm-tracks/hurdat2-1851-2019-052520.txt');
line=1;
if strcmp(stormName,'all')
    count = 0;
    
    while line ~= -1
        line = fgetl(fid);
        if strcmp(line(1),'A')
            commaInd = strfind(line,',');
            nameWithSpaces = line(commaInd(1)+1:commaInd(2)-1);
            spaceInd = strfind(nameWithSpaces,' ');
            name = nameWithSpaces(spaceInd(end)+1:end);
            nameYear = str2double(line(commaInd(1)-4:commaInd(1)-1));
            disp(['found ' name ' ' num2str(nameYear)])
            line = fgetl(fid);
            stormCount = 0;
            while line(1) ~= 'A' & line ~= -1
                count = count + 1;
                stormCount = stormCount + 1;
                %storm.name = stormName;
                
                year = str2double(line(1:4));
                month = str2double(line(5:6));
                day = str2double(line(7:8));
                hour = str2double(line(11:12));
                minute = str2double(line(13:14));
                
                storm.time(count)    = datenum([year month day hour minute 0]); %time
                
                storm.lat(count)     = str2double(line(24:27));
                % storm.latHemi        = line(28);
                storm.lon(count)     = -str2double(line(31:35));
                % storm.lonHemi        = line(36);
                storm.maxWind(count) = str2double(line(39:41));
                storm.minPres(count) = str2double(line(44:47));
                storm.r34NE(count)   = str2double(line(50:53));
                storm.r34SE(count)   = str2double(line(56:59));
                storm.r34SW(count)   = str2double(line(62:65));
                storm.r34NW(count)   = str2double(line(68:71));
                storm.r50NE(count)   = str2double(line(74:77));
                storm.r50SE(count)   = str2double(line(80:83));
                storm.r50SW(count)   = str2double(line(86:89));
                storm.r50NW(count)   = str2double(line(92:95));
                storm.r64NE(count)   = str2double(line(98:101));
                storm.r64SE(count)   = str2double(line(104:107));
                storm.r64SW(count)   = str2double(line(110:113));
                storm.r64NW(count)   = str2double(line(116:119));
                
                storm.index(count) = stormCount;
                
                
                if stormCount == 1
                    storm.dir(count) = NaN;
                    storm.vel(count) = NaN;
                else
                    
                    [X,Y,D] = latlon2xy(storm.lat(count),storm.lon(count),storm.lat(count-1),storm.lon(count-1));
                    timeDiff = (storm.time(count)-storm.time(count-1))*24*60*60; %convert to seconds
                    storm.dir(count) = deg360(90 - (180/pi).*atan2(Y,X));
                    storm.vel(count) = 1000*D./timeDiff; %m/s
                end
                % convert metadata info to numeric code
                % 1  C ? Closest approach to a coast, not followed by a landfall G ? Genesis
                % 2  I ? An intensity peak in terms of both pressure and wind
                % 3  L ? Landfall (center of system crossing a coastline)
                % 4  P ? Minimum in central pressure
                % 5  R ? Provides additional detail on the intensity of the cyclone when rapid changes are underway
                % 6  S ? Change of status of the system
                % 7  T ? Provides additional detail on the track (position) of the cyc
                meta    = line(17);
                switch meta
                    
                    case 'C'
                        storm.meta(count) = 1;
                    case 'I'
                        storm.meta(count) = 2;
                    case 'L'
                        storm.meta(count) = 3;
                    case 'P'
                        storm.meta(count) = 4;
                    case 'R'
                        storm.meta(count) = 5;
                    case 'S'
                        storm.meta(count) = 6;
                    case 'T'
                        storm.meta(count) = 7;
                    otherwise
                        storm.meta(count) = NaN;
                end
                % convert status info to numeric code
                % 1  TD ? Tropical cyclone of tropical depression intensity (< 34 knots)
                % 2  TS ? Tropical cyclone of tropical storm intensity (34-63 knots)
                % 3  HU ? Tropical cyclone of hurricane intensity (> 64 knots)
                % 4  EX ? Extratropical cyclone (of any intensity)
                % 5  SD ? Subtropical cyclone of subtropical depression intensity (< 34 knots)
                % 6  SS ? Subtropical cyclone of subtropical storm intensity (> 34 knots)
                % 7  LO ? A low that is neither a tropical cyclone, a subtropical cyclone, nor an extratropical cyclone (of any intensity)
                % 8  WV ? Tropical Wave (of any intensity)
                % 9  DB ? Disturbance (of any intensity)
                
                status    = line(20:21);
                switch status
                    
                    case 'TD'
                        storm.status(count) = 1;
                    case 'TS'
                        storm.status(count) = 2;
                    case 'HU'
                        storm.status(count) = 3;
                    case 'EX'
                        storm.status(count) = 4;
                    case 'SD'
                        storm.status(count) = 5;
                    case 'SS'
                        storm.status(count) = 6;
                    case 'LO'
                        storm.status(count) = 7;
                    case 'WV'
                        storm.status(count) = 8;
                    case 'DB'
                        storm.status(count) = 9;
                    otherwise
                        storm.status(count) = NaN;
                end
                %subsample
                
                
                
                line = fgetl(fid);
            end
        end
    end
else
    while line ~= -1
        line = fgetl(fid);
        if strcmp(line(1),'A')
            commaInd = strfind(line,',');
            nameWithSpaces = line(commaInd(1)+1:commaInd(2)-1);
            spaceInd = strfind(nameWithSpaces,' ');
            name = nameWithSpaces(spaceInd(end)+1:end);
            nameYear = str2double(line(commaInd(1)-4:commaInd(1)-1));
            if strcmp(name,stormName) && nameYear == stormYear
                disp('found it!')
                count = 0;
                line = fgetl(fid);
                while line(1) ~= 'A'
                    count = count + 1;
                    % storm.name = stormName;
                    
                    year = str2double(line(1:4));
                    month = str2double(line(5:6));
                    day = str2double(line(7:8));
                    hour = str2double(line(11:12));
                    minute = str2double(line(13:14));
                    
                    storm.time(count)    = datenum([year month day hour minute 0]); %time
                    storm.lat(count)     = str2double(line(24:27));
                    % storm.latHemi        = line(28);
                    storm.lon(count)     = -str2double(line(31:35));
                    % storm.lonHemi        = line(36);
                    
                    storm.maxWind(count) = str2double(line(39:41));
                    storm.minPres(count) = str2double(line(44:47));
                    storm.r34NE(count)   = str2double(line(50:53));
                    storm.r34SE(count)   = str2double(line(56:59));
                    storm.r34SW(count)   = str2double(line(62:65));
                    storm.r34NW(count)   = str2double(line(68:71));
                    storm.r50NE(count)   = str2double(line(74:77));
                    storm.r50SE(count)   = str2double(line(80:83));
                    storm.r50SW(count)   = str2double(line(86:89));
                    storm.r50NW(count)   = str2double(line(92:95));
                    storm.r64NE(count)   = str2double(line(98:101));
                    storm.r64SE(count)   = str2double(line(104:107));
                    storm.r64SW(count)   = str2double(line(110:113));
                    storm.r64NW(count)   = str2double(line(116:119));
                    
                    if count == 1
                        storm.dir(count) = NaN;
                        storm.vel(count) = NaN;
                    else
                        [X,Y,D] = latlon2xy(storm.lat(count),storm.lon(count),storm.lat(count-1),storm.lon(count-1));
                        timeDiff = (storm.time(count)-storm.time(count-1))*24*60*60; %convert to seconds
                        storm.dir(count) = deg360(90 - (180/pi).*atan2(Y,X));
                        storm.vel(count) = 1000*D./timeDiff; %m/s
                    end
                    
                    % convert metadata info to numeric code
                    % 1  C ? Closest approach to a coast, not followed by a landfall G ? Genesis
                    % 2  I ? An intensity peak in terms of both pressure and wind
                    % 3  L ? Landfall (center of system crossing a coastline)
                    % 4  P ? Minimum in central pressure
                    % 5  R ? Provides additional detail on the intensity of the cyclone when rapid changes are underway
                    % 6  S ? Change of status of the system
                    % 7  T ? Provides additional detail on the track (position) of the cyc
                    meta    = line(17);
                    switch meta     
                        case 'C'
                            storm.meta(count) = 1;
                        case 'I'
                            storm.meta(count) = 2;
                        case 'L'
                            storm.meta(count) = 3;
                        case 'P'
                            storm.meta(count) = 4;
                        case 'R'
                            storm.meta(count) = 5;
                        case 'S'
                            storm.meta(count) = 6;
                        case 'T'
                            storm.meta(count) = 7;
                        otherwise
                            storm.meta(count) = NaN;
                    end
                    % convert status info to numeric code
                    % 1  TD ? Tropical cyclone of tropical depression intensity (< 34 knots)
                    % 2  TS ? Tropical cyclone of tropical storm intensity (34-63 knots)
                    % 3  HU ? Tropical cyclone of hurricane intensity (> 64 knots)
                    % 4  EX ? Extratropical cyclone (of any intensity)
                    % 5  SD ? Subtropical cyclone of subtropical depression intensity (< 34 knots)
                    % 6  SS ? Subtropical cyclone of subtropical storm intensity (> 34 knots)
                    % 7  LO ? A low that is neither a tropical cyclone, a subtropical cyclone, nor an extratropical cyclone (of any intensity)
                    % 8  WV ? Tropical Wave (of any intensity)
                    % 9  DB ? Disturbance (of any intensity)
                    
                    status    = line(20:21);
                    switch status
                        case 'TD'
                            storm.status(count) = 1;
                        case 'TS'
                            storm.status(count) = 2;
                        case 'HU'
                            storm.status(count) = 3;
                        case 'EX'
                            storm.status(count) = 4;
                        case 'SD'
                            storm.status(count) = 5;
                        case 'SS'
                            storm.status(count) = 6;
                        case 'LO'
                            storm.status(count) = 7;
                        case 'WV'
                            storm.status(count) = 8;
                        case 'DB'
                            storm.status(count) = 9;
                        otherwise
                            storm.status(count) = NaN;
                    end
                    storm.index(count) = count;
                    
                    line = fgetl(fid);
                end
            end
        end
    end
end

fclose(fid);


fields = fieldnames(storm);
for idx = 1:length(fields)
    nanIndex = storm.((fields{idx}))==-999;
    storm.(fields{idx})(nanIndex) = NaN;
end
