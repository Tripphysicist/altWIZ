function storm = getStormTrack(stormYear,stormName)
% function storm = getStormTrack(stormYear,stormName)
% Grab storm data from HURDAT2 file

fid = fopen('/Users/tripp/D/Analysis/altimeterComparison/storm-tracks/hurdat2-1851-2019-052520.txt');
line=1;
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
                commaInd = strfind(line,',');
                year = str2double(line(1:4));
                month = str2double(line(5:6));
                day = str2double(line(7:8));
                hour = str2double(line(commaInd(1)+2:commaInd(1)+3));
                minute = str2double(line(commaInd(1)+4:commaInd(1)+5));
                storm.time(count) = datenum([year month day hour minute 0]);
                latSpace = line(commaInd(4)+1:commaInd(5)-1);
                spaceInd = strfind(latSpace,' ');
                storm.lat(count) = str2double(latSpace(spaceInd(end)+1:end-1));
                lonSpace = line(commaInd(5)+1:commaInd(6)-1);
                spaceInd = strfind(lonSpace,' ');
                storm.lon(count) = -str2double(lonSpace(spaceInd(end)+1:end-1));
                line = fgetl(fid);
            end
        end
    end
end
fclose(fid);