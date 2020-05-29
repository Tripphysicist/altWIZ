function [LONobs LATobs TIMEobs HSobsMean HSobsStd HSobsNoSam] = bubbleMethod(obsLon, obsLat, obsTime, obsHs, gridLon,gridLat,gridTime,maxDistance,maxTimeDiff,gridStatus,minNumberObs) 
% [LONobs LATobs TIMEobs HSobs] = bubbleMethod(obsLon, obsLat, obsTime, obsHs, gridLon,gridLat,gridTime,maxDistance,maxTimeDiff,gridStatus,minNumberObs) 
%
%troubleshooting
% obsLon = LONobsNA; obsLat = LATobsNA; obsTime = TIMEobsNA; 
% obsHs = HSobsNA; gridLon = mdTest.lon; gridLat = mdTest.lat;
% gridTime = mdTest.time; gridStatus = mdTest.gridstatus;
% 
%
%
%TO DO: put in a nearest option so it grabs just one altimeter sample
%instead of an average

if sum([ismember(gridLon,359.5); ismember(gridLon,0)]) == 2
    for i = 1:length(gridLon)
        if gridLon(i) > 180
            gridLon(i) = gridLon(i) - 360;
        end       
    end
    for i = 1:length(obsLon)
        if obsLon(i) > 180
            obsLon(i) = obsLon(i) - 360;
        end       
    end
end

mdLonLength = length(gridLon);
mdLatLength = length(gridLat);

skipTime = ceil(maxTimeDiff/(gridTime(2)-gridTime(1)));
mdTimeLength = length(gridTime(1:skipTime:end));

meanHsByDistance = NaN(mdLonLength,mdLatLength,mdTimeLength);
stdHsByDistance = NaN(mdLonLength,mdLatLength,mdTimeLength);
noSamHsByDistance = NaN(mdLonLength,mdLatLength,mdTimeLength);
% tic
for i = 1:mdLonLength
    for j = 1:mdLatLength
        if gridStatus(i,j)==1
            for k = 1:skipTime:mdTimeLength
                [xdummy,ydummy,distance] = latlon2xy(obsLat,obsLon,gridLat(j),gridLon(i));
                [meanIndex dum2] = find(distance<=maxDistance &...
                    abs(gridTime(k)- obsTime) <= maxTimeDiff);                
                if length(meanIndex) >= minNumberObs
                    currentData = obsHs(meanIndex);
                    meanHsByDistance(i,j,k) = nanmean(currentData);
                    stdHsByDistance(i,j,k) = nanstd(currentData);
                    noSamHsByDistance(i,j,k) = length(currentData);                                        
                end
            end
        end
    end
end
% timerM1 = toc
meanHsByDistance = meanHsByDistance(:);
stdHsByDistance  = stdHsByDistance(:);
noSamHsByDistance = noSamHsByDistance(:);

indNaN = ~isnan(meanHsByDistance);

HSobsMean  = meanHsByDistance(indNaN);
HSobsStd   = stdHsByDistance(indNaN);
HSobsNoSam = noSamHsByDistance(indNaN);

[LON, LAT, TIME] = ndgrid(gridLon,gridLat,gridTime(1:skipTime:end));
LON  = LON(:);
LAT  = LAT(:);
TIME = TIME(:);

LONobs = LON(indNaN);
LATobs = LAT(indNaN);
TIMEobs = TIME(indNaN);

if sum(LONobs<0) > 0
    for i = 1:length(LONobs)
        if LONobs(i) < 0
            LONobs(i) = LONobs(i) + 360;
        end
    end
end