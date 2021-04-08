function [LONobs, LATobs, TIMEobs, HSmean, HSobsSTD, HSnearest, WINDmean, WINDobsSTD, WINDnearest, indNaN] = bubbleMethod(obsLon, obsLat, obsTime, obsHs, obsWind, buoyLon, buoyLat, buoyTime,maxDistance,maxTimeDiff,minNumberObs)
% [LONobs LATobs TIMEobs HSobs] = bubbleMethod(obsLon, obsLat, obsTime, obsHs, gridLon,gridLat,gridTime,maxDistance,maxTimeDiff,gridStatus,minNumberObs)

% for debugging
% obsLon = LONobsNA;
% obsLat = LATobsNA;
% obsTime = TIMEobsNA;
% obsHs = HSobsNA;
% obsWind = WINDobsNA;
% buoyLon = buoyData.lon;
% buoyLat = buoyData.lat;
% buoyTime = buoyData.time;


buoyTimeLength = length(buoyTime);
meanHsByDistance = NaN(buoyTimeLength,1);
meanWindByDistance = NaN(buoyTimeLength,1);
stdHsByDistance = NaN(buoyTimeLength,1);
stdWindByDistance = NaN(buoyTimeLength,1);
nearestHs = NaN(buoyTimeLength,1);
nearestWind = NaN(buoyTimeLength,1);

% tic
% debugCount = 0;
tic

for i = 1:buoyTimeLength
    %    [xdummy,ydummy,distance]=latlon2xy(obsLat,obsLon,buoyLat(i),buoyLon(i));
    distance = latlon2dist(obsLat,obsLon,buoyLat(i),buoyLon(i));
    [meanIndex , ~] = find(distance<=maxDistance &  abs(buoyTime(i) - obsTime) <= maxTimeDiff);
    if length(meanIndex) >= minNumberObs
        currentHsData = obsHs(meanIndex);
        meanHsByDistance(i) = mean(currentHsData,'omitnan');
        stdHsByDistance(i) = std(currentHsData,'omitnan');        
        currentWindData = obsWind(meanIndex);
        meanWindByDistance(i) = mean(currentWindData,'omitnan');
        stdWindByDistance(i) = std(currentWindData,'omitnan');
        %         debugCount = debugCount +1
        currentDistance = distance(meanIndex);
        [~ , sortIndex] = sort(currentDistance);
        nearestHs(i) = currentHsData(sortIndex(1));
        nearestWind(i) = currentWindData(sortIndex(1));        
        disp([num2str(100*i./buoyTimeLength) ' %'])
    end
end
% timerM1 = toc
% stdHsByBin = stdHsByBin(:);

indNaN = ~isnan(meanHsByDistance);
HSmean = meanHsByDistance(indNaN);
HSnearest = nearestHs(indNaN);
WINDmean = meanWindByDistance(indNaN);
WINDnearest = nearestWind(indNaN);
HSobsSTD = stdHsByDistance(indNaN);
WINDobsSTD = stdWindByDistance(indNaN);
LONobs = buoyLon(indNaN);
LATobs = buoyLat(indNaN);
TIMEobs = buoyTime(indNaN);