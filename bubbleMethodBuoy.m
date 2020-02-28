function [LONobs LATobs TIMEobs HSobs, WINDobs, indNaN] = bubbleMethod(obsLon, obsLat, obsTime, obsHs, obsWind, buoyLon, buoyLat, buoyTime,maxDistance,maxTimeDiff,minNumberObs) 
% [LONobs LATobs TIMEobs HSobs] = bubbleMethod(obsLon, obsLat, obsTime, obsHs, gridLon,gridLat,gridTime,maxDistance,maxTimeDiff,gridStatus,minNumberObs) 

% for debugging
% obsLon = LONobsNA;
% obsLat = LATobsNA;
% obsTime = TIMEobsNA;
% obsHs = HSobsNA;
% buoyLon = buoyTest.lon;
% buoyLat = buoyTest.lat;
% buoyTime = buoyTest.time;


mdLonLength = length(buoyLon);
mdLatLength = length(buoyLat);

skipTime = ceil(maxTimeDiff/(buoyTime(2)-buoyTime(1)));
buoyTimeLength = length(buoyTime(1:skipTime:end));

meanHsByDistance = NaN(buoyTimeLength,1);
meanWindByDistance = NaN(buoyTimeLength,1);
% tic
% debugCount = 0;
for i = 1:skipTime:buoyTimeLength
    [xdummy,ydummy,distance]=latlon2xy(obsLat,obsLon,buoyLat(i),buoyLon(i));
    [meanIndex dum2] = find(distance<=maxDistance &  abs(buoyTime(i) - obsTime) <= maxTimeDiff);
    if length(meanIndex) >= minNumberObs
        currentHsData = obsHs(meanIndex);
        meanHsByDistance(i) = nanmean(currentHsData);
        currentWindData = obsWind(meanIndex);
        meanWindByDistance(i) = nanmean(currentWindData);
%         debugCount = debugCount +1        
    end
end
% timerM1 = toc
% stdHsByBin = stdHsByBin(:);

indNaN = ~isnan(meanHsByDistance);
HSobs = meanHsByDistance(indNaN);
WINDobs = meanWindByDistance(indNaN); 

LONobs = buoyLon(indNaN);
LATobs = buoyLat(indNaN);
TIMEobs = buoyTime(indNaN);