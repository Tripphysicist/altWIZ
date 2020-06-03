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
% tic
% debugCount = 0;
tic
R = 6371; % Earth's radius in km
lat1=obsLat*pi/180;
lon1=obsLon*pi/180;

for i = 1:buoyTimeLength
    %    [xdummy,ydummy,distance]=latlon2xy(obsLat,obsLon,buoyLat(i),buoyLon(i));
    lat2=buoyLat(i)*pi/180;
    lon2=buoyLon(i)*pi/180;
    deltaLat=lat2-lat1;
    deltaLon=lon2-lon1;
    a=sin((deltaLat)/2).^2 + cos(lat1).*cos(lat2) .* sin(deltaLon/2).^2;
    c=2*atan2(sqrt(a),sqrt(1-a));
    distance = R.*c;    %Haversine distance    
    [meanIndex , ~] = find(distance<=maxDistance &  abs(buoyTime(i) - obsTime) <= maxTimeDiff);
    if length(meanIndex) >= minNumberObs
        currentHsData = obsHs(meanIndex);
        meanHsByDistance(i) = nanmean(currentHsData);
        stdHsByDistance(i) = nanstd(currentHsData);        
        currentWindData = obsWind(meanIndex);
        meanWindByDistance(i) = nanmean(currentWindData);
        stdWindByDistance(i) = nanstd(currentWindData);
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