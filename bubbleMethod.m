openfunction [LONobs LATobs TIMEobs HSobs] = bubbleMethod(obsLon, obsLat, obsTime, obsHs, gridLon,gridLat,gridTime,maxDistance,maxTimeDiff,gridStatus,minNumberObs) 
% [LONobs LATobs TIMEobs HSobs] = bubbleMethod(obsLon, obsLat, obsTime, obsHs, gridLon,gridLat,gridTime,maxDistance,maxTimeDiff,gridStatus,minNumberObs) 

mdLonLength = length(gridLon);
mdLatLength = length(gridLat);

skipTime = ceil(maxTimeDiff/(gridTime(2)-gridTime(1)));
mdTimeLength = length(gridTime(1:skipTime:end));

meanHsByDistance = NaN(mdLonLength,mdLatLength,mdTimeLength);
% tic
for i = 1:mdLonLength
    for j = 1:mdLatLength
        if gridStatus(i,j)==1
            for k = 1:skipTime:mdTimeLength
                [xdummy,ydummy,distance]=latlon2xy(obsLat,obsLon,gridLat(j),gridLon(i));
                [meanIndex dum2] = find(distance<=maxDistance &  abs(gridTime(k)- obsTime) <= maxTimeDiff);
                if length(meanIndex) >= minNumberObs
                    currentData = obsHs(meanIndex);
                    meanHsByDistance(i,j,k) = nanmean(currentData);
                    
                end
            end
        end
    end
end
% timerM1 = toc
meanHsByDistance = meanHsByDistance(:);
% stdHsByBin = stdHsByBin(:);

indNaN = ~isnan(meanHsByDistance);

HSobs = meanHsByDistance(indNaN);


[LON, LAT, TIME] = ndgrid(gridLon,gridLat,gridTime(1:skipTime:end));
LON  = LON(:);
LAT  = LAT(:);
TIME = TIME(:);

LONobs = LON(indNaN);
LATobs = LAT(indNaN);
TIMEobs = TIME(indNaN);