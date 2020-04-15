function [LONobs LATobs TIMEobs HSobs] = boxMethod(obsLon, obsLat, obsTime, obsHs, gridLon, gridLat, maxTimeDiff, gridStatus, minNumberObs) 
% [LONobs LATobs TIMEobs HSobs] = boxMethod(obsLon, obsLat, obsTime, obsHs, gridLon, gridLat,maxTimeDiff,gridStatus,minNumberObs) 

% %for debugging
% obsLon = LONobsNA;
% obsLat = LATobsNA;
% obsTime = TIMEobsNA;
% obsHs = HSobsNA;
% gridLon = mdTest.lon;
% gridLat = mdTest.lat;
% gridStatus = mdTest.gridStatus;

% There is an issue with longitude jumps 359:0 at the prime meirdian, need
% to switch to -180:180
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

%define bin spacing
lonBinSpacing  = .5; %degrees
latBinSpacing  = .5; %degrees
timeBinSpacing = maxTimeDiff*2; %days

%define bin edges

lonEdgeMin = (floor(min(obsLon)*2))/2;
lonEdgeMax = (ceil(max(obsLon)*2))/2;
lonEdges  = lonEdgeMin-lonBinSpacing/2:lonBinSpacing:lonEdgeMax+lonBinSpacing/2;

latEdgeMin = (floor(min(obsLat)*2))/2;
latEdgeMax = (ceil(max(obsLat)*2))/2;
latEdges  = latEdgeMin-latBinSpacing/2:latBinSpacing:latEdgeMax+latBinSpacing/2;

timeEdgeStart = (floor(min(obsTime)*(24*(60./(timeBinSpacing*24*60)))))/(24*(60./(timeBinSpacing*24*60))); %round down to closest time space unit
timeEdgeEnd = (ceil(max(obsTime)*(24*(60./(timeBinSpacing*24*60)))))/(24*(60./(timeBinSpacing*24*60))); %round up to closest time space unit
timeEdges = timeEdgeStart-maxTimeDiff:timeBinSpacing:timeEdgeEnd+maxTimeDiff;%+timeBinSpacing;

%get bin indexes
lonBinidx  = discretize(obsLon,lonEdges); %dependency
latBinidx  = discretize(obsLat,latEdges); %dependency
timeBinidx = discretize(obsTime,timeEdges); %dependency

%get bin centers
lonBinCenter  = lonEdges(1:end-1) + diff(lonEdges)/2;
latBinCenter  = latEdges(1:end-1) + diff(latEdges)/2;
timeBinCenter = timeEdges(1:end-1) + diff(timeEdges)/2;

%lengths
lonLength = length(lonBinCenter);
latLength = length(latBinCenter);
timeLength = length(timeBinCenter);

%allocate
meanHsByBin = NaN(lonLength,latLength,timeLength);
stdHsByBin =  NaN(lonLength,latLength,timeLength);

%loop through each dimension to average each cubic bin
% tic
for i = 1:lonLength;
    if sum(gridLon == lonBinCenter(i))
        for j = 1:latLength;
            if sum(gridLat == latBinCenter(j))
                if gridStatus(find(gridLon==lonBinCenter(i)),find(gridLat==latBinCenter(j)))==1 %check grid status
                    for k = 1:timeLength;
                        currentDataBin = obsHs(lonBinidx == i & latBinidx == j & timeBinidx == k);
                        if length(currentDataBin) >= minNumberObs
                            meanHsByBin(i,j,k) = nanmean(currentDataBin);
                        end
                        %            stdHsByBin(i,j,k)  = nanstd(currentDataBin); Any other
                        %            function
                    end
                end
            end
        end
    end
end
% timerM2 = toc
meanHsByBin = meanHsByBin(:);
% stdHsByBin = stdHsByBin(:);

indNaN = ~isnan(meanHsByBin);

HSobs = meanHsByBin(indNaN);

[LON, LAT, TIME] = ndgrid(lonBinCenter,latBinCenter,timeBinCenter);

%if we converted to negative longitudes, convert back to 0-365

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
