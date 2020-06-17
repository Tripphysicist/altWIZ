% TC - Altimeter Script
%
%
clear
close all
distCount=0;
for maxDistance = 250;
    %cmaxDistance = 500;
    maxTimeDiff = .5/24; %half hour
    distCount = distCount + 1;
    % load in raw altimeter data results for the model run
    %% IRMA
    % pairedDataPath = '/Users/tripp/D/Analysis/altimeterComparison/data/Irma_wis_saccs_basin_l1_SACS_T_none.mat';
    % mdPath = '/Users/tripp/D/Analysis/altimeterComparison/modelData/SACS/Irma/SACS_TP_0023_HIS_basin_l1_field_T/wis_saccs_basin_l1_SACS_T.nc';
    % stormPath = '/Users/tripp/D/Analysis/altimeterComparison/storm-tracks/2017/irma/irma_best_track_no_header.csv';
    %% MATTHEW
    % pairedDataPath = '/Users/tripp/D/Analysis/altimeterComparison/data/Matthew_wis_saccs_basin_l1_SACS_T_none.mat';
    % mdPath = '/Users/tripp/D/Analysis/altimeterComparison/modelData/SACS/Matthew/SACS_TP_0021_HIS_basin_l1_field/wis_saccs_basin_l1_SACS_T.nc';
    % stormPath = '/Users/tripp/D/Analysis/altimeterComparison/storm-tracks/2016/matthew/matthew_best_track_no_header.csv';
    %% BONNIE
%     pairedDataPath = '/Users/tripp/D/Analysis/altimeterComparison/data/atlantic2016MayBonnie.mat';
%     mdPath = '/Users/tripp/D/Analysis/altimeterComparison/modelData/tau6AtlanticYearRun/ww3.201605.nc';
%     stormTemp = getStormTrack(2016,'BONNIE');

    %% NICOLE
     pairedDataPath = '/Users/tripp/D/Analysis/altimeterComparison/data/atlantic_2016_Oct_Nicole.mat';
     mdPath = '/Users/tripp/D/Analysis/altimeterComparison/modelData/Nicole/ww3.201610.nc';
     stormTemp = getStormTrack(2016,'NICOLE');
    %%
    %load raw altimeter paired with model
    load(pairedDataPath,'pData');
    
    %load original model data
    try
        md.lat = ncread(mdPath, 'latitude');
    catch SACSfile
        md.lat = ncread(mdPath, 'lat');
    end
    md.lat = double(md.lat);
    try
        md.lonO = ncread(mdPath, 'longitude');
    catch
        md.lonO = ncread(mdPath, 'lon');
    end
    md.lon = mod(double(md.lonO),360); % 0 - 360 degrees
    md.lonO = md.lon;
    md.lonO(md.lon>180) = md.lon(md.lon>180)-360; %-180 - 180
    
    
    try
        md.hs = ncread(mdPath, 'hs');
    catch
        md.hs = ncread(mdPath, 'waveHs');
    end
    try
        md.gridStatus = ncread(mdPath, 'MAPSTA');
    catch
        md.gridStatus = ncread(mdPath, 'mask');
    end
    %becareful here, the code is not robust to different time stamps, may
    %need to make some changes yourself
    if exist('SACSfile')
        timeInSeconds = ncread(mdPath, 'time');
        md.time = double(timeInSeconds)/(60*60*24) + datenum([1970 01 01 00 00 00]);
    else
        md.time = ncread(mdPath, 'time') + datenum([1990 01 01 00 00 00]);
    end
    
    %load in storm track
    % stormData = dlmread(stormPath);
    % date  = num2str(stormData(:,1));
    % year  = str2num(date(:,1:4));
    % month = str2num(date(:,5:6));
    % day   = str2num(date(:,7:8));
    % hour  = str2num(date(:,9:10));
    % stormTime  = datenum([year month day hour zeros(length(date),1) zeros(length(date),1)]);
    % stormLon   = stormData(:,7)+360; %from negative W to E poisitive 0 - 360
    % stormLat   = stormData(:,6);
    stormTime = stormTemp.time;
    stormLon = stormTemp.lon + 360;
    stormLat = stormTemp.lat;
    %interp storm to model
    storm.time = interp1(stormTime, stormTime, md.time);
    storm.lon  = interp1(stormTime, stormLon, md.time);
    storm.lat  = interp1(stormTime, stormLat, md.time);
    %get rid on NaNs
    storm.time = storm.time(~isnan(storm.time));
    storm.lon = storm.lon(~isnan(storm.lon));
    storm.lat = storm.lat(~isnan(storm.lat));
    %keep data only close to the storm
    stormTimeLength = length(storm.time);
    % tic
    count = 0;
    for i = 1:stormTimeLength
        [xdummy,ydummy,distance] = latlon2xy(pData.lat,pData.lon,storm.lat(i),storm.lon(i));
        [nearStormIndex dum2] = find(distance <= maxDistance &...
            abs(storm.time(i) - pData.time) <= maxTimeDiff);
        if ~isempty(nearStormIndex)
            count = count + 1;
            stormObs(count).lat   = pData.lat(nearStormIndex);
            stormObs(count).lon   = pData.lon(nearStormIndex);
            stormObs(count).time  = pData.time(nearStormIndex);
            stormObs(count).altHs = pData.altHs(nearStormIndex);
            stormObs(count).mdHs  = pData.mdHs(nearStormIndex);
            %         plot(stormObs(count).time,stormObs(count).altHs,'.')
            %         hold on
            %         pause
            %         shg
            %         datetick('x')
        end
    end
    %New structure with only storm related paired data
    stormCat.lat = vertcat(stormObs(:).lat);
    stormCat.lon = vertcat(stormObs(:).lon);
    stormCat.time = vertcat(stormObs(:).time);
    stormCat.altHs = vertcat(stormObs(:).altHs);
    stormCat.mdHs = vertcat(stormObs(:).mdHs);
    %get rid of NaNs
    stormCat.lat   = stormCat.lat(~isnan(stormCat.altHs));
    stormCat.lon   = stormCat.lon(~isnan(stormCat.altHs));
    stormCat.time  = stormCat.time(~isnan(stormCat.altHs));
    stormCat.mdHs  = stormCat.mdHs(~isnan(stormCat.altHs));
    stormCat.altHs = stormCat.altHs(~isnan(stormCat.altHs));
    %get rid of bad data
    [goodModelIndex k]= find(stormCat.mdHs > 0);
    stormCat.lat   = stormCat.lat(goodModelIndex);
    stormCat.lon   = stormCat.lon(goodModelIndex);
    stormCat.time  = stormCat.time(goodModelIndex);
    stormCat.mdHs  = stormCat.mdHs(goodModelIndex);
    stormCat.altHs = stormCat.altHs(goodModelIndex);
    %sort by time
    [stormCat.time, timeIndex] = sort(stormCat.time);
    stormCat.lon = stormCat.lon(timeIndex);
    stormCat.lat = stormCat.lat(timeIndex);
    stormCat.mdHs = stormCat.mdHs(timeIndex);
    stormCat.altHs = stormCat.altHs(timeIndex);
    % simple plots
    lonLatTimeSeries(stormCat.lon,stormCat.lat,stormCat.time,stormCat.altHs,stormCat.mdHs)
    colorScatter(stormCat.altHs,stormCat.mdHs,100)
    L = Ocstatsp(stormCat.altHs,stormCat.mdHs,0.1)
    bias(distCount) = L.bias;
    R(distCount) = L.R;
    rmse(distCount) = L.rmse;
    %% interpolate storm track to paired data time
    stormLon = interp1(storm.time,storm.lon,stormCat.time);
    stormLat = interp1(storm.time,storm.lat,stormCat.time);
    stormLon = stormLon(~isnan(stormLon));
    stormLat = stormLat(~isnan(stormLat));
    
    for i=1:length(stormLon)-1
        [X(i),Y(i),D(i)] = latlon2xy(stormLat(i+1),stormLon(i+1),stormLat(i),stormLon(i));
    end
    X(length(stormLon)) = NaN;
    Y(length(stormLon)) = NaN;
    D(length(stormLon)) = NaN;
    stormDir = deg360(90 - (180/pi).*atan2(Y,X));
    [x, y, r] = latlon2xy(stormCat.lat,stormCat.lon,stormLat,stormLon);
    altAngle =  deg360(90 - (180/pi).*atan2(y,x));
    stormRefAngle = deg360(270-(altAngle + stormDir'));
    %
    stormCat.lon = deg180(stormCat.lon);
    minLON1 = min(stormCat.lon);
    minLON2 = min(deg180(storm.lon));
    maxLON1 = max(stormCat.lon);
    maxLON2 = max(deg180(storm.lon));
    minLAT1 = min(stormCat.lat);
    minLAT2 = min(storm.lat);
    maxLAT1 = max(stormCat.lat);
    maxLAT2 = max(storm.lat);
    
    minLON = min(minLON1,minLON2);
    maxLON = max(maxLON1,maxLON2);
    minLAT = min(minLAT1,minLAT2);
    maxLAT = max(maxLAT1,maxLAT2);
    
    %%
    theta = linspace(0,2*pi,1000);
    
    % figure
    % subplot(1,2,1)
    % plot(x,y,'o')
    % axis('equal')
    % hold on
    % plot(100*cos(theta),100*sin(theta),'--')
    % plot(250*cos(theta),250*sin(theta),'--')
    % plot(500*cos(theta),500*sin(theta),'--')
    % plot(750*cos(theta),750*sin(theta),'--')
    % grid on
    %
    % subplot(1,2,2)
    % plot(r.*cosd(stormRefAngle),r.*sind(stormRefAngle),'o')
    % axis('equal')
    % hold on
    % theta = linspace(0,2*pi,1000);
    % plot(100*cos(theta),100*sin(theta),'--')
    % plot(250*cos(theta),250*sin(theta),'--')
    % plot(500*cos(theta),500*sin(theta),'--')
    % plot(750*cos(theta),750*sin(theta),'--')
    % grid on
    %%
    
    figure
    subplot(2,1,1)
    m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
    m_plot(stormCat.lon,stormCat.lat,'.','markerSize',10)
    hold on
    m_coast('patch',[.8 .8 .8]);
    m_plot(deg180(storm.lon),storm.lat,'ko-')
    m_grid('box','fancy','tickdir','out');
    title('Hurricane and Altimeter Tracks')
    fontsize(16,16,16,16)
    grid on
    
    plot(stormLon,stormLat,'ko-')
    hold on
    plot(stormCat.lon,stormCat.lat,'o')
    subplot(2,1,2)
    plot(100*cos(theta),100*sin(theta),'--')
    hold on
    axis('equal')
    plot(200*cos(theta),200*sin(theta),'--')
    plot(300*cos(theta),300*sin(theta),'--')
%    plot(400*cos(theta),400*sin(theta),'--')
%    plot(500*cos(theta),500*sin(theta),'--')
%    plot(600*cos(theta),600*sin(theta),'--')
%    plot(700*cos(theta),700*sin(theta),'--')
%    plot(800*cos(theta),800*sin(theta),'--')
%    plot(900*cos(theta),900*sin(theta),'--')
%    plot(1000*cos(theta),1000*sin(theta),'--')
    % plot(r.*cosd(stormRefAngle),r.*sind(stormRefAngle),'o');
    scatter(r.*cosd(stormRefAngle),r.*sind(stormRefAngle),1000,stormCat.altHs-stormCat.mdHs,'.');
    xlabel('distance from eye [km]')
    ylabel('distance from eye [km]')
    grid on
    h = colorbar;
    ylabel(h, 'wave height difference [m]')
    axis('square')
    shg
    fontsize(20,20,20,20)
end

%%
theta = linspace(0,2*pi,1000);

plot(100*cos(theta),100*sin(theta),'--')
hold on
axis('equal')
plot(200*cos(theta),200*sin(theta),'--')
plot(300*cos(theta),300*sin(theta),'--')
xlabel('distance from eye [km]')
ylabel('distance from eye [km]')
grid on
h = colorbar;
ylabel(h, 'wave height difference [m]')
axis('square')
shg
fontsize(20,20,20,20)

for i = 1:length(stormObs)
    if isempty(stormObs(i).hs)
        continue
    end    
plot(stormObs(i).distance2center.*cosd(stormObs(i).stormRefAngle),stormObs(i).distance2center.*sind(stormObs(i).stormRefAngle),'o');
% scatter(r.*cosd(stormRefAngle),r.*sind(stormRefAngle),1000,stormCat.altHs-stormCat.mdHs,'.');
shg
end
%%
blockSize = 10; %[degrees]
minLON = min(pData.lon);
maxLON = max(pData.lon);
minLAT = min(pData.lat);
maxLAT = max(pData.lat);
statGridEdgesx = minLON:blockSize:maxLON;
statGridEdgesLat = minLAT:blockSize:maxLAT;
statGridCenterLon = statGridEdgesLon(1:end-1) + diff(statGridEdgesLon)/2;
statGridCenterLat = statGridEdgesLat(1:end-1) + diff(statGridEdgesLat)/2;

for i = 1:length(statGridCenterLon)
    for j = 1:length(statGridCenterLat)
        [index dumy] = find(pData.lon >= statGridEdgesLon(i) & pData.lon <= statGridEdgesLon(i+1) &...
            pData.lat >= statGridEdgesLat(j) & pData.lat <= statGridEdgesLat(j+1));
        meanAlt(i,j) = nanmean(pData.altHs(index));
        meanMod(i,j) = nanmean(pData.mdHs(index));        
        dHS = pData.altHs(index)-pData.mdHs(index);
        bias(i,j) = meanMod(i,j) - meanAlt(i,j);
        nbias(i,j) = bias(i,j)./meanAlt(i,j);
        rmse(i,j) = sqrt(mean(dHS.^2));
        nrmse(i,j) = sqrt(mean(dHS.^2)/meanAlt(i,j).^2);
        scatIn(i,j) = rmse(i,j)/meanAlt(i,j);      
    end
end
