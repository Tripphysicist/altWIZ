% altWW3CompPlotting
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Script to make plots for comparing WW3 output to altimeter measurements %
% of significant wave height and wind speed                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Author:       Clarence Olin Collins III, Ph.D.                          %
% Affiliation:  ERDC - FRF                                                %
% created:      02/21/2020                                                %
% version:      1.0                                                       %
% updates:                                                                %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %


%% PART ?: calculate statistics

%% PART ?: make simple plots
figure
subplot(1,3,1)
plot(TIMEobsCat,HSobsCat,'.')
hold on
plot(TIMEobsCat,HSmdCat,'.')
grid on
datetick
subplot(1,3,2)
plot(LONobsCat,HSobsCat,'.')
hold on
plot(LONobsCat,HSmdCat,'.')
grid on
subplot(1,3,3)
plot(LATobsCat,HSobsCat,'.')
hold on
plot(LATobsCat,HSmdCat,'.')
grid on

packfig(1,3)

%subplot(1,2,2)
L = Ocstatsp(HSobsCat,HSmdCat,0.1)

% distributions, histograms
% q-q plots
% scatter plot
% taylor diagram

%% map of stats

minLON = min(LONobsCat);
maxLON = max(LONobsCat);
minLAT = min(LATobsCat);
maxLAT = max(LATobsCat);
blockSize = 1; %[degrees]
statGridEdgesLon = minLON:blockSize:maxLON;
statGridEdgesLat = minLAT:blockSize:maxLAT;
statGridCenterLon = statGridEdgesLon(1:end-1) + diff(statGridEdgesLon)/2;
statGridCenterLat = statGridEdgesLat(1:end-1) + diff(statGridEdgesLat)/2;

for i = 1:length(statGridCenterLon)
    for j = 1:length(statGridCenterLat)
        [index dumy] = find(LONobsCat >= statGridEdgesLon(i) & LONobsCat <= statGridEdgesLon(i+1) &...
            LATobsCat >= statGridEdgesLat(j) & LATobsCat <= statGridEdgesLat(j+1));
        dHS = HSobsCat(index)-HSmdCat(index);
        bias(i,j) = mean(dHS);
        rmse(i,j) = sqrt(mean(dHS.^2));
    end
end


figure
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,rmse')
shading('interp')
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
title('rmse')
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('jet')
caxis([-2 2])

figure
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,bias')
shading('interp')
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
title('bias')
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('lansey')
caxis([0 3])

%% months and seasons
date = datevec(TIMEobsCat);

for month = 1:12
    monthIndex = date(:,2) == month;
    eval(['LON' num2str(month) '= LONobsCat(monthIndex);'])
    eval(['LAT' num2str(month) '= LATobsCat(monthIndex);'])
    eval(['TIME' num2str(month) '= TIMEobsCat(monthIndex);'])
    eval(['HSobs' num2str(month) '= HSobsCat(monthIndex);'])
    eval(['HSmd' num2str(month) '= HSmdCat(monthIndex);'])           
end

for i = 1:length(statGridCenterLon)
    for j = 1:length(statGridCenterLat)
        for month = 1:12
            LON = eval(['LON' num2str(month)]);
            LAT = eval(['LAT' num2str(month)]);
            HSobs = eval(['HSobs' num2str(month)]);
            HSmd = eval(['HSmd' num2str(month)]);
            [index dumy] = find(LON >= statGridEdgesLon(i) & LON <= statGridEdgesLon(i+1) &...
            LAT >= statGridEdgesLat(j) & LAT <= statGridEdgesLat(j+1));
            dHS = HSobs(index)-HSmd(index);
            eval(['bias' num2str(month) '(i,j) = mean(dHS);'])
            eval(['rmse' num2str(month) '(i,j) = sqrt(mean(dHS.^2));'])
        end
    end
end

LONs1 = [LON1; LON2; LON3];
LONs2 = [LON4; LON5; LON6];
LONs3 = [LON7; LON8; LON9];
LONs4 = [LON10; LON11; LON12];
LATs1 = [LAT1; LAT2; LAT3];
LATs2 = [LAT4; LAT5; LAT6];
LATs3 = [LAT7; LAT8; LAT9];
LATs4 = [LAT10; LAT11; LAT12];
HSobss1 = [HSobs1; HSobs2; HSobs3];
HSobss2 = [HSobs4; HSobs5; HSobs6];
HSobss3 = [HSobs7; HSobs8; HSobs9];
HSobss4 = [HSobs10; HSobs11; HSobs12];
HSmds1 = [HSmd1; HSmd2; HSmd3];
HSmds2 = [HSmd4; HSmd5; HSmd6];
HSmds3 = [HSmd7; HSmd8; HSmd9];
HSmds4 = [HSmd10; HSmd11; HSmd12];

for i = 1:length(statGridCenterLon)
    for j = 1:length(statGridCenterLat)
        for s = 1:4
            LON = eval(['LONs' num2str(s)]);
            LAT = eval(['LATs' num2str(s)]);
            HSobs = eval(['HSobss' num2str(s)]);
            HSmd = eval(['HSmds' num2str(s)]);
            [index dumy] = find(LON >= statGridEdgesLon(i) & LON <= statGridEdgesLon(i+1) &...
                LAT >= statGridEdgesLat(j) & LAT <= statGridEdgesLat(j+1));
            dHS = HSobs(index)-HSmd(index);
            eval(['biasS' num2str(s) '(i,j) = mean(dHS);'])
            eval(['rmseS' num2str(s) '(i,j) = sqrt(mean(dHS.^2));'])
        end
    end
end
%%
figure
subplot(2,2,1)
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,rmseS1')
shading('interp')
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('jet')
caxis([-2 2])
subplot(2,2,2)
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,rmseS2')
shading('interp')
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('jet')
caxis([-2 2])
subplot(2,2,3)
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,rmseS3')
shading('interp')
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('jet')
caxis([-2 2])
subplot(2,2,4)
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,rmseS4')
shading('interp')
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('jet')
caxis([-2 2])
packfig(2,2)

figure
subplot(2,2,1)
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,biasS1')
shading('interp')
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
title('bias')
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('lansey')
%caxis([0 3])
subplot(2,2,2)
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,biasS2')
shading('interp')
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('lansey')
%caxis([0 3])
subplot(2,2,3)
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,biasS3')
shading('interp')
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('lansey')
%caxis([0 3])
subplot(2,2,4)
m_proj('miller','lon',[minLON maxLON],'lat',[minLAT maxLAT]);
m_pcolor(statGridCenterLon,statGridCenterLat,biasS4')
shading('interp')
m_coast('patch',[.8 .8 .8]);
m_grid('box','fancy','tickdir','out');
fontsize(16,16,16,16)
h = colorbar;
ylabel(h, '[m]')
colormap('lansey')
%caxis([0 3])