% altimeterBuoyComparisonPlots

%% %% PART ?: collocation 
Ocstatsp(buoyTest.hs(indNaN),HSobs,.12)

% F = scatteredInterpolant(LONobs,LATobs,TIMEobs,HSobs,'linear','none');
% HSalt = F(papaLon,papaLat,papaTime);
% indNaN = ~isnan(HSalt);
% HSalt = HSalt(indNaN);

% phase 1: interpolation

%interp3d
%{
        LONobs  = vertcat(obs(:).lon);
        LATobs  = vertcat(obs(:).lat);
        TIMEobs = vertcat(obs(:).time);
        HSobs   = vertcat(obs(:).hsKcal);
        HSmd = interp1(papaTime, mdTest.hs, TIMEobs);
        
        indNaN = ~isnan(HSmd);
        HSmd = HSmd(indNaN);
        
        HSobs  = HSobs(indNaN);
        LONobs = LONobs(indNaN);
        LATobs = LATobs(indNaN);
        TIMEobs = TIMEobs(indNaN);
%}

% phase 2: kd tree
figure
plot(buoyTest.time, buoyTest.hs,'.')
hold on
errorbar(TIMEobs, HSobs, HSobsSTD,'.k')
grid on
datetick('x','mm-yy')