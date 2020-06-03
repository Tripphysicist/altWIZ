function distance = latlon2dist(lat1,lon1,lat2,lon2);
% function distance = latlon2dist(lat1,lon1,lat2,lont2);
% Haversine distance on Earth
R = 6371; % Earth's radius in km
lat1=lat1*pi/180;
lon1=lon1*pi/180;
lat2=lat2*pi/180;
lon2=lon2*pi/180;
deltaLat=lat2-lat1;
deltaLon=lon2-lon1;
a=sin((deltaLat)/2).^2 + cos(lat1).*cos(lat2) .* sin(deltaLon/2).^2;
c=2*atan2(sqrt(a),sqrt(1-a));
distance = R.*c;    %Haversine distance    