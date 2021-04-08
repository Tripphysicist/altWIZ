function [distance, bearing] = latlon2dist(lat1,lon1,lat2,lon2)
% function [distance, bearing] = latlon2dist(lat1,lon1,lat2,lont2);
% Haversine distance on Earth
% Initial Bearing Formula
R = 6371; % Earth's radius in km
lat1=lat1*pi/180; %radial lat
lon1=lon1*pi/180; %radial lon
lat2=lat2*pi/180;
lon2=lon2*pi/180;
deltaLat=lat2-lat1; %diff lat
deltaLon=lon2-lon1; %diff lon
a=sin((deltaLat)/2).^2 + cos(lat1).*cos(lat2) .* sin(deltaLon/2).^2;
c=2*atan2(sqrt(a),sqrt(1-a));
distance = R.*c;    %Haversine distance

% Initial bearing formula (https://www.movable-type.co.uk/scripts/latlong.html)
% For small distances, the initial and final bearing are similar
y = sin(deltaLon).*cos(lat2);
x = cos(lat1).*sin(lat2) - sin(lat1).*cos(lat2).*cos(deltaLon);
theta = atan2(y, x);
bearing = ((theta.*(180./pi) + 360)-180); %flopped the initial and end points need to -180