function loadSatList = defineSatList(buoyTime)

% Altimeter	  Freq.-Band Latitude-coverage Initial-Date Final-Date
% GEOSAT	  Ku	     -73 to 72	       31/03/1985	31/12/1989
% ERS-1	      Ku	     -81.5 to 81.5	   01/08/1991	02/06/1996
% TOPEX	      Ku C	     -66 to 66	       25/09/1992	08/10/2005
% ERS-2	      Ku	     -81.5 to 81.5	   29/04/1995	11/05/2009
% GFO	      Ku	     -73 to 72	       07/06/2000	07/09/2008
% JASON-1	  Ku C	     -66.15 to 66.15   15/01/2002	21/06/2013
% ENVISAT	  Ku S	     -82 to 82	       14/05/2002	08/04/2012
% JASON-2	  Ku C	     -66.15 to 66.15   04/07/2008	Ongoing
% CRYOSAT-2	  Ku	     -88 to 88	       14/07/2010	Ongoing
% HY-2A	      Ku C	     -81 to 80	       01/10/2011	06/06/2018
% SARAL	      Ka	     -81.49 to 81.49   14/03/2013	Ongoing
% JASON-3	  Ku C	     -66.15 to 66.15   12/02/2016	Ongoing
% SENTINEL-3A Ku C	     -78 to 81	       01/03/2016	Ongoing

%what to do with altimeters with 2 bands?

altPath = 'd:\Datasets\Satellite\Altimeter\Ribal_Young_2019\';
satList  = dir(altPath);
satList = rmfield (satList,{'date','bytes','isdir','datenum'});
satList = satList(3:end);

for i = 1:length(satList)
    satFileList = dir([altPath satList(i).name]);
    satFileList = rmfield (satFileList,{'date','bytes','isdir','datenum'});
    satFileList = satFileList(3:end);
    info = ncinfo([altPath satList(i).name '\' satFileList(1).name]);
    %    timeOffsetMeta(i) = info.Variables(1).Attributes(3);
    altTime = ncread([altPath satList(i).name '\' satFileList(1).name],'TIME');
    satList(i).timeStart = altTime(1)+datenum([1985 01 01 00 00 00]);
    satList(i).timeEnd = altTime(end)+datenum([1985 01 01 00 00 00]);
end

% use model time to neglect some of the altimeter data and build a list
count1 = 0;
for i = 1:length(satList)
    satSpan = floor(satList(i).timeStart):ceil(satList(i).timeEnd);
    mdSpan  = floor(min(buoyTime)):ceil(max(buoyTime));
    
    if  intersect(satSpan,mdSpan)
        count1 = count1 +1;
        loadSatList(count1) = satList(i);
    end
end
