function loadSatList = defineSatListESA(buoyTime,altPath)
% loadSatList = defineSatListESA(buoyTime,altPath)
%
% Altimeter	  Freq.-Band Latitude-coverage Initial-Date Final-Date
% ERS-1	      Ku	     -81.5 to 81.5	   01/08/1991	02/06/1996
% TOPEX	      Ku C	     -66 to 66	       25/09/1992	08/10/2005
% ERS-2	      Ku	     -81.5 to 81.5	   29/04/1995	11/05/2009
% GFO	      Ku	     -73 to 72	       07/06/2000	07/09/2008
% JASON-1	  Ku C	     -66.15 to 66.15   15/01/2002	21/06/2013
% ENVISAT	  Ku S	     -82 to 82	       14/05/2002	08/04/2012
% JASON-2	  Ku C	     -66.15 to 66.15   04/07/2008	Ongoing
% CRYOSAT-2	  Ku	     -88 to 88	       14/07/2010	Ongoing
% SARAL	      Ka	     -81.49 to 81.49   14/03/2013	Ongoing
% JASON-3	  Ku C	     -66.15 to 66.15   12/02/2016	Ongoing

%what to do with altimeters with 2 bands?
[glyph , ~] = giveGlyph;
timeOffset = datenum([1981 01 01 00 00 00]); 

satList  = dir(altPath);
satList = rmfield (satList,{'date','bytes','isdir','datenum'});
satList = satList(3:end);

for i = 1:length(satList)
    yearList = dir([altPath satList(i).name glyph 'v1.1' glyph]);
    yearList = rmfield (yearList,{'date','bytes','isdir','datenum'});
    yearList = yearList(3:end);
    
    monListStart = dir([altPath satList(i).name glyph 'v1.1' glyph...
        yearList(1).name]);
    monListStart = rmfield (monListStart,{'date','bytes','isdir','datenum'});
    monListStart = monListStart(3:end);
    
    monListEnd = dir([altPath satList(i).name glyph 'v1.1' glyph...
        yearList(end).name]);
    monListEnd= rmfield (monListEnd,{'date','bytes','isdir','datenum'});
    monListEnd = monListEnd(3:end);
    
    dayListStart = dir([altPath satList(i).name glyph 'v1.1' glyph...
        yearList(1).name glyph monListStart(1).name]);
    dayListStart = rmfield (dayListStart,{'date','bytes','isdir','datenum'});
    dayListStart = dayListStart(3:end);
    
    dayListEnd = dir([altPath satList(i).name glyph 'v1.1' glyph...
        yearList(end).name glyph monListEnd(end).name]);
    dayListEnd = rmfield (dayListEnd,{'date','bytes','isdir','datenum'});
    dayListEnd = dayListEnd(3:end);
        
    fileListStart = dir([altPath satList(i).name glyph 'v1.1' glyph...
        yearList(1).name glyph monListStart(1).name glyph...
        dayListStart(1).name]);
    fileListStart = rmfield (fileListStart,{'date','bytes','isdir','datenum'});
    fileListStart = fileListStart(3:end);
    
    fileListEnd = dir([altPath satList(i).name glyph 'v1.1' glyph...
        yearList(end).name glyph monListEnd(end).name glyph...
        dayListEnd(end).name]);
    fileListEnd = rmfield (fileListEnd,{'date','bytes','isdir','datenum'});
    fileListEnd = fileListEnd(3:end);
    
    altTimeStart = ncread([altPath satList(i).name glyph 'v1.1' glyph...
        yearList(1).name glyph monListStart(1).name glyph...
        dayListStart(1).name glyph fileListStart(1).name], 'time');

    altTimeEnd = ncread([altPath satList(i).name glyph 'v1.1' glyph...
        yearList(end).name glyph monListEnd(end).name glyph...
        dayListEnd(end).name glyph fileListEnd(end).name], 'time');
    
    satList(i).timeStart = altTimeStart(1)./(24*60*60) + timeOffset;    
    satList(i).timeEnd = altTimeEnd(end)./(24*60*60) + timeOffset;
end

% use model time to build a list of active altimeters
count1 = 0;
for i = 1:length(satList)
    satSpan = floor(satList(i).timeStart):ceil(satList(i).timeEnd);
    mdSpan  = floor(min(buoyTime)):ceil(max(buoyTime));
    
    if  intersect(satSpan,mdSpan)
        count1 = count1 + 1;
        loadSatList(count1) = satList(i);
    end
end

if count1 == 0;
    loadSatList = [];
end
