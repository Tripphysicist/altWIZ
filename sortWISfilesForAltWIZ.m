% function sortWISfilesForAltWIZ(modelPath)
% not repeatable, so switched to script
%% application #1, all field grids and field files
% modelPath = '';
[glyph , ~] = giveGlyph;

if modelPath(end)~=glyph
    modelPath = [modelPath glyph];
end

modelDir = dir(modelPath);
modelDir = modelDir(3:end);

grids = dir([modelPath modelDir(1).name]);
grids = grids(3:end);

for i = 1:numel(grids)
    mkdir([modelPath grids(i).name])
end

for i = 1:numel(modelDir)
    for j = 1:numel(grids)
        currentPath = [modelPath modelDir(i).name glyph grids(j).name glyph];
        currentDir = dir(currentPath);
        currentDir = currentDir(3:end);
        if numel(currentDir)==1
            untar([currentPath currentDir(1).name],currentPath)
        end
        currentDir = dir(currentPath);
        currentDir = currentDir(3:end);
        fileName = [currentPath currentDir(2).name];
        destination = [modelPath grids(j).name];
        movefile(fileName,destination);
    end
end


%% application # 2, level 1 field files
modelPath = '/Users/tripp/D/Analysis/altimeterComparison/modelData/atlanticTropicalCyclones/WIS/';

[glyph , ~] = giveGlyph;

if modelPath(end)~=glyph
    modelPath = [modelPath glyph];
end

modelDir = dir(modelPath);
modelDir = modelDir(3:end);

for i = 1:numel(modelDir)
    filePath = [modelPath modelDir(i).name];
    untar(filePath)
end

