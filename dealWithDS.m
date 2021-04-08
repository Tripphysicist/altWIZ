function newDir = dealWithDS(oldDir);
% DEAL WITH GD .DS files that MAC creates
if    contains(oldDir(1).name,'.DS_Store')
    newDir = oldDir(2:end);
end
