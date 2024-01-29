function [folderNames, fileNames] = analyzeDir()
% analyzeDir.m
% [folderNames, fileNames] = analyzeDir()
%   - Works in a directory for extraxt file(s) and folder(s) names.
    myDir = dir;
    folderNames = repmat("",size(myDir,1),1);
    fileNames = repmat("",size(myDir,1),1);
    n=1;
    while n <= size(myDir,1)
        if strcmp(myDir(n).name,'.') 
            myDir(n) = []; %clear entry
            continue;
        elseif strcmp(myDir(n).name,'..')
            myDir(n) = []; %clear entry
            continue;
        elseif strcmp(myDir(n).name,'Analysis')
            myDir(n) = []; %clear entry
            continue;
        elseif myDir(n).isdir == 1
            folderNames(n) = string(myDir(n).name);
            n=n+1;
        elseif myDir(n).isdir == 0
            fileNames(n) = string(myDir(n).name);
            n=n+1;
        else
            n=n+1;
        end
    end
    fileNames(strcmp(fileNames,"")) = [];
    folderNames(strcmp(folderNames,"")) = [];
end
