% checkExpData
% scripe checkExpData, finds all data files in the current folder, and
% check their duration and imaging depth.

format short
clc
[~,fileNames] = analyzeDir();
for n = 1 : length(fileNames)
    load(fileNames(n))
    expDur = Data.expStopTime - Data.expStartTime;
    disp(strcat("Checking ",fileNames(n)));
    disp(strcat(num2str(expDur)," seconds (Exp Duration)"));
    disp(strcat(num2str(Data.imagingDepth)," um imaging depth."));
    disp(strcat(num2str(length(unique(Data.stimONTimes)))," unique stim ON times (",...
        num2str(Data.stimONTimes(1)),'-',num2str(Data.stimONTimes(end)),')'));
    disp(strcat(num2str(length(unique(Data.stimOnFrameNum)))," unique stim frame numbers (",...
        num2str(Data.stimOnFrameNum(1)),'-',num2str(Data.stimOnFrameNum(end)),')'));
    disp('---------------------------------------------------------------');
end