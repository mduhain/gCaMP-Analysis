%extrace the red frames to make a rfp trace for all files in a directory
function [rfpTrace] = rfpTraceFromAllFiles()
    inputFiles = dir(fullfile('*.tif')); %find all .tif files in dir
    fileNames = {inputFiles.name};
    imInfo = imfinfo(fileNames{1}); %get info on first tif
    numFrames = length(inputFiles); %Num of available frames
    rfpTrace = zeros(numFrames,1);
    disp("extracting RFP frames...");
    tic;
    for k = 1 : numFrames
        try
            thisFileName = fileNames{k};
            rfpImg = imread(thisFileName,2);
            rfpTrace(k) = mean(rfpImg,'all');
        catch
            disp(strcat("Error on frame # ",num2str(k)));
        end
    end
    toc;
end