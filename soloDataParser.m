%solo Data Parser

% id directory
startingDir = pwd;

%find im2p files and Z000 depth folders
[depths, fileNames] = analyzeDir();

%% Match variables from file path
bDir = strsplit(startingDir,'\');
for n=1:length(bDir)
    if ~isempty(regexp(bDir{n},'[m]\d{3}', 'once')) %'m000'
        mouseNum = bDir{n};
    end
    if ~isempty(regexp(bDir{n},'(\d{4})-(\d{2})-(\d{2})', 'once')) %'yyy-mm-dd'
       ddate = bDir{n};
    end
    if ~isempty(regexp(bDir{n},'[exp]\d{3}', 'once')) %'exp000'
        expName = bDir{n};
    end
end

%% Check to make sure no im2p files exist already
if any(contains(fileNames,'im2p'))
    disp(strcat("IM2P found, skipping ",mouseDir,'\',dates(nd),'\',expNames(ne)));
    return
end     

%% Parse Data
for k = 1 : length(depths) %loop through different depths
    %change to Z000 directory (depth)
    cd(strcat(startingDir,'\',depths(k)));

    %extract gfp signals, resize 512->256 && 3frame Avg
    gfpStack = gfpFromDir();

    %extract rfp signal (trace)
    rfpTrace = rfpTraceFromAllFiles();

    %save parameters into im2p
    a = im2p;

    % back out 1 directory to 'Exp000'
    cd(startingDir);

    % add other params. rfp imgs, vasculature 
    a.gfp3FrameAvg = gfpStack;
    a.rfpExpTrace = rfpTrace;

    % save updated im2p
    fileOutName = strcat(mouseNum,'-',ddate,'-',expName,'-',depths(k),'-im2p.mat');
    save(fileOutName,'a','-v7.3');

    %clear workspace
    clear a gfpStack rfpTrace

end

%% Nested functions
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
%--------------------------------------------------------------------------


