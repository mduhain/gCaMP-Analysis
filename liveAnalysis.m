% LIVE ANALYSIS
% first pass of script to complete analysis on folder of .tif 2P images.
% 2023-02-02


%% Constants

frameRate = 30; %capture rate of 2P
preStimWin = 1; %seconds
postStimWin = 6; %seconds


%% File Identification

% file names
inputFiles = dir(fullfile('*.tif')); %find all .tif files in dir
fileNames = {inputFiles.name};

% information of frame(s)
imInfo = imfinfo(fileNames{1}); %get info on first tif
Xpx = imInfo(1).Width; %X pixels
Ypx = imInfo(1).Height; %Y pixels
numFrames = length(inputFiles); %Num of available frames
rfpTrace = zeros(numFrames,1);
rfpImg = zeros(Xpx,Ypx);
gfpImg = zeros(Xpx,Ypx);


%% Extract Stim Times (RFP)

disp("calculating stim times...");
tic;
for k = 1 : numFrames
    thisFileName = fileNames{k};
    rfpStack(:,:,k) = imread(thisFileName,1);
    rfpImg = imread(thisFileName,2);
    rfpTrace(k) = mean(rfpImg,'all');
end


%% Process stim times

dtRfp = detrend(rfpTrace); % detrend data
meanRfp = mean(dtRfp); % identify mean
rfpThresh = meanRfp + 2*std(rfpTrace); % threshold mean + 2std
% find frame numbers of above threshold flashes
I = find(dtRfp > rfpThresh);
dI = diff(I); dI = [0; dI]; %frames between events
I(dI == 1) = []; %remove any incidents within 1 frame of another (noise)
I(find(I<60)) = []; %#ok<FNDSB> %remove any events during first 60 frames
toc;


%% Loop through and extract GFP traces

disp("extracting gfp signals...");
tic;
gfpTraces = zeros(4,frameRate*(preStimWin+postStimWin),size(I,1));
% [4 quadrants, window frames, number of events]

for j = 1 : size(I,1)
    firstF = I(j) - (frameRate*preStimWin-1);
    lastF = I(j) + frameRate*postStimWin;
    frameList = [firstF:lastF];
    if lastF > numFrames
        break;
    end
    %loop through and extract gfp frames in this window
    for h = 1 : size(frameList,2)
        thisFileName = fileNames{frameList(h)};
        gfpImg = imread(thisFileName,1);
        gfpTraces(1,h,j) = mean(gfpImg(1:Xpx/2,1:Ypx/2),'all'); %upper left
        gfpTraces(2,h,j) = mean(gfpImg(1:Xpx/2,Ypx/2+1:end),'all'); %upper right
        gfpTraces(3,h,j) = mean(gfpImg(Xpx/2+1:end,1:Ypx/2),'all'); %lower left
        gfpTraces(4,h,j) = mean(gfpImg(Xpx/2+1:end,Ypx/2+1:end),'all'); %lower right
    end
end

xs = [-1*preStimWin+1/frameRate:1/frameRate:postStimWin];
toc;

%Plotting section
meanTrace1 = squeeze(mean(gfpTraces(1,:,:),3));
meanTrace1 = meanTrace1 - mean(meanTrace1(1:30));
semTrace1 = std(squeeze(gfpTraces(1,:,:))')/sqrt(size(gfpTraces,3));

meanTrace2 = squeeze(mean(gfpTraces(2,:,:),3));
meanTrace2 = meanTrace2 - mean(meanTrace2(1:30));
semTrace2 = std(squeeze(gfpTraces(2,:,:))')/sqrt(size(gfpTraces,3));

meanTrace3 = squeeze(mean(gfpTraces(3,:,:),3));
meanTrace3 = meanTrace3 - mean(meanTrace3(1:30));
semTrace3 = std(squeeze(gfpTraces(3,:,:))')/sqrt(size(gfpTraces,3));

meanTrace4 = squeeze(mean(gfpTraces(4,:,:),3));
meanTrace4 = meanTrace4 - mean(meanTrace4(1:30));
semTrace4 = std(squeeze(gfpTraces(4,:,:))')/sqrt(size(gfpTraces,3));

figure;
shadedErrorBar(xs,meanTrace1,semTrace1); hold on;
shadedErrorBar(xs,meanTrace2,semTrace2);
shadedErrorBar(xs,meanTrace3,semTrace3);
shadedErrorBar(xs,meanTrace4,semTrace4);
ca = gca;
plot([0 0],[ca.YLim(1) ca.YLim(2)],'r-'); 
title("Quick PSTH, 4 quadrants in field (Mean + SEM)");
xlabel("Time (sec.)");
ylabel("Norm. DeltaF/F");
hold off;









