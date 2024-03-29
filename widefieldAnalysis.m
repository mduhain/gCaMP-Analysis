%function [gfpStack] = widefieldAnalysis(~)

%% Constants

frameRate = 30; %capture rate of 2P
preStimWin = 1; %seconds
postStimWin = 6; %seconds
subdivSize = 128; %number of sections to divide window into


%% File Identification

% file names

inputFiles = dir(fullfile('*.tif')); %find all .tif files in dir
fileNames = {inputFiles.name};


% information of frame(s)
[gfpStack, rfpTrace, scanInfo] = loadExp(fileNames{2});
% imInfo = imfinfo(fileNames{1}); %get info on first tif
% Xpx = imInfo(1).Width; %X pixels
% Ypx = imInfo(1).Height; %Y pixels
% numFrames = floor(length(inputFiles)); %Num of available frames
% rfpTrace = zeros(numFrames,1);
% gfpStack = zeros(subdivSize,subdivSize,numFrames);
% rfpImg = zeros(Xpx,Ypx);
% gfpImg = zeros(Xpx,Ypx);   


%% Extract image data

% disp("Extracting image data...");
% tic;
% parfor k = 1 : numFrames
%     thisFileName = fileNames{k};
%     rfpImg = imread(thisFileName,2);
%     gfpImg = imread(thisFileName,1);
%     rfpTrace(k) = mean(rfpImg,'all');
%     %index through voxels
%     xV = Xpx/subdivSize;
%     yV = Ypx/subdivSize;
%     for nx = 1 : subdivSize
%         xs = (nx*xV-(xV-1)):(nx*xV);
%         for ny = 1 : subdivSize
%             ys = (ny*yV-(yV-1)):(ny*yV);
%             gfpStack(nx,ny,k) = mean(gfpImg(xs,ys),'all');
%         end
%     end
% end
% gfpMean = squeeze(mean(gfpStack,3));
% toc;


%% Subdivide the gfpStack into voxel sections
disp("Extracting trace from each voxel...");
tic;
[Xpx,Ypx,numFrames] = size(gfpStack);
xV = Xpx/subdivSize;
yV = Ypx/subdivSize;
gfpStackVox = zeros(subdivSize,subdivSize,numFrames);
for nx = 1 : subdivSize
    xs = (nx*xV-(xV-1)):(nx*xV);
    for ny = 1 : subdivSize
        ys = (ny*yV-(yV-1)):(ny*yV);
        gfpStackVox(nx,ny,:) = mean(gfpStack(xs,ys,:),[1 2]);
    end
end
toc;


%% Process stim times
disp("Processing Stim Times...");
tic;
dtRfp = detrend(rfpTrace); % detrend data
meanRfp = mean(dtRfp); % identify mean
rfpThresh = meanRfp + 2*std(rfpTrace); % threshold mean + 2std
% find frame numbers of above threshold flashes
I = find(dtRfp > rfpThresh);
dI = diff(I); dI = [0; dI]; %frames between events
I(dI == 1) = []; %remove any incidents within 1 frame of another (noise)
I(find(I<60)) = []; %#ok<FNDSB> %remove any events during first 60 frames
toc;


%% Calculate responsivity (dF/F)

disp("Calculating dF/F...");
tic;
responseMap = zeros(subdivSize,subdivSize);
fprintf("Loading ");
for nx = 1 : subdivSize
    for ny = 1 : subdivSize
        postStimVals = zeros(length(I),1);
        thisTrace = squeeze(mean(gfpStack(nx,ny,:),[1 2]));
        for ni = 1 : length(I)
            baseline = mean(thisTrace(I(ni)-(frameRate-1):I(ni)));
            postStimVals(ni) = mean(thisTrace(I(ni):I(ni)+frameRate)) - baseline;
        end
        responseMap(nx,ny) = mean(postStimVals);
    end
    if rem(nx,floor(subdivSize/10)) == 0
        fprintf(". ");
    end
end
fprintf("\n");
toc;


%% Plotting Section

figure;
imagesc(responseMap);
title("dF/F per Voxel (8x8px)");
colorbar;

%% Save Image

imwrite(rescale(responseMap),'response_map.png');


%end