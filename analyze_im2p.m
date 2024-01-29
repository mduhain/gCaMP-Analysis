%Analysis script
%
%
%
%

% [512 x 512] @30fps .tif images -->  [256 x 256] @10fps uint16 array
compressImgs = 1; % 1 for yes, 0 for no.

%% Load im2p file
currDir = pwd;
disp("Select im2p file to load in...");
[fileName1,filePath1] = uigetfile; % get user input of .im2p file
cd(filePath1); %changes to the folder where the file was selected from
disp(strcat("Loading ",fileName1,"...")); tic;  % prints loading mesage
load(fileName1); toc;                           % loads selected file
gfpStack = uint16(a.gfp3FrameAvg);          % stores gfp frames from file
gfpStack_mean = uint16(mean(gfpStack,3));   % determine avg from all frames


%% Load RFP field image
disp("Select rfp .tif stack to load in...");
[fileName2,filePath2] = uigetfile('*.tif');
disp(strcat("Loading ",fileName2,"...")); tic;
[~,rfp_stack] = loadImg(fileName2); toc;
a.rfpPVImgPre = imresize(squeeze(mean(rfp_stack,3)),0.5);
rfpImg = a.rfpPVImgPre;


%% Run CDeisters Image Processing Code

importer;           

% roiMaker;

% extractor;

cimgBad = 0;

%% IF CIMG LOOKS BAD TRY:
cimgBad = 1;
figure; imagesc(gfpStack_mean); title('GFP Mean');
%%figure; imagesc(cimg); title('XCorr GFP Img');
gfpStack_median = median(gfpStack,3);
figure; imagesc(gfpStack_median); title('GFP Median');

if cimgBad == 1
    figure; imagesc(squeeze(mean(gfpStack(:,:,1:1000),3)));
    figure; imagesc(squeeze(mean(gfpStack(:,:,end-999:end),3)));
    gfpStack_mean = squeeze(mean(gfpStack(:,:,1:1000),3));
end

%% save image-related variables to im2p
a.gfpMotionCorrected = gfpStack;
if exist('meanProj_gfpStack')
    a.templateImg4MotionCorrection = meanProj_gfpStack;
else
    a.templateImg4MotionCorrection = gfpStack_mean;
end
if exist('rfpStack')
    rfpFlag = 1;
end
a.gfpXCorrImg = cimg;
a.somaticF = somaticF;
a.redSomaticF = redSomaticF;
a.somaticFDT = zeros(size(somaticF,1),size(somaticF,2));
for n = 1 : size(somaticF,1)
    a.somaticFDT(n,:) = detrend(somaticF(n,:)) + mean(somaticF(n,:));
end
a.redSomaticFDT = zeros(size(redSomaticF,1),size(redSomaticF,2));
for n = 1 : size(redSomaticF,1)
    a.redSomaticFDT(n,:) = detrend(redSomaticF(n,:)) + mean(redSomaticF(n,:));
end
%Somatic IDs
a.somaticROI_PixelLists = somaticROI_PixelLists;
a.somaticROIBoundaries = somaticROIBoundaries;
a.somaticROICenters = somaticROICenters;
a.somaticRoiCounter = somaticRoiCounter;
a.somaticROIs = somaticROIs;
%Red Somatic IDs
a.redSomaticROI_PixelLists = redSomaticROI_PixelLists;
a.redSomaticROIBoundaries = redSomaticROIBoundaries;
a.redSomaticROICenters = redSomaticROICenters;
a.redSomaticRoiCounter = redSomaticRoiCounter;
a.redSomaticROIs = redSomaticROIs;
%combined red / green image
if cimgBad == 1
    a.overlayImg = cat(3,rescale(rfpImg),rescale(gfpStack_mean),zeros(256,256));
else
    a.overlayImg = cat(3,rescale(rfpImg),rescale(cimg),zeros(256,256));
end
disp('Finished with img components!');

%% get experiment metadata
currDir = pwd;
disp("Select experiment metadata file to load in...");
[fileNameExp,filePathExp] = uigetfile;
cd(filePathExp);
expMeta = load(fileNameExp);
cd(currDir);

% assign metadata to .im2p
a.stimDurs = expMeta.Data.stimDurations;
a.stimFreqs = expMeta.Data.stimFrequencies;
a.stimAmps = expMeta.Data.stimAmplitudes;
a.stimOnTimes = expMeta.Data.stimONTimes;
a.stimOnFrameNum = expMeta.Data.stimOnFrameNum;
a.ledTimes = expMeta.Data.ledTimes;
a.audioTrials = expMeta.Data.audioTrials;
a.lickTimes = expMeta.Data.lickTimes;
a.imgDepth = expMeta.Data.imagingDepth;

%clean up & save grasping data
numGraspON = max(size(expMeta.Data.graspTimesRight));
numGraspOFF = max(size(expMeta.Data.graspTimesRightOff));
graspTimes = zeros(max(numGraspON,numGraspOFF),2); %horzcat graspON & graspOFF
graspTimes(1:numGraspON,1) = expMeta.Data.graspTimesRight;
graspTimes(1:numGraspOFF,2) = expMeta.Data.graspTimesRightOff;
graspDurs = zeros(min(numGraspON,numGraspOFF),1);
graspDurs = graspTimes(1:min(numGraspON,numGraspOFF),2) - graspTimes(1:min(numGraspON,numGraspOFF),1);
a.graspTimes = graspTimes;
a.graspDurs = graspDurs;

disp('Finished assigning metadata.');

%% remove bulky files from -analyzed im2p
a.gfp3FrameAvg = NaN;
a.gfpMotionCorrected = NaN;
a.rfpExpTrace = NaN;
disp('removed bulky files!');

% save updated im2p file
fileNameOut = [fileName1(1:end-4),'-analyzed.mat']
save(fileNameOut,'a','-v7.3'); 

% export overlayed photo
imwrite(a.overlayImg,strcat('pseudo-color-field-z',num2str(a.imgDepth),'.png')); 

%% grab the vasculature images
photoExtractor;

%% Finish Up
clear all; close all; clc; 
