% im2pAnalysis
% - A semi-automatic script to analyze im1p data from Tactile Stim
% Experiments in haptics lab mice.
%
% update: 2022-08-29
% First script to work with im2p files.
%
% update: 2022-09-16
% Added output to PDF functionality
%
%
%
%
%
%



%% Load in im2p file
a = im2p;
makeReport = 0;
disp('Please select the experiment folder to analyze...');
fileDir = uigetdir;
cd(fileDir);

%%
[~, fileNames] = analyzeDir();

%%
fprintf('Loading im2p file...');
load(fileNames(contains(fileNames,'im2p')));
if exist('output')
    a = output;
else

end
a = 
clear output;
toc;

gfpStack = uint16(a.gfp3FrameAvg);
gfpStack_mean = uint16(mean(gfpStack,3));
importer

%% Create plots and figs

%View all RFP images
figure;
imagesc(a.rfpPVImgPre); hold on;
title("RFP Img Stack Avg"); hold off;

%View composite image from red + green channels
gfp10 = mean(a.gfp3FrameAvg(:,:,1:floor(size(a.gfp3FrameAvg,3)/10)),3);
a.imDual(gfp10,imresize(a.rfpPVImgPre,size(gfp10,1)/size(a.rfpPVImgPre,1)))
title("Pseudo-color Field Image, PV(rfp) gCaMP(gfp)");


roiMaker

if size(a.rfpPVImgPre,1) == 512
    rfpImg = imresize(a.rfpPVImgPre,0.5);
end

%replot GFP information
if exist('cimg_gfpStack','var') == 1
    h2 = figure; hold on;
    imagesc(cimg_gfpStack);
else
    h2 = a.imStack(a.gfp3FrameAvg,1,floor(size(a.gfp3FrameAvg,3)/10));
end
title("GFP ROIs");



for n = 1 : size(somaticF,1)
    plot(h2.CurrentAxes,somaticROICenters{1,n}.Centroid(1),somaticROICenters{1,n}.Centroid(2),'ko')
    text(somaticROICenters{1,n}.Centroid(1)+3,somaticROICenters{1,n}.Centroid(2)-5,strcat('#',num2str(n)));
end

a.movie(gfpStack,10)

%make PSTH from ROI
a.psth2(a.somaticF,a.stimOnFrameNum,a.stimFreqs);


%% Save analysis to new variable
a.gfpMotionCorrected = gfpStack;
if exist('meanProj_gfpStack')
    a.templateImg4MotionCorrection = meanProj_gfpStack;
else
    a.templateImg4MotionCorrection = gfpStack_mean;
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
%Check if stim frame times is empty (if so, manually calculate)
if sum(a.stimOnFrameNum ~= 0) == 0
    smallRFPtrace = a.rfpExpTrace;
    smallRFPtrace = smallRFPtrace + 1; 
    figure;
    plot(smallRFPtrace);
    hold on;
    threshVal = 2;
    plot([1 size(smallRFPtrace,1)],[threshVal*mean(smallRFPtrace) threshVal*mean(smallRFPtrace)],'k-');
    I = find(smallRFPtrace > threshVal*mean(smallRFPtrace));
    dI = diff(I);
    dI = [0; dI];
    I(dI == 1) = [];
    I(1:3) = []; %remove first three flashes (exp start id)
    a.stimOnFrameNum = I;
end
%Remove bulky raw image files
a.gfpMotionCorrected = NaN;
a.gfp3FrameAvg = NaN;
a.rfpExpTrace = NaN;
%file output
newFileName = strsplit(fileNames(contains(fileNames,'im2p')),'.');
newFileName = strcat(newFileName(1),'-Analyzed.',newFileName(2))
save(newFileName,'a','-v7.3');

%%Remove NaN's


%% LOCAL FUNCTIONS

%% ------------------------------------------------------------------------
function [img, images] = addToPDF(images, fig, titleA)
    import mlreportgen.dom.*
    % Set figure size, recommended
    values = [3.5, 2.5];
    fig.PaperSize = values;
    fig.PaperPosition = [0 0 values];
    fig.Units = 'inches';
    fig.Position(3:4) = values;

    % Add the plot to the document
    name = sprintf('%s.svg', titleA);
    print(fig, name, '-dsvg');
    img = Image(name);
    delete(fig) %delete plot figure window
    images = [images {img}];
end

%% ------------------------------------------------------------------------
function [folderNames, fileNames] = analyzeDir()
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
