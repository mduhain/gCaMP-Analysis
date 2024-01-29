% CLEAN PHOTO EXAMPLES FROM MICE FOR MANUSCRIPT
% 2024/01/22

%% Load im2p file
currDir = pwd;
disp("Select im2p file to load in...");
[fileName1,filePath1] = uigetfile; % get user input of .im2p file
cd(filePath1); %changes to the folder where the file was selected from
disp(strcat("Loading ",fileName1,"...")); tic;  % prints loading mesage
load(fileName1); toc;                           % loads selected file
rfpImg = a.rfpPVImgPre;
overlayImg = a.overlayImg;
gfpImg = a.templateImg4MotionCorrection;
cimg = a.gfpXCorrImg;

roiMaker


%% 

% CLEAN EXAMPLE IMAGES
dir = strcat(currDir,'\savedImages');
cd(dir);

%Loading
CIMG = rescale(mean(imread("expImg_1.png"),3));
GFP = rescale(mean(imread("expImg_2.png"),3));
RFP = rescale(mean(imread("expImg_3.png"),3));

%% CIMG 
% Comparison
figure('Color',[1,1,1])
subplot(1,2,1); imagesc(CIMG); axis square
subplot(1,2,2); imagesc(CIMG_2); axis square
% Alone
figure('Color',[1,1,1])
imagesc(CIMG); axis square
ax1 = gca; ax1.XTick = []; ax1.YTick = [];
fig1 = gcf; 
fig1.Colormap = cat(2,linspace(0,1,256)',linspace(0,1,256)',linspace(0,1,256)');


%% GFP
% Comparison
figure('Color',[1,1,1])
subplot(1,2,1); imagesc(GFP); axis square
subplot(1,2,2); imagesc(GFP); axis square

% Alone
figure('Color',[1,1,1])
imagesc(GFP); axis square
ax1 = gca; ax1.XTick = []; ax1.YTick = [];
fig2 = gcf; 
fig2.Colormap = cat(2,zeros(256,1),linspace(0,1,256)',zeros(256,1));


%% RFP
% Comparison
figure('Color',[1,1,1])
subplot(1,2,1); imagesc(RFP);
title("RFP");
subplot(1,2,2); imagesc(GFP);
title("GFP");

% Alone
figure('Color',[1,1,1])
imagesc(RFP); axis square
ax1 = gca; ax1.XTick = []; ax1.YTick = [];
fig3 = gcf; 
fig3.Colormap = cat(2,linspace(0,1,256)',zeros(256,1),linspace(0,0,256)');


%% COMPLETE
figure('Color',[1,1,1])
% GFP
subplot(1,3,1);
imagesc(GFP); axis square
ax1 = gca; ax1.XTick = []; ax1.YTick = [];
fig4 = gcf; 
fig4.Colormap = cat(2,zeros(256,1),linspace(0,1,256)',zeros(256,1));
title("GFP");
% RFP
subplot(1,3,2);
imagesc(RFP); axis square
ax1 = gca; ax1.XTick = []; ax1.YTick = [];
fig5 = gcf; 
fig5.Colormap = cat(2,linspace(0,1,256)',zeros(256,1),linspace(0,0,256)');
title("RFP");
% CIMG
subplot(1,3,3);
imagesc(CIMG); axis square
ax1 = gca; ax1.XTick = []; ax1.YTick = [];
fig6 = gcf; 
fig6.Colormap = cat(2,linspace(0,1,256)',linspace(0,1,256)',linspace(0,1,256)');
title("X-Corr IMG");


%% COMPOSITE IMGS

% CIMG and RFP
figure('Color',[1,1,1])
composite1 = cat(3,rescale(RFP),rescale(CIMG),zeros(size(GFP,1),size(GFP,2)));
imshow(composite1);
imwrite(composite1,'clean-XCORR-RFP.png');

% GFP and RFP
figure('Color',[1,1,1])
composite2 = cat(3,rescale(RFP),rescale(GFP),zeros(size(GFP,1),size(GFP,2)));
imshow(composite2);
imwrite(composite2,'clean-GFP-RFP.png');

% ALL 
figure('Color',[1,1,1])
composite3 = cat(3,rescale(RFP),rescale(GFP),rescale(CIMG));
imshow(composite3);
imwrite(composite3,'clean-RFP-GFP-XCorr.png');

% GFP and XCORR
figure('Color',[1,1,1])
composite4 = cat(3,zeros(size(GFP,1),size(GFP,2)),rescale(GFP),rescale(CIMG));
imshow(composite4);
imwrite(composite4,'clean-GFP-XCorr.png');












