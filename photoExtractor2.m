% photoExtractor2
%
% The portions of photoextractor that are used in widefield mapping
%
%

% Load in depth frames
[~,fileNames] = analyzeDir();
names = fileNames(contains(fileNames,'Z'));
allRFPimgs = [];
allGFPimgs = [];
for n = 1 : length(names)
    [gfpStack,rfpStack] = loadImg(names(n));
    rfpImg = mean(rfpStack,3);
    gfpImg = mean(gfpStack,3);
    if n == 1
        allRFPimgs = rfpImg;
        allGFPimgs = gfpImg;
    else
        allRFPimgs = cat(3,allRFPimgs,rfpImg);
        allGFPimgs = cat(3,allGFPimgs,gfpImg);
    end
    fileNameHere = char(names(n));
    fileNameHere(end-3:end) = [];
    %imwrite(rescale(rfpImg),[fileNameHere '.png']);
    figure('Color',[1 1 1]); imagesc(rfpImg);
    title(strcat("RFP ",fileNameHere));
    set(gca,'DataAspectRatio',[1 1 1])
    colormap gray
    figure('Color',[1 1 1]); imagesc(gfpImg);
    title(strcat("GFP ",fileNameHere));
    set(gca,'DataAspectRatio',[1 1 1])
    colormap gray
end

%Creat avg RFP image of all depths
wdfRFP = sum(allRFPimgs,3);
figure('Color',[1 1 1]); imagesc(wdfRFP);
set(gca,'DataAspectRatio',[1 1 1]);
colormap gray
imwrite(rescale(rfpImg),'widefieldAvgRFP.png');

%Creat avg GFP image of all depths
wdfGFP = sum(allGFPimgs,3);
figure('Color',[1 1 1]); imagesc(wdfGFP);
set(gca,'DataAspectRatio',[1 1 1]);
colormap gray
imwrite(rescale(gfpImg),'widefieldAvgGFP.png');

%% Image processing RFP
%visualize red image
figure('Color',[1 1 1]); imagesc(wdfRFP);
title("Draw polygon outside window boundary for background subtraction");
%Trace ROI of background RFP signal from outside the window
h = drawpolygon('FaceAlpha',0);
%Create mask
mask0 = poly2mask(h.Position(:,1),h.Position(:,2),512,512);
mask0 = ~mask0;
%Image cropping
wdfRFPcrop = wdfRFP;
wdfRFPcrop(mask0) = NaN;
rfpBackground = mean(wdfRFPcrop,'all','omitnan');
imRFP = wdfRFP - rfpBackground;
imRFP(imRFP < 0) = 0;
figure; imagesc(imRFP);
title("Background Corrected RFP Image");


%% Image processing GFP

%visualize green image
figure('Color',[1 1 1]); imagesc(wdfGFP);
title("Draw polygon inside usable section within window ");
%Trace ROI of usable GFP signal from inside window
h = drawpolygon('FaceAlpha',0);
%Create mask
mask = poly2mask(h.Position(:,1),h.Position(:,2),512,512);
mask = ~mask;
%Image cropping
wdfGFPcrop = wdfGFP;
wdfGFPcrop(mask) = NaN;
%subtract mean and remove negative values
imSM = wdfGFPcrop - mean(wdfGFPcrop,'all','omitnan');
imSM(imSM < 0) = 0;
% %Image denoise: using an empirical Bayesian method
% imDN = wdenoise2(wdfGFPcrop);
% %apply gauss filter
% imGF = imgaussfilt(wdfGFPcrop,1);

% figure('Color',[1 1 1]); imshow(cat(3,rescale(wdfRFP),rescale(imSM),rescale(z100resp)));

figure('Color',[1 1 1]); imshow(cat(3,rescale(wdfRFP),rescale(imSM),zeros(512,512)));


%% Calculate reponse map from Img Dir
currentDir = pwd;
cd(strcat(currentDir,'\Z110'));
widefieldAnalysis();
respMap = imread('response_map.png');
cd(currentDir);
respMap = double(imresize(respMap,4));

%% Load in Z depth dF/F maps
currentDir = pwd;
cd(strcat(currentDir,'\Z110'));
respMap1 = imread('response_map.png');
cd(currentDir);
cd(strcat(currentDir,'\Z73'));
respMap2 = imread('response_map.png');
cd(currentDir);
zResp = double(imresize(respMap1,4)) + double(imresize(respMap2,4));
figure('Color',[1 1 1]); imagesc(zResp);
% h2 = drawpolygon('FaceAlpha',0);
% %Create mask
% mask2 = poly2mask(h2.Position(:,1),h2.Position(:,2),512,512);
% mask2 = ~mask2;
%Image cropping
zRespCrop = zResp;
zRespCrop(mask) = NaN;

zRespSM = zRespCrop - mean(zRespCrop,'all','omitnan');
zRespSM(zRespSM < 0) = 0;

%replace NaN values with 0
imSM(isnan(imSM)) = 0;
zRespSM(isnan(zRespSM)) = 0;

%save out finished pseudo-colored image
imFull = cat(3,rescale(imRFP),rescale(imSM),rescale(zRespSM));
figure('Color',[1 1 1]); imshow(imFull);
imwrite(imFull,'completeMap.png');


