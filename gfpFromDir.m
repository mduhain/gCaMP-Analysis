function [gfpStack] = gfpFromDir()
    inputFiles = dir(fullfile('*.tif')); %find all .tif files in dir
    fileNames = {inputFiles.name};
    imInfo = imfinfo(fileNames{1}); %get info on first tif
    Xpx = imInfo(1).Width; %X pixels
    Ypx = imInfo(1).Height; %Y pixels
    numFrames = length(inputFiles); %Num of available frames
    if Xpx == 512 && Ypx == 512
        disp('GFP data is too large... Down-sampling...');
        gfpStack = zeros(Xpx/2,Ypx/2,floor(numFrames/3));
        tic;
        parfor k = 1 : floor(numFrames/3)
%             if rem(k,1000) == 0
%                 fprintf(strcat(num2str(k),'.'));
%             end
            gfp1 = imresize(imread(fileNames{3*k-2},1),0.5);
            gfp2 = imresize(imread(fileNames{3*k-1},1),0.5);
            gfp3 = imresize(imread(fileNames{3*k},1),0.5);
            gfpStack(:,:,k) = squeeze(mean(cat(3,gfp1,gfp2,gfp3),3));
        end
        toc;
    end
end