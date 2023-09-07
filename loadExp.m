function [gfpStack, rfpTrace, scanInfo] = loadExp(fileName,opt)
    % function loadExp.m
    % [gfpStack, rfpTrace, scanInfo] = loadExp(fileName,opt)
    %
    % from mduhain (2023-02-27)
    % loads in, and resizes, a large .tif file (gfp/rfp) into MAT array
    % 
    % opt == 1 (extracts both rfp trace & gfpStack)
    %
    
    imResized = 0;
    imInfo = imfinfo(fileName);
    scanInfo.num_frames = size(imInfo,1);
    scanInfo.Xpx = imInfo(1).Width;
    scanInfo.Ypx = imInfo(1).Height;
    if opt == 1 
        %prepopulate array to store GFP Frames
        try
            gfpStack = zeros(scanInfo.Xpx,scanInfo.Ypx,floor(scanInfo.num_frames/6));
        catch ME
            if strcat(ME.identifier,'MATLAB:array:SizeLimitExceeded')
                fprintf("Resizing...");
                imResized = 1;
                newX = floor(scanInfo.Xpx/2);
                newY = floor(scanInfo.Ypx/2);
                gfpStack = zeros(newX,newY,floor(scanInfo.num_frames/6));
            end
        end
        
        %prepopulate array to store RFP frames
        rfpTrace = zeros(floor(scanInfo.num_frames/2),1);
        
        %populate GFP Matrix
        for n = 1 : floor(scanInfo.num_frames/6)
            if imResized == 1
                img1 = imresize(double(imread(fileName,6*n-5,'Info',imInfo)),0.5);
                img2 = imresize(double(imread(fileName,6*n-3,'Info',imInfo)),0.5);
                img3 = imresize(double(imread(fileName,6*n-1,'Info',imInfo)),0.5);
                gfpStack(:,:,n) = mean(cat(3,img1,img2,img3),3);
            else
                img1 = imread(fileName,6*n-5,'Info',imInfo);
                img2 = imread(fileName,6*n-3,'Info',imInfo);
                img3 = imread(fileName,6*n-1,'Info',imInfo);
                gfpStack(:,:,n) = mean(cat(3,img1,img2,img3),3);
            end
        end
        
        for n = 1 : floor(scanInfo.num_frames/2)
           rfpTrace(n) = mean(mean(imread(fileName,2*n,'Info',imInfo))); 
        end
    elseif opt == 0
        gfpStack = zeros(scanInfo.Xpx,scanInfo.Ypx,floor(scanInfo.num_frames/3)+1);
        rfpTrace = [];
        for n = 1 : floor(scanInfo.num_frames/3)
            img1 = imread(fileName,3*n-2,'Info',imInfo);
            img2 = imread(fileName,3*n-1,'Info',imInfo);
            img3 = imread(fileName,3*n,'Info',imInfo);
            gfpStack(:,:,n) = mean(cat(3,img1,img2,img3),3);
        end
    else 
        error("var opt not recognized! must be 1 or 0");
    end
    clear img1 img2 img3 
end