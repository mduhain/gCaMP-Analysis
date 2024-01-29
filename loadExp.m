%Load in GFP images and RFP trace from Experiment 
function [gfpStack, rfpTrace, scanInfo] = loadExp(fileName)

    %start loading process
    disp(strcat("<",string(datetime("now")),"> ","Loading imfinfo on large tiff..."));
    tic;
    imInfo = imfinfo(fileName);
    toc;
    
    %get stack dimenstions
    scanInfo.num_frames = size(imInfo,1); %interleaved (1.Green,2.Red,3.Green,4.Red...etc)
    scanInfo.Xpx = imInfo(1).Width;
    scanInfo.Ypx = imInfo(1).Height;
    %downsample images from 512x512 px to 256x256 px
    newX = floor(scanInfo.Xpx/2);
    newY = floor(scanInfo.Ypx/2);
    %downsample time 30fps to 10fps
    gfpStack = zeros(newX,newY,floor(scanInfo.num_frames/6));
    
    %populate GFP Matrix
    disp(strcat("<",string(datetime("now")),"> ","Extracting gfp frames..."));
    tic;
    fprintf("Loading frames ");
    for n = 1 : floor(scanInfo.num_frames/6)
        img1 = imresize(double(imread(fileName,6*n-5,'Info',imInfo)),0.5);
        img2 = imresize(double(imread(fileName,6*n-3,'Info',imInfo)),0.5);
        img3 = imresize(double(imread(fileName,6*n-1,'Info',imInfo)),0.5);
        gfpStack(:,:,n) = mean(cat(3,img1,img2,img3),3);
        if rem(n,floor(floor(scanInfo.num_frames/6)/10)) == 0
            fprintf(". ");
        end
    end
    fprintf("\n");
    toc;
    
    %prepopulate array to store RFP frame avg trace
    rfpTrace = zeros(floor(scanInfo.num_frames/6),1);

    %populate rfp trace
    disp(strcat("<",string(datetime("now")),"> ","Extracting RFP trace..."));
    tic;
    fprintf("Loading frames ");
    for n = 1 : floor(scanInfo.num_frames/6)
        img1 = double(imread(fileName,6*n-4,'Info',imInfo));
        img2 = double(imread(fileName,6*n-2,'Info',imInfo));
        img3 = double(imread(fileName,6*n,'Info',imInfo));
        rfpTrace(n) = mean(cat(1,img1,img2,img3),'all');
        if rem(n,floor(floor(scanInfo.num_frames/6)/10)) == 0
            fprintf(". ");
        end
    end
    fprintf("\n");
    toc;

end