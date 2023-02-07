%function [gfp_stack,rfp_stack] = quickAnalysis(fileName)
fileName = 'testFile_00001.tif';
%Loads in a .tif/.tiff image from the 2P Microscope
% - Alternates Green / Red frames
% - Returns 2 image stacks, one for each channel
    imInfo = imfinfo(fileName);
    frame_rate = 30; % frames per second
    time_window = 120; % in seconds, to be analyzed
    num_frames = frame_rate * time_window;
    gfp_fps = 10;

    Xpx = imInfo(1).Width;
    Ypx = imInfo(1).Height;

    rfp_img = zeros(Xpx,Ypx);
    gfp_img = zeros(Xpx,Ypx);

    rfp_trace = zeros(num_frames,1);
    gfp_trace = zeros(floor(num_frames/(frame_rate/gfp_fps)),4); %downsample 10hz, seprate 4 quadrants

    for n = 1 : num_frames
        rfp_img = imread(fileName,2*n,'Info',imInfo);
        rfp_trace(n) = mean(rfp_img,'all');
    end

    for n = 1 : floor(num_frames/3)
        gfp_img = imread(fileName,6*n-5,'Info',imInfo);
        gfp_trace(n,1) = mean(gfp_img(1:Xpx/2,1:Ypx/2),'all'); %upper left
        gfp_trace(n,2) = mean(gfp_img(1:Xpx/2,Ypx/2+1:end),'all'); %upper right
        gfp_trace(n,3) = mean(gfp_img(Xpx/2+1:end,1:Ypx/2),'all'); %lower left
        gfp_trace(n,4) = mean(gfp_img(Xpx/2+1:end,Ypx/2+1:end),'all'); %lower right
    end

    %plotting section
    figure;
    xRfp = [1/frame_rate:1/frame_rate:time_window];
    xGfp = [1/gfp_fps:1/gfp_fps:time_window];
    plot(xGfp,gfp_trace,'g');
    hold on;
    plot(xRfp,rfp_trace,'r');
    workingDir = pwd;
    title(workingDir);
    ylabel("Avg. Brightness / Frame");
    xlabel("Time (sec)");
    hold off;
%end