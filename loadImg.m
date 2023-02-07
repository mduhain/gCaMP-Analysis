function [gfp_stack,rfp_stack] = loadImg(fileName)
%Loads in a .tif/.tiff image from the 2P Microscope
% - Alternates Green / Red frames
% - Returns 2 image stacks, one for each channel
    imInfo = imfinfo(fileName);
    num_frames = size(imInfo,1);
    Xpx = imInfo(1).Width;
    Ypx = imInfo(1).Height;
    rfp_stack = zeros(Xpx,Ypx,floor(num_frames/2));
    gfp_stack = zeros(Xpx,Ypx,floor(num_frames/2));
    for n = 1 : floor(num_frames/2)
        gfp_stack(:,:,n) = imread(fileName,2*n-1,'Info',imInfo);
        rfp_stack(:,:,n) = imread(fileName,2*n,'Info',imInfo);
    end
end