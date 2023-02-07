function [rfpImage] = loadRfpImg(fileName,opt)
%opt defines the frame skip feature, interleaved green & red frames.
    imInfo = imfinfo(fileName);
    num_frames = size(imInfo,1);
    Xpx = imInfo(1).Width;
    Ypx = imInfo(1).Height;
    if opt == 1
        rfp_stack = zeros(Xpx,Ypx,floor(num_frames/2));
        for n = 1 : floor(num_frames/2)
            rfp_stack(:,:,n) = imread(fileName,2*n,'Info',imInfo);
        end
    elseif opt == 0
        rfp_stack = zeros(Xpx,Ypx,num_frames);
        for n = 1 : num_frames
            rfp_stack(:,:,n) = imread(fileName,n,'Info',imInfo);
        end
    else
        error("var opt not recognized! must be 1 or 0");
    end
    rfpImage = mean(rfp_stack,3);  
end