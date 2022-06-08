function [outputArg1,outputArg2] = ImNeuron(file,path)
%ImNeuron: Identify, mask, and track neuron activity in a tif stack.
%   Accepts file name input ('filename.tif') in working dir
%   Also accepts two inputs of a filename & directory of file
%   Assumes .tif stacks are [GFP,RFP,...,GFP,RFP]

%   0. Load in file / prepare workspace.
%   1. Identify Neuron locations (GFP / RFP)
%   2. Track frame drift through stack
%   3. Build masks (square array 1.5x size of neuron radius, pad nans)
%   4. Apply masks through stack

% Creation June 07 2022 - @mduhain
%
%

% CHECK INPUTS AND PREP WORKSPACE
if nargin == 0
    %input blank, just run script
    if ismac
        cd('/Users/mduhain/Desktop/PROJECT/2P-Ca-day1');
        file = "Co-labelled neuron 3 stim 31122.tif";
        im_info = imfinfo(file);
        num_frames = length(im_info)/2;
        im_width = im_info(1).Width;
        im_height = im_info(1).Height;
    else
        %labcomputer
    end
elseif nargin == 1
    if exist(file,'file')
        im_info = imfinfo(file);
        num_frames = length(im_info)/2;
        im_width = im_info(1).Width;
        im_height = im_info(1).Height;
    else
        error('File not found in current directory.');
    end
elseif nargin ==2
    %change dir
    wr = pwd;
    cd(path);
    %get image info
    im_info = imfinfo(file);
    num_frames = length(im_info)/2; %for each channel
    im_width = im_info(1).Width;
    im_height = im_info(1).Height;
    %return to previous dir
    cd(wd);
end

if rem(num_frames,2)
    error('Your .tif stack has an odd number of frames');
end

%LOAD IN FILE
tic;
im_gfp = zeros(im_height,im_width,num_frames);
im_rfp = zeros(im_height,im_width,num_frames);
fprintf(strcat("Loading ",num2str(num_frames)," frames"));
dift = zeros(im_height-1,im_width,5);
for n = 1:num_frames
    im_gfp(:,:,n) = mat2gray(imread(file,n*2-1));
    im_rfp(:,:,n) = mat2gray(imread(file,n*2));
    if rem(n,num_frames/10) == 0
        fprintf(".");
    end
end
fprintf(strcat("in ",num2str(toc)," sec \n"));


%%
%im_diff = diff(rfp_stack(:,:,1));
%im_diff(im_diff ~= 0) = 1;
%dift(:,:,n) = im_diff;

%stuff for gif making
figure;
h1 = gcf;
h1.Position = [100 100 800 300];
filename = 'new_cell_tracking.gif';

step_size = 10; %frames to average across
for ns=1:step_size:num_frames
    im_overlay = mean(im_rfp(:,:,ns:ns+(step_size-1)),3);
    im_avg = mean(im_overlay,'all');
    im_std = std2(im_overlay);
    im_overlay(im_overlay < (im_avg + 4*im_std)) = 0;
    %REMOVE ISOLATED PEAKS
    for n=1:im_width %from columns
        x = im_overlay(:,n)'; %flip
        y = conv(x,[1,1,1],'same'); %find convolution
        ind = find(x==y); %index isolated peaks
        x(ind(2:end-1)) = 0; %remove isolated peaks (no edges!)
        im_overlay(:,n) = x'; %replace original vector
    end
    for n=1:im_height %from rows
        x = im_overlay(n,:);
        y = conv(x,[1,1,1],'same'); %find convolution
        ind = find(x==y); %index isolated peaks
        x(ind(2:end-1)) = 0; %remove isolated peaks (no edges!)
        im_overlay(n,:) = x; %replace original vector
    end
    for n=1:im_width %from columns
        x = im_overlay(:,n)'; %flip
        y = conv(x,[1,1,1],'same'); %find convolution
        ind = find(x==y); %index isolated peaks
        x(ind(2:end-1)) = 0; %remove isolated peaks (no edges!)
        im_overlay(:,n) = x'; %replace original vector
    end
    
    %check if multiple cells (pixel clusters) remain
    [~,L] = bwboundaries(im_overlay,'noholes');
    num_clust = max(max(L));
    if num_clust > 1
        bright = zeros(1,num_clust);
        %two or more clusters, find greatest brightness
        for n=1:max(max(L)) %loop through clusters
            mask = L; %duplicate array
            mask(mask~=n) = 0; %erase all other clusters
            im_mask = mask.*im_overlay; %mask signal
            im_mask(im_mask == 0) = []; %remove zeros
            bright(n) = sum(im_mask); %average brightness
        end
        [~,max_clust] = max(bright); %brightest cluster
        L(L~=max_clust) = 0; 
        L(L==max_clust) = 1; %L is now correct mask of brightest cluster
        subplot(1,2,1);
        imagesc(im_overlay.*L);
        pause(0.05);
        disp(ns);
    else
        subplot(1,2,1);
        imagesc(im_overlay)
        pause(0.05);
        disp(ns);
    end
    subplot(1,2,2);
    imagesc(mean(im_rfp(:,:,ns:ns+(step_size-1)),3));
          
    
    
end












%Prepare outputs
outputArg1 = NaN;
outputArg2 = NaN;
end

