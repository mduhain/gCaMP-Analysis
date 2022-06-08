figure;
h1 = gcf;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated1.gif';
%file_in = "Co-labelled neurons 2 stim 31122.tif"; %OG
file_in = "exp_01.tif";
t_info = imfinfo(file_in);
num_frames = length(t_info)/2;
Ca_data_pv1 = zeros(1,num_frames);
Ca_data_pv2 = zeros(1,num_frames);
Ca_data_py1 = zeros(1,num_frames);
Ca_data_py2 = zeros(1,num_frames);

max_gfp = zeros(16384,1);

%PV With PV correlation
% PV with PY Correlation
% Ask jingyu for forepaw brush for verifying low mag response in forepaw


disp(strcat(num2str(length(t_info))," total frames."));
for n = 1:length(t_info)
    % Draw plot for y = x.^n
    %x = 0:0.01:1;
    %y = x.^n;
    %plot(x,y) 
    %drawnow 
      % Capture the plot as an image 
      if rem(n,2) %odd
          temp_im_green = mat2gray(imread(file_in,n));
          %max_gfp = (max_gfp(temp + temp_im_green) ./ 2;
          %max_gfp = temp_im_green(temp_im_green > max_gfp);
      else
          temp_im_red = mat2gray(imread(file_in,n));
          rgb_im = cat(3,temp_im_red,temp_im_green,zeros(128,128));
          
          %VIEW PER FRAME IMAGE
          subplot(2,2,1);
          imshow(cat(3,zeros(128,128),temp_im_green,zeros(128,128)));
          subplot(2,2,[2 4]);
          imshow(rgb_im);
          subplot(2,2,3);
          imshow(cat(3,temp_im_red,zeros(128,128),zeros(128,128)));
          h1.Position=[100 100 800 400];
          %input(num2str(n/2));
          
          if n==2 %|| rem(n,20) == 0 %on first, and every 10 frames)
              %find ROIs
%               [centers,radii] = imfindcircles(temp_im_red,[2 10],'ObjectPolarity','bright','Sensitivity',0.86);
%               %calculate image size
%               dims = size(temp_im_red);
%               imageSizeX = dims(1);
%               imageSizeY = dims(2);
%               [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
%               %ROI components
%               % main PV neuron (PV1)
%               centerX = centers(1,1);
%               centerY = centers(1,2);
%               radius = radii(1);
%               
%               %other PV Neuro n
%               PV2_cent_x = 93;
%               PV2_cent_y = 41;
%               PV2_rad = 4.0;
%               
%               %main GFP neuron N1
%               N1_cent_x = 22;
%               N1_cent_y = 31;
%               N1_rad = 4.2;
%               
%               %other GFP neuron N2
%               N2_cent_x = 34;
%               N2_cent_y = 105;
%               N2_rad = 4.2
              
              %create mask for isolating out calcium signal
%               mask1 = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= radius.^2; %PV1
%               mask2 = (rowsInImage - N1_cent_y).^2 + (columnsInImage - N1_cent_x).^2 <= N1_rad.^2; %N1
%               mask3 = (rowsInImage - N2_cent_y).^2 + (columnsInImage - N2_cent_x).^2 <= N2_rad.^2; %N2
%               mask4 = (rowsInImage - PV2_cent_y).^2 + (columnsInImage - PV2_cent_x).^2 <= PV2_rad.^2; %PV2
          end  
          
                
          %Calculating Calcium intensity for frame
%           %PV1
%           masked_im_1 = temp_im_green .* mask1;
%           masked_im_1(masked_im_1 == 0) = []; %clear away 0s
%           Ca_data_pv1(n/2) = mean(masked_im_1);
% 
%           %PV2
%           masked_im_4 = temp_im_green .* mask4;
%           masked_im_4(masked_im_4 == 0) = []; %clear away 0s
%           Ca_data_pv2(n/2) = mean(masked_im_4);
%           
%           %N1
%           masked_im_2 = temp_im_green .* mask2;
%           masked_im_2(masked_im_2 == 0) = []; %clear away 0s
%           Ca_data_py1(n/2) = mean(masked_im_2);
% %           
% %           %N2
%           masked_im_3 = temp_im_green .* mask3;
%           masked_im_3(masked_im_3 == 0) = []; %clear away 0s
%           Ca_data_py2(n/2) = mean(masked_im_3);
          
          %Big Counter
          disp(num2str(n));
          
          %visualizing changes

%           h = viscircles(centers,radii);
          frame = getframe(h1); 
          im = frame2im(frame); 
          [imind,cm] = rgb2ind(im,256); 
          % Write to the GIF File 
          if n == 2 
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
          else 
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1); 
          end 
      end
end

%


%% cross-correlations from MGR
%figure;
fs = 10;
step_size = 2; %sec
length_signal = num_frames;
bin_size = step_size * fs;
num_bins = round(length_signal/bin_size);
frame_corrs = zeros(num_bins,2*bin_size-1,3);
pv1_norm = zeros(num_bins,bin_size);
pv2_norm = zeros(num_bins,bin_size);
py1_norm = zeros(num_bins,bin_size);
py2_norm = zeros(num_bins,bin_size);
start_frame = 1;
end_frame = bin_size;

for n=1:num_bins
    
    py1_norm(n,:) = Ca_data_py1(start_frame:end_frame) - mean(Ca_data_py1(start_frame:start_frame+3));
    py2_norm(n,:) = Ca_data_py2(start_frame:end_frame) - mean(Ca_data_py2(start_frame:start_frame+3));
    pv1_norm(n,:) = Ca_data_pv1(start_frame:end_frame) - mean(Ca_data_pv1(start_frame:start_frame+3));
    pv2_norm(n,:) = Ca_data_pv2(start_frame:end_frame) - mean(Ca_data_pv2(start_frame:start_frame+3));

    %frame_corrs(start_frame,:) = xcorr(py1_norm,pv1_norm,'unbiased');
    %frame_corrs(start_frame,:) = xcorr(py1_norm,py2_norm,'unbiased');
    frame_corrs(n,:,1) = xcorr(py1_norm(n,:),py2_norm(n,:),'unbiased'); 
    frame_corrs(n,:,2) = xcorr(py1_norm(n,:),py2_norm(n,:),'none'); 
    frame_corrs(n,:,3) = xcorr(py1_norm(n,:),py2_norm(n,:),'norm'); 
    %try this with 'none' pr 'biased' and then the normalized function
    
    start_frame = start_frame + bin_size;
    end_frame = end_frame + bin_size;
    %plot(frame_corrs(n,:));
    %title(strcat("correlation frames ",num2str(n),"-",num2str(n+99)));
    %input("...");
end

num_shuffles = 5000;
shuff_corrs = zeros(num_shuffles,2*bin_size-1,3);
for n = 1 : num_shuffles
    rand_bins = randi(num_bins,1,2); % 2 random bins
    shuff_corrs(n,:,1) = xcorr(py1_norm(rand_bins(1),:),py2_norm(rand_bins(2),:),'unbiased');
    shuff_corrs(n,:,2) = xcorr(py1_norm(rand_bins(1),:),py2_norm(rand_bins(2),:),'none');
    shuff_corrs(n,:,3) = xcorr(py1_norm(rand_bins(1),:),py2_norm(rand_bins(2),:),'norm');
end



%add corrections, shuffle frames and make correlations, 
% signal 1:24 and signals 48:72

% compare random bins of (1:120), 5000 times, average and subtract from
% avereaged frame_corrs

  
%% Finding ROI's

%take the first set of images (R + G)
red_ch = mat2gray(imread(file_in,2));
green_ch = mat2gray(imread(file_in,1));
rgb_im = cat(3,red_ch,green_ch,zeros(128,128));
imshow(rgb_im);

%Get user input on one cell of interest
d=drawline;

%Calculate size
pos = d.Position;
diffPos = diff(pos);
diameter = hypot(diffPos(1),diffPos(2))

%Find Cells
[centers,radii] = imfindcircles(red_ch,[2 10],'ObjectPolarity','bright','Sensitivity',0.9)


%capturing calcium dynamics
val1 = [];
for n=1:2:1000
val1 = [val1 mean(mean(mat2gray(imread(file_in,n))))];
end
val2 = [];
for n=1001:2:2001
val2 = [val2 mean(mean(mat2gray(imread(file_in,n))))];
end
val3 = [];
for n=2001:2:3001
val3 = [val3 mean(mean(mat2gray(imread(file_in,n))))];
end
figure;
subplot(2,1,1);
plot(val1);
subplot(2,1,2);
plot(val2);




function mat = gauss2d(mat, sigma, center)
gsize = size(mat);
[R,C] = ndgrid(1:gsize(1), 1:gsize(2));
mat = gaussC(R,C, sigma, center);
end