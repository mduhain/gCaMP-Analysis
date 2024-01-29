% Large scale forepaw mapping
%
%

%forepaw image path
startingPath = "C:\Users\ramirezlab\Documents\Michael\ForepawMap";
cd(startingPath);

[~,mapNames] = analyzeDir();


img = imread(mapNames(2));
figure; subplot(1,2,1); imagesc(img); hold on;
colormap gray

startDraw();






cd(strcat(startingPath,'\',mouseList(3)));
img = imread('Full_window_m810_rotated.tif');
img810 = img(:,:,2);
figure; imagesc(img810); hold on;
colormap gray
title(mouseList(3));

cd(strcat(startingPath,'\',mouseList(1)));
img = imread('Full_window_m792_rotated.tif');
img792 = img(:,:,2);
figure; subplot(1,2,1); imagesc(img792); hold on;
colormap gray
title(mouseList(1)); 
hold off;







function [out] = startDraw()
    [filename,pathname] = uigetfile;
    cd(pathname)
    data = imread(filename);
    title("Allignment");
    subplot(1,2,2)
    imagesc(data);
    input('Finished...? hit [Enter]');
    


end
