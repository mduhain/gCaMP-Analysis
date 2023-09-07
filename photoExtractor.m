
function [] = photoExtractor()
    % photoExtractor 
    % - Identifies "Vasculature_Nx_00001.tif" files and converts them to PNG
    % - Auto adjusts pixel values of PNG output to normalized range from .tif (uint16)

    startingDir = pwd;
    myDir = dir;
    %remove '.' '..' directories
    if strcmp(myDir(1).name,'.')
        myDir(1) = [];
    end
    if strcmp(myDir(1).name,'..')
        myDir(1) = [];
    end

    %% extract vasculature photos
    [vasc1GfpStack,vasc1RfpStack] = loadImg('Vasculature_1x_00001.tif');
    vasc1 = mean(vasc1GfpStack,3);
    imwrite(rescale(vasc1),'vasc1.png');
    figure; imagesc(vasc1);
    set(gca,'DataAspectRatio',[1 1 1])
    colormap gray

    [vasc2GfpStack,vasc2RfpStack] = loadImg('Vasculature_2x_00001.tif');
    vasc2 = mean(vasc2GfpStack,3);
    imwrite(rescale(vasc2),'vasc2.png');
    figure; imagesc(vasc2);
    set(gca,'DataAspectRatio',[1 1 1])
    colormap gray

    [vasc3GfpStack,vasc3RfpStack] = loadImg('Vasculature_3x_00001.tif');
    vasc3 = mean(vasc3GfpStack,3);
    imwrite(rescale(vasc3),'vasc3.png');
    figure; imagesc(vasc3);
    set(gca,'DataAspectRatio',[1 1 1])
    colormap gray
    
    %% 

    [~,fileNames] = analyzeDir();
    redNames = fileNames(contains(fileNames,'RFP'));
    for n = 1 : length(redNames)
        [~,rfpStack] = loadImg(redNames(n));
        rfpImg = mean(rfpStack,3);
        fileNameHere = char(redNames(n));
        fileNameHere(end-3:end) = [];
        imwrite(rescale(rfpImg),[fileNameHere '.png']);
        figure; imagesc(rfpImg);
        set(gca,'DataAspectRatio',[1 1 1])
        colormap gray
    end
    
    %% 

    [~,fileNames] = analyzeDir();
    names = fileNames(contains(fileNames,'Z'));
    allRFPimgs = [];
    for n = 1 : length(names)
        [gfpStack,rfpStack] = loadImg(names(n));
        rfpImg = mean(rfpStack,3);
        if n == 1
            allRFPimgs = rfpImg;
        else
            allRFPimgs = cat(3,allRFPimgs,rfpImg);
        end
        gfpImg = mean(gfpStack,3);
        fileNameHere = char(names(n));
        fileNameHere(end-3:end) = [];
        %imwrite(rescale(rfpImg),[fileNameHere '.png']);
        figure; imagesc(rfpImg);
        title(strcat("RFP ",fileNameHere));
        set(gca,'DataAspectRatio',[1 1 1])
        colormap gray
        figure; imagesc(gfpImg);
        title(strcat("GFP ",fileNameHere));
        set(gca,'DataAspectRatio',[1 1 1])
        colormap gray
    end
    widefieldAvgImg = sum(allRFPimgs,3);
    figure; imagesc(widefieldAvgImg);
    set(gca,'DataAspectRatio',[1 1 1]);
    colormap gray
    imwrite(rescale(rfpImg),'widefieldAvg.png');
    
    


end

