
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
%     [vasc1GfpStack,vasc1RfpStack] = loadImg('Vasculature_1x_00001.tif');
%     vasc1 = mean(vasc1GfpStack,3);
%     imwrite(rescale(vasc1),'vasc1.png');
%     figure; imagesc(vasc1);
%     set(gca,'DataAspectRatio',[1 1 1])
%     colormap gray
% 
%     [vasc2GfpStack,vasc2RfpStack] = loadImg('Vasculature_2x_00001.tif');
%     vasc2 = mean(vasc2GfpStack,3);
%     imwrite(rescale(vasc2),'vasc2.png');
%     figure; imagesc(vasc2);
%     set(gca,'DataAspectRatio',[1 1 1])
%     colormap gray
% 
%     [vasc3GfpStack,vasc3RfpStack] = loadImg('Vasculature_3x_00001.tif');
%     vasc3 = mean(vasc3GfpStack,3);
%     imwrite(rescale(vasc3),'vasc3.png');
%     figure; imagesc(vasc3);
%     set(gca,'DataAspectRatio',[1 1 1])
%     colormap gray
    

    %%
    
    [~,fileNames] = analyzeDir();
    vascNames = fileNames(contains(fileNames,'Vasc'));
    for n = 1 : length(vascNames)
        [~,vascStack] = loadImg(vascNames(n));
        vascImg = mean(vascStack,3);
        fileNameHere = char(vascNames(n));
        fileNameHere(end-3:end) = [];
        imwrite(rescale(vascImg),[fileNameHere '.png']);
        figure; imagesc(vascImg);
        set(gca,'DataAspectRatio',[1 1 1])
        colormap gray
    end



    %% 

    [~,fileNames] = analyzeDir();
    redNames = fileNames(contains(fileNames,'RFP','IgnoreCase',true));
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

end

