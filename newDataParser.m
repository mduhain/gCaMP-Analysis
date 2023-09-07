% NEW DATAPARSER

testAllFiles = 1; %FLAG: Load each im2p file and check for accuracy;
if ismac() == 1
    %startingDir = 'C:\Users\ramirezlab\Desktop\2P-Data';
else
    startingDir = 'C:\dataParser\2P-Data';
    expMetaDir = 'C:\dataParser\Exp-Metadata';
end

cd(startingDir);
[mouseNum, ~] = analyzeDir();

for n = 1 : length(mouseNum) %mouse number loop
    mouseDir = strcat(startingDir,'\',mouseNum(n));
    cd(mouseDir);
    [dates, ~] = analyzeDir();
    for nd = 1 : length(dates) %date loop
        cd(strcat(mouseDir,'\',dates(nd)));
        [expNames,~] = analyzeDir();
        for ne = 1 : length(expNames) %experiment loop
            cd(strcat(mouseDir,'\',dates(nd),'\',expNames(ne)));
            [depths, fileNames] = analyzeDir();
            if any(contains(fileNames,'im2p'))
                disp(strcat("IM2P found, skipping ",mouseDir,'\',dates(nd),'\',expNames(ne)));
                continue
            end     
            %analyze
            for k = 1 : length(depths) %loop through different depths
                %change to Z000 directory (depth)
                cd(strcat(mouseDir,'\',dates(nd),'\',expNames(ne),'\',depths(k)));
                
                %extract gfp signals, resize 512->256 && 3frame Avg
                gfpStack = gfpFromDir();
                
                %extract rfp signal (trace)
                rfpTrace = rfpFromAllFiles();
                
                %save parameters into im2p
                a = im2p;
                
                % back out 1 directory to 'Exp000'
                cd(strcat(mouseDir,'\',dates(nd),'\',expNames(ne)));
                
                % add other params. rfp imgs, vasculature 
                a.gfp3FrameAvg = gfpStack;
                a.rfpExpTrace = rfpTrace;
                
                % save updated im2p
                fileOutName = strcat(mouseNum(n),'-',dates(nd),'-',expNames(ne),'-',depths(k),'-im2p.mat');
                save(fileOutName,'a','-v7.3');
                
                %clear workspace
                clear a gfpStack rfpTrace
                
            end % end of depth loop
            
            
        end
        
    end
    
end


%--------------------------------------------------------------------------
%extrace the red frames to make a rfp trace for all files in a directory
function [rfpTrace] = rfpFromAllFiles()
    inputFiles = dir(fullfile('*.tif')); %find all .tif files in dir
    fileNames = {inputFiles.name};
    imInfo = imfinfo(fileNames{1}); %get info on first tif
    Xpx = imInfo(1).Width; %X pixels
    Ypx = imInfo(1).Height; %Y pixels
    numFrames = length(inputFiles); %Num of available frames
    rfpTrace = zeros(numFrames,1);
    disp("extracting RFP frames...");
    tic;
    for k = 1 : numFrames
        try
            thisFileName = fileNames{k};
            rfpImg = imread(thisFileName,2);
            rfpTrace(k) = mean(rfpImg,'all');
        catch
            disp(strcat("Error on frame # ",num2str(k)));
        end
    end
    toc;
end
%--------------------------------------------------------------------------




