% dataParser.m
% Large script to go through all 2P data, aggregate .tiff stacks into .mat
% Containing:
%   1. Raw RFP Pre-Experiment signal (512 x 512 px image)
%   2. RFP trace from experiment, [1 x NumFrrames] array.
%   3. GFP data with 3 frame averaged images [xpx, ypx, numFrames/6]
%
% mduhain 2022-08-14
%

testAllFiles = 1; %FLAG: Load each im2p file and check for accuracy;
if ismac() == 1
    %startingDir = 'C:\Users\ramirezlab\Desktop\2P-Data';
else
    startingDir = 'C:\Users\skich\Documents\dataParser\2P-Data';
    expMetaDir = 'C:\Users\skich\Documents\dataParser\Exp-Metadata';
    %startingDir = 'V:\Haptics-2P\2P-Data';
    %expMetaDir = 'V:\Haptics-2P\Exp-Metadata';
    %startingDir = 'C:\Users\scanimage\Desktop\dataParser\2P-Data';
    %expMetaDir = 'C:\Users\scanimage\Desktop\dataParser\Exp-Metadata';
end

specificLoadPath = 0;
if specificLoadPath == 1
    mouseNum = "m806";
    dates = "2022-11-02";
    dataPath2P = strcat("D:\Haptics-2P\2P-Data\",mouseNum,"\",dates);
    dataPathMeta = strcat("D:\Haptics-2P\Exp-Metadata\",mouseNum,"\",dates);
else
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
            [~, fileNames] = analyzeDir();
            if sum(contains(fileNames,'im2p')) > 0
                disp(strcat(mouseNum(n),'-',dates(nd),'-',expNames(ne)," im2p found."));
                if testAllFiles == 0
                    disp('skip to next...');
                    continue;
                end

                %Load in im2p File
                fprintf("Loading im2p file..."); tic;
                load(fileNames(contains(fileNames,'im2p')));
                fprintf(strcat(num2str(toc)," sec.\n"))

                %Extract Metadata -----------------------------------------
                %extractMetadata(expMetaDir);

                %find and add experiment metadata
                cd(strcat(expMetaDir,'\',mouseNum(n),'\',dates(nd),'\'));
                [~, metaFileNames] = analyzeDir();
                cellData = load(metaFileNames(ne));
                cellData = cellData.Data;
                if iscell(cellData)
                    cleanStruct = decodeCell(cellData);
                else
                    cleanStruct = cellData;
                end

                %Assign new values into the IM2P file
                output.stimOnTimes = cleanStruct.stimONTimes;
                output.stimOnFrameNum = cleanStruct.stimOnFrameNum;
                output.stimFreqs = cleanStruct.stimFrequencies;
                output.stimDurs = cleanStruct.stimDurations;
                output.stimAmps = cleanStruct.stimAmplitudes;
                output.audioTrials = cleanStruct.audioTrials;
                output.imgDepth = cleanStruct.imagingDepth;

                %clean up & save grasping data
                numGraspON = max(size(cellData.graspTimesRight));
                numGraspOFF = max(size(cellData.graspTimesRightOff));
                graspTimes = zeros(max(numGraspON,numGraspOFF),2); %horzcat graspON & graspOFF
                graspTimes(1:numGraspON,1) = cellData.graspTimesRight;
                graspTimes(1:numGraspOFF,2) = cellData.graspTimesRightOff;
                graspDurs = zeros(min(numGraspON,numGraspOFF),1);
                graspDurs = graspTimes(1:min(numGraspON,numGraspOFF),2) - graspTimes(1:min(numGraspON,numGraspOFF),1);
                output.graspTimes = graspTimes;
                output.graspDurs = graspDurs;

                % - -------------------------------------------------------
                
                %save updated IM2P file
                cd(strcat(mouseDir,'\',dates(nd),'\',expNames(ne)));
                fprintf(strcat("Saving updated im2p file ",mouseNum(n),'-',dates(nd),'-',expNames(ne),'...'))
                tic;
                fileOutName = strcat(mouseNum(n),'-',dates(nd),'-',expNames(ne),'-im2p.mat');
                save(fileOutName,'output','-v7.3');
                clear output cellData %purge memory
                fprintf(strcat(" saved! (",num2str(toc)," sec)\n"))
                
                continue; %skip file extraction, im2p file exists
            end

            for nf = 1 : length(fileNames)
                if contains(mouseNum(n),'747')
                    % Remove blacklisted files of m747
                    if strcmp(dates(nd),'2022-08-08') && contains(fileNames(nf),'TactileStim')
                        fprintf(strcat("Loading ",mouseNum(n),'-',dates(nd),'-',fileNames(nf),'...'))
                        tic; 
                        [gfpStack, ~, scanInfo] = loadExp(fileNames(nf),0);
                        fprintf(strcat(" Took ",num2str(toc)," seconds!\n"))
                    elseif strcmp(dates(nd),'2022-08-11') && contains(fileNames(nf),'PV')
                        fprintf(strcat("Loading ",mouseNum(n),'-',dates(nd),'-',fileNames(nf),'...'))
                        tic;
                        if contains(fileNames(nf),'pre')
                            rfpImagePre = loadRfpImg(fileNames(nf),0);
                        elseif contains(fileNames(nf),'post')
                            rfpImagePost = loadRfpImg(fileNames(nf),0);
                        end
                        fprintf(strcat(" Took ",num2str(toc)," seconds!\n"))
                    elseif strcmp(dates(nd),'2022-08-11') && contains(fileNames(nf),'TactileStim')
                        fprintf(strcat("Loading ",mouseNum(n),'-',dates(nd),'-',fileNames(nf),'...'))
                        tic;
                        [gfpStack, ~, scanInfo] = loadExp(fileNames(nf),0);
                        fprintf(strcat(" Took ",num2str(toc)," seconds!\n"))
                    else %non-blacklisted experiment of m747
                        if contains(fileNames(nf),'TactileStim')
                            fprintf(strcat("Loading ",mouseNum(n),'-',dates(nd),'-',fileNames(nf),'...'))
                            tic;
                            [gfpStack, rfpTrace, scanInfo] = loadExp(fileNames(nf),1);
                            fprintf(strcat(" Took ",num2str(toc)," seconds!\n"))
                        elseif contains(fileNames(nf),'rfp-pre')
                            fprintf(strcat("Loading ",mouseNum(n),'-',dates(nd),'-',fileNames(nf),'...'))
                            tic;
                            rfpImagePre = loadRfpImg(fileNames(nf),1);
                            fprintf(strcat(" Took ",num2str(toc)," seconds!\n"))
                        elseif contains(fileNames(nf),'rfp-post')
                            fprintf(strcat("Loading ",mouseNum(n),'-',dates(nd),'-',fileNames(nf),'...'))
                            tic;
                            rfpImagePost = loadRfpImg(fileNames(nf),1);
                            fprintf(strcat(" Took ",num2str(toc)," seconds!\n"))
                        end
                    end
                else
                    if contains(fileNames(nf),'TactileStim')
                        fprintf(strcat("Loading ",mouseNum(n),'-',dates(nd),'-',fileNames(nf),'...'))
                        tic;
                        [gfpStack, rfpTrace, scanInfo] = loadExp(fileNames(nf),1);
                        fprintf(strcat(" Took ",num2str(toc)," seconds!\n"))
                    elseif contains(fileNames(nf),'rfp-pre')
                        fprintf(strcat("Loading ",mouseNum(n),'-',dates(nd),'-',fileNames(nf),'...'))
                        tic;
                        rfpImagePre = loadRfpImg(fileNames(nf),1);
                        fprintf(strcat(" Took ",num2str(toc)," seconds!\n"))
                    elseif contains(fileNames(nf),'rfp-post')
                        fprintf(strcat("Loading ",mouseNum(n),'-',dates(nd),'-',fileNames(nf),'...'))
                        tic;
                        rfpImagePost = loadRfpImg(fileNames(nf),1);
                        fprintf(strcat(" Took ",num2str(toc)," seconds!\n"))
                    end
                end
            end %end of file loop

            %find and add experiment metadata
            if specificLoadPath == 1
                %find and add experiment metadata
                cd(strcat(expMetaDir,'\',mouseNum(n),'\',dates(nd),'\'));
                [~, metaFileNames] = analyzeDir();
                cellData = load(metaFileNames(ne));
                cellData = cellData.Data;
                if iscell(cellData)
                    cleanStruct = decodeCell(cellData);
                else
                    cleanStruct = cellData;
                end

                %Assign new values into the IM2P file
                output.stimOnTimes = cleanStruct.stimONTimes;
                output.stimOnFrameNum = cleanStruct.stimOnFrameNum;
                output.stimFreqs = cleanStruct.stimFrequencies;
                output.stimDurs = cleanStruct.stimDurations;
                output.stimAmps = cleanStruct.stimAmplitudes;
                output.audioTrials = cleanStruct.audioTrials;
                output.imgDepth = cleanStruct.imagingDepth;

                %clean up & save grasping data
                numGraspON = max(size(cellData.graspTimesRight));
                numGraspOFF = max(size(cellData.graspTimesRightOff));
                graspTimes = zeros(max(numGraspON,numGraspOFF),2); %horzcat graspON & graspOFF
                graspTimes(1:numGraspON,1) = cellData.graspTimesRight;
                graspTimes(1:numGraspOFF,2) = cellData.graspTimesRightOff;
                graspDurs = zeros(min(numGraspON,numGraspOFF),1);
                graspDurs = graspTimes(1:min(numGraspON,numGraspOFF),2) - graspTimes(1:min(numGraspON,numGraspOFF),1);
                output.graspTimes = graspTimes;
                output.graspDurs = graspDurs;
            end

            %Aggregare and save data collected for this experiment
            cd(strcat(mouseDir,'\',dates(nd),'\',expNames(ne)) );
            fprintf(strcat("Saving ",mouseNum(n),'-',dates(nd),'-',expNames(ne),'...'))
            tic;
            fileOutName = strcat(mouseNum(n),'-',dates(nd),'-',expNames(ne),'-im2p.mat');
            if exist('rfpImagePre','var')
                output.rfpPVImgPre = rfpImagePre;
            end
            if exist('rfpImagePost','var')
                output.rfpPVImgPost = rfpImagePost;
            end
            if exist('rfpTrace','var')
                output.rfpExpTrace = rfpTrace;
            end
            output.gfp3FrameAvg = gfpStack;
            output.expMetaData = []; %TODO
            output.numFrames = scanInfo.num_frames;
            save(fileOutName,'output','-v7.3');
            clear output rfpImagePre rfpImagePost rfpTrace gfpStack scanInfo %purge memory
            fprintf(strcat(" saved! (",num2str(toc)," sec)\n"))
        end
    end
end

disp('Data Parser complete.');
    
    

%% visualization of gfp and rfp channels
% img1 = imread(fileNames(nf),1);
% img2 = imread(fileNames(nf),2);
% img3 = imread(fileNames(nf),3);
% img4 = imread(fileNames(nf),4);
% img5 = imread(fileNames(nf),5);
% img6 = imread(fileNames(nf),6);
% imgOdd = mean(cat(3,img1, img3, img5),3);
% imgEven = mean(cat(3,img2, img4, img6),3);
% f = figure;
% f.Position = [1056 578 850 400];
% subplot(1,2,1);
% imagesc(imgOdd);
% title(strcat(mouseNum(n)," ",dates(nd)," ",fileNames(nf)," Odd (gfp)"));
% subplot(1,2,2);
% imagesc(imgEven);
% title("Even Frames (RFP?)");
% close all;



%% Local Functions
%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
%Load in GFP images and RFP trace from Experiment 
function [gfpStack, rfpTrace, scanInfo] = loadExp(fileName,opt)
    imResized = 0;
    imInfo = imfinfo(fileName);
    scanInfo.num_frames = size(imInfo,1);
    scanInfo.Xpx = imInfo(1).Width;
    scanInfo.Ypx = imInfo(1).Height;
    if opt == 1 %skip even frames (RFP info)
        %prepopulate array to store GFP Frames
        try
            gfpStack = zeros(scanInfo.Xpx,scanInfo.Ypx,floor(scanInfo.num_frames/6));
        catch ME
            if strcat(ME.identifier,'MATLAB:array:SizeLimitExceeded')
                fprintf("Resizing...");
                imResized = 1;
                newX = floor(scanInfo.Xpx/2);
                newY = floor(scanInfo.Ypx/2);
                gfpStack = zeros(newX,newY,floor(scanInfo.num_frames/6));
            end
        end
        
        %prepopulate array to store RFP frames
        rfpTrace = zeros(floor(scanInfo.num_frames/2),1);
        
        %populate GFP Matrix
        for n = 1 : floor(scanInfo.num_frames/6)
            if imResized == 1
                img1 = imresize(double(imread( fileName,6*n-5,'Info',imInfo)),0.5);
                img2 = imresize(double(imread(fileName,6*n-3,'Info',imInfo)),0.5);
                img3 = imresize(double(imread(fileName,6*n-1,'Info',imInfo)),0.5);
                gfpStack(:,:,n) = mean(cat(3,img1,img2,img3),3);
            else
                img1 = imread(fileName,6*n-5,'Info',imInfo);
                img2 = imread(fileName,6*n-3,'Info',imInfo);
                img3 = imread(fileName,6*n-1,'Info',imInfo);
                gfpStack(:,:,n) = mean(cat(3,img1,img2,img3),3);
            end
        end
        
        for n = 1 : floor(scanInfo.num_frames/2)
           rfpTrace(n) = mean(mean(imread(fileName,2*n,'Info',imInfo))); 
        end
    elseif opt == 0
        gfpStack = zeros(scanInfo.Xpx,scanInfo.Ypx,floor(scanInfo.num_frames/3)+1);
        rfpTrace = [];
        for n = 1 : floor(scanInfo.num_frames/3)
            img1 = imread(fileName,3*n-2,'Info',imInfo);
            img2 = imread(fileName,3*n-1,'Info',imInfo);
            img3 = imread(fileName,3*n,'Info',imInfo);
            gfpStack(:,:,n) = mean(cat(3,img1,img2,img3),3);
        end
    else 
        error("var opt not recognized! must be 1 or 0");
    end
    clear img1 img2 img3 
end

%--------------------------------------------------------------------------
%Analyze a Directory for File(s) and Folder(s) names
function [folderNames, fileNames] = analyzeDir()
    myDir = dir;
    folderNames = repmat("",size(myDir,1),1);
    fileNames = repmat("",size(myDir,1),1);
    n=1;
    while n <= size(myDir,1)
        if strcmp(myDir(n).name,'.') 
            myDir(n) = []; %clear entry
            continue;
        elseif strcmp(myDir(n).name,'..')
            myDir(n) = []; %clear entry
            continue;
        elseif strcmp(myDir(n).name,'Analysis')
            myDir(n) = []; %clear entry
            continue;
        elseif myDir(n).isdir == 1
            folderNames(n) = string(myDir(n).name);
            n=n+1;
        elseif myDir(n).isdir == 0
            fileNames(n) = string(myDir(n).name);
            n=n+1;
        else
            n=n+1;
        end
    end
    fileNames(strcmp(fileNames,"")) = [];
    folderNames(strcmp(folderNames,"")) = [];
end

%--------------------------------------------------------------------------
%Decode unlabeled cell array of experiment metadata prior to 2022/08/08
function [cleanStruct] = decodeCell(cellData)
    cellLength = size(cellData,2);
    cleanStruct = struct();
    if cellLength == 8
        %fill structure
        cleanStruct.stimDurations = cellData{1,1};
        cleanStruct.stimFrequencies = cellData{1,2};
        cleanStruct.stimAmplitudes = cellData{1,3};
        cleanStruct.stimONTimes = cellData{1,4};
        cleanStruct.stimOnFrameNum = cellData{1,5};
        cleanStruct.ledTimes = cellData{1,6};
        cleanStruct.expStartTime = cellData{1,7};
        cleanStruct.expStopTime = cellData{1,8};
    end
end

%--------------------------------------------------------------------------














