% dataSurfer1.m
%
% mduhain 2023-03-25
% - Surf through all files in V:/Haptics-2P/2P-Data/ and create freq.
%   response maps from the gfpStack and Exp-Metadata
%
% mduhain 2023-04-25
% - Go through all prepared mouse data and and calculate non-ROI quadrant
% responses, and avg reesponse of all roi's within a quadrant.
%
%--------------------------------------------------------------------------

% Directory assignment
startingDir = 'C:\Users\scanimage\Desktop\prepMouseData';
%rawImgDir = 'V:\Haptics-2P\2P-Data\';
rawImgDir = 'C:\dataParser\2P-Data\';

cd(startingDir);
[mouseNum, ~] = analyzeDir();

for nm = 1 : length(mouseNum) %mouse number loop
    mouseDir = strcat(startingDir,'\',mouseNum(nm));
    cd(mouseDir);
    [dates, ~] = analyzeDir();
    for nd = 1 : length(dates) %date loop, nd = date number
        cd(strcat(mouseDir,'\',dates(nd)));
        [expNames,~] = analyzeDir();
        for ne = 1 : length(expNames) %experiment loop, ne = exp number
            cd(strcat(mouseDir,'\',dates(nd),'\',expNames(ne)));
            [depths, fileNames] = analyzeDir();
            if any(contains(fileNames,'im2p')) && any(contains(fileNames,'analyzed'))
                %Load in im2p file for Exp Metadata
                targetList = fileNames(contains(fileNames,'im2p-analyzed'));
                for nf = 1 : length(targetList)
                    load(targetList(nf)); %loads in im2p, name a0
                    a0 = a;
                    clear a;
                    
                    %create roiMask
                    for nr = 1 : a0.somaticRoiCounter
                        if nr == 1
                            roiMask = a0.somaticROIs{nr};
                        else
                            roiMask = roiMask + a0.somaticROIs{nr};
                        end
                    end
                    roiMask(roiMask == 2) = 1;
                    roiMask = (roiMask .* -1) + 1;
                    figure; imagesc(roiMask);
                    
                    %grab gfp stack
                    cd(strcat(rawImgDir,mouseNum(nm),'\',dates(nd),'\',expNames(ne)));
                    tempStr = char(strjoin(strsplit(targetList(nf),'-analyzed')));
                    tempStr(tempStr == ' ') = [];
                    load(tempStr);
                    gfpStack = a.gfp3FrameAvg;
%                     if you want to do motion correction
%                     gfpStack16 = uint16(a.gfp3FrameAvg);          % stores gfp frames from file
%                     gfpStack_mean = uint16(mean(gfpStack16,3));   % determine avg from all frames
%                     importer
%                     end
                    clear a;
                    
                    %apply mask & isolate quadrant F traces
                    roiMask(roiMask == 0) = NaN;
                    gfpStackNoROI = gfpStack .* roiMask;
                    
                    q1 = gfpStackNoROI(1:end/2,1:end/2,:);
                    q1F = squeeze(mean(q1,[1 2],'omitnan'));
                    q2 = gfpStackNoROI(end/2+1:end,1:end/2,:);
                    q2F = squeeze(mean(q2,[1 2],'omitnan'));
                    q3 = gfpStackNoROI(1:end/2,end/2+1:end,:);
                    q3F = squeeze(mean(q3,[1 2],'omitnan'));
                    q4 = gfpStackNoROI(end/2+1:end,end/2+1:end,:);
                    q4F = squeeze(mean(q4,[1 2],'omitnan'));
                    
                    clear q1 q2 q3 q4
                    
                    %plot the traces form each quadrant background signal
                    figure;
                    a0.getTraces(q1F',a0.stimOnFrameNum,a0.stimFreqs,a0.audioTrials,1);
                    figure;
                    a0.getTraces(q2F',a0.stimOnFrameNum,a0.stimFreqs,a0.audioTrials,1);
                    figure;
                    a0.getTraces(q3F',a0.stimOnFrameNum,a0.stimFreqs,a0.audioTrials,1);
                    figure;
                    a0.getTraces(q4F',a0.stimOnFrameNum,a0.stimFreqs,a0.audioTrials,1);
                  
                    %average all ROI's in quadrant
                    backMask = roiMask;
                    backMask(backMask == 1) = 0;
                    backMask(isnan(backMask)) = 1;
                    backMask(backMask == 0) = NaN;
                    gfpStackROI = gfpStack .* backMask;
                    
                    q1 = gfpStackROI(1:end/2,1:end/2,:);
                    q1roiF = squeeze(mean(q1,[1 2],'omitnan'));
                    q2 = gfpStackROI(end/2+1:end,1:end/2,:);
                    q2roiF = squeeze(mean(q2,[1 2],'omitnan'));
                    q3 = gfpStackROI(1:end/2,end/2+1:end,:);
                    q3roiF = squeeze(mean(q3,[1 2],'omitnan'));
                    q4 = gfpStackROI(end/2+1:end,end/2+1:end,:);
                    q4roiF = squeeze(mean(q4,[1 2],'omitnan'));
                    clear q1 q2 q3 q4
                    
                end    
                    
            end 
            
            
        end
        
    end
    
end


