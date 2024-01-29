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
% mduhain 2023-08-08
% - Reformatting this code to act as a general converter from a Directory
% of im2p files into a large metaDataStruct:
%
% 1. CELLID
%
% 2. ROI X Y
%
% 3. ROI MASK
%
% 4. RESPONSE STRUCT (from getTraces)
%
% 5. RAW TRACE
%
% 6. EXP DEPTH
%
% 7. MOUSE EXP ID
%
% 8. IS_RESPONSIVE
%
% 9. IS_SELECTIVE
%
% 10. PREFERED FREQUENCY

%% ------------------------------------------------------------------------

% Directory assignment
% startingDir = 'C:\Users\ramirezlab\Desktop\prepSomData';
% 
clear all; close all; clc
startingDir = pwd;
cd(startingDir);
[mouseNum, reportNames] = analyzeDir();

% Find and load most recent report
appendingReportFlag = 0;
if isempty(reportNames)
    disp("No previous reports found: Starting Analysis...");
    allFileNames = ["blank","test1"];
else
    disp("Which report should be loaded...?");
    [selectedReport,~] = uigetfile();
    disp(strcat("Loading: ",selectedReport," ..."));
    R = load(selectedReport);
    appendingReportFlag = 1;
    % loop through to identify dates in report
    allFileNames = repmat("",length(R.bigOut));
    for n = 1 : length(R.bigOut)
        allFileNames(n) = R.bigOut{n,4}.fileSource;
    end
    allFileNames = unique(allFileNames);
end

% Output Structure Creation
bigOut = cell(10000,10);
boc = 1; %bigOut Counter;

% Hardcoded Variables
freqList = [100;300;500;700;900;1100];


%% MAIN LOOP 

for nm = 1 : length(mouseNum) %mouse number loop
    mouseDir = strcat(startingDir,'\',mouseNum(nm));
    cd(mouseDir);
    [dates, ~] = analyzeDir();
    mouseUniqueExps = 1;
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
                    %check if experiment is in the report
                    if any(contains(allFileNames,targetList(nf)))
                        disp(strcat("Skipping: ",targetList(nf)));
                        continue;
                    end
                    %Experient is not in the report, proceed with analysis
                    load(targetList(nf)); %loads in im2p, name a0
                    a0 = a;
                    clear a;
                    imgSizeFlag = 0;
                    for ng = 1 : a0.somaticRoiCounter % green ROI loop, ng = green ROI number 
                        %% CELL ID / ROI X Y / ROI MASK
                        thisG = a0.somaticROIs{1,ng};
                        for nr = 1 : a0.redSomaticRoiCounter %red ROI loop, nr = red ROI number
                            if ~iscell(a0.redSomaticROIs) %no red ROIS for this Exp
                                bigOut{boc,1} = 0; %GFP only
                                break
                            end
                            thisR = a0.redSomaticROIs{1,nr};
                            %check for size discrepancies between images
                            if size(thisG) > size(thisR)
                                thisG = imresize(thisG,0.5);
                                bigOut{boc,2} = a0.somaticROICenters{ng}.Centroid./2;
                                bigOut{boc,3} = thisG;
                                imgSizeFlag = 1;
                            end
                            %sometimes RFP image hasnt yet been resized 
                            if size(thisR) > size(thisG)
                                thisR = imresize(thisR,0.5);
                            end
                            %if overlap, calculate % overlap
                            if any((thisR + thisG) > 1,'all')  %ROI overlap found
                                bigOut{boc,1} = mean([nnz(thisR&thisG)/nnz(thisG),nnz(thisR&thisG)/nnz(thisR)]); %percent overlap
                                break
                            end
                            if nr == a0.redSomaticRoiCounter %final loop, no overlap
                                bigOut{boc,1} = 0;
                            end
                        end
                        
                        if all(isnan(a0.somaticFDT(ng,:)))
                            disp(strcat(" Skipping ROI ",num2str(ng),", it's empty..."));
                            continue
                        end
    
                        if imgSizeFlag == 0
                            %not assigned earlier, assign now
                            bigOut{boc,2} = a0.somaticROICenters{ng}.Centroid;
                            bigOut{boc,3} = a0.somaticROIs{ng};
                        end

                        %% RESP STRUCT & RAW TRACE
                        if ~isnan(a0.somaticFDT) %FDT present
                            if ~all(isnan(a0.somaticFDT(ng,:)))
                                out = a0.getTraces(a0.somaticFDT(ng,:),a0.stimOnFrameNum,a0.stimFreqs,a0.audioTrials,1);
                                bigOut{boc,4} = out;
                                pause(1);
                                close all;
                                bigOut{boc,5} = a0.somaticFDT(ng,:);
                            end
                        else
                            %whole field missing, use regular flourescence trace non-de-trended
                            if ~all(isnan(a0.somaticF(ng,:)))
                                bigOut{boc,4} = a0.getTraces(a0.somaticF(ng,:),a0.stimOnFrameNum,a0.stimFreqs,a0.audioTrials,1);
                                pause(1);
                                close all;
                                bigOut{boc,5} = a0.somaticF(ng,:);
                            end
                        end
    
                        %% EXP DEPTH
                        bigOut{boc,6} = a0.imgDepth;

                        %% MOUSE EXP ID
                        bigOut{boc,7} = strcat(mouseNum(nm),'_',num2str(mouseUniqueExps));
                        bigOut{boc,4}.fileSource = targetList(nf);   
                        bigOut{boc,4}.stimOnFrameNum = a0.stimOnFrameNum;
                        bigOut{boc,4}.stimFreqs = a0.stimFreqs;

                        %% RESPONSIVE / SELECTIVE / PREF FREQ
                        isRespF = out.respFreqZP(:,2)<0.05;
                        if any(isRespF)
                            bigOut{boc,8} = 1; %neuron is responsive
                            %now calculate if any 2 responsive frequencies
%                             if sum(isRespF) > 1 % if 2 or more are responsive
%                             respFreqs = find(isRespF);
%                             selFlag = 0; %selectivity found flag

%                             for n1 = 1 : length(respFreqs)
%                                 firstF = strcat('f',num2str(freqList(respFreqs(n1))));
%                                 resps1 = out.(firstF).avgRespPostStim(:);
%                                 for n2 = n1+1 : length(respFreqs)
%                                     secondF = strcat('f',num2str(freqList(respFreqs(n2))));
%                                     resps2 = out.(secondF).avgRespPostStim(:);
%                                     pVal = ranksum(resps1,resps2)
% %                                     figure; plot([1 2], [resps2 resps1],'k.'); hold on;
% %                                     ca = gca; ca.XLim = [0 3]; hold off;
%                                     if pVal < 0.05
%                                         bigOut{boc,9} = 1; %neuron is selective
%                                         selFlag = 1;
%                                         break
%                                     end
%                                 end
%                                 if selFlag == 1
%                                     break
%                                 end
%                             end
                            if out.isSelectiveP < 0.05
                                bigOut{boc,9} = 1; %neuron is selective
                            else
                                bigOut{boc,9} = 0; %neuron is non-selective
                            end
%                             end
                        else
                             bigOut{boc,8} = 0; %neuron is non-responsive
                        end
                        

                        %calculate frequency with best average response
                        avgResps = bigOut{boc,4}.avgRespPost;
                        bigOut{boc,10} = freqList(avgResps == max(avgResps));

                        %% LOOP CLEANUP
                        %increment up bigOut counter
                        boc = boc + 1;

                        %command window counter
                        if boc == 1000
                            clc;
                        else
                            disp(strcat("roi",num2str(boc)));
                        end

                    end % ng
                    mouseUniqueExps = mouseUniqueExps + 1;
                end %nf  
                    
            end %if contains im2p-analyzed
            
        end %ne
        
    end %nd
    
end %nm

% append bigOut onto R.bigOut
if appendingReportFlag == 1
    boc = size(R.bigOut,1)+1;
    for n = 1 : size(bigOut,1)
        for m = 1 : size(bigOut,2)
            R.bigOut{boc,m} = bigOut{n,m};
        end
        boc = boc + 1;
    end
end






%% SAVE COPY OF bigOut

disp("Exporting large cell array to starting Dir...");
tic;
cd(startingDir);
dt = datetime('now');
fileOutName = strcat('surfReport_',num2str(dt.Year),'-',num2str(dt.Month),...
    '-',num2str(dt.Day),'_',num2str(dt.Hour),'-',num2str(dt.Minute),'-',...
    num2str(round(dt.Second)),'.mat');
save(fileOutName,'bigOut');
toc;


%% for temp conversion

% freqNames = ["f100","f300","f500","f700","f900","f1100"];
% load('dataSurfed2023-08-09.mat');
% for n = 1 : size(bigOut,1)
%     if ~isempty(bigOut{n,1}) %populated trial
%         isRespF = bigOut{n,4}.respFreqZP(:,2)<0.05;
%         if any(isRespF)
%             bigOut{boc,8} = 1; %neuron is responsive
%             %now calculate if any 2 responsive frequencies
%             indxF = find(isRespF);
%             for x1 = 1 : sum(isRespF)
%                 for x2 = x1 + 1 : sum(isRespF)
%                     [P,H] = ranksum(bigOut{n,4}.(freqNames(x1)).avgRespPostStim,...
%                         bigOut{n,4}.(freqNames(x2)).avgRespPostStim); 
%                 end
%             end
% 
%             %have a significant difference
%         end
%     end
% end
% 
% % temp plotting section
% figure; hold on
% plot(1,bigOut{1,4}.f100.avgRespPostStim,'k.');
% plot(2,bigOut{1,4}.f300.avgRespPostStim,'k.');
% plot(3,bigOut{1,4}.f500.avgRespPostStim,'k.');
% plot(4,bigOut{1,4}.f700.avgRespPostStim,'k.');
% plot(5,bigOut{1,4}.f900.avgRespPostStim,'k.');
% plot(6,bigOut{1,4}.f1100.avgRespPostStim,'k.');
% plot(1,mean(bigOut{1,4}.f100.avgRespPostStim),'ro');
% plot(2,mean(bigOut{1,4}.f300.avgRespPostStim),'ro');
% plot(3,mean(bigOut{1,4}.f500.avgRespPostStim),'ro');
% plot(4,mean(bigOut{1,4}.f700.avgRespPostStim),'ro');
% plot(5,mean(bigOut{1,4}.f900.avgRespPostStim),'ro');
% plot(6,mean(bigOut{1,4}.f1100.avgRespPostStim),'ro');
% %plot([1:6],bigOut{1,4}. 
% ax=gca; ax.XLim=[0,7]; hold off;



%% Quick Analysis 2023-09-11

isSom = zeros(length(bigOut),1);
for n = 1 : length(bigOut)
    if ~isempty(bigOut{n,1})
        if string(bigOut{n,1}) == "Overlap"
            isSom(n) = 1;
        end
    end
    if isempty(bigOut{n,9})
        bigOut{n,9} = 0;
    end
end

isSom = logical(isSom);
numSom = sum(isSom);
numRespSom = sum([bigOut{isSom,8}]);
numNonRespSom = numSom - numRespSom;
numSelSom = sum([bigOut{isSom,9}]);
numNonSelSom = numSom - numSelSom;
numNonSom = sum(isSom == 0);
numRespNonSom = sum([bigOut{~isSom,8}]);
numNonRespNonSom = numNonSom - numRespNonSom;
numSelNonSom = sum([bigOut{~isSom,9}]);
numNonSelNonSom = numNonSom - numSelNonSom;

f = figure; hold on;
b0 = bar(1,[numRespSom/numSom numNonRespSom/numSom],'grouped');
b0(1,1).FaceColor = [0 0.4470 0.7410];
b0(1,2).FaceColor = [0 0.4470 0.7410] .* 0.5;
stdVal0 = mannyPropSTD(numRespSom,numSom);
text(b0(1,1).XEndPoints,b0(1,1).YEndPoints+stdVal0,string(strcat("n=",num2str(numRespSom))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b0(1,2).XEndPoints,b0(1,2).YEndPoints+stdVal0,string(strcat("n=",num2str(numNonRespSom))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

b1 = bar(2,[numSelSom/numSom numNonSelSom/numSom],'grouped');
b1(1,1).FaceColor = [0.4660 0.6740 0.1880];
b1(1,2).FaceColor = [0.4660 0.6740 0.1880] .* 0.5;
stdVal1 = mannyPropSTD(numSelSom,numSom);
text(b1(1,1).XEndPoints,b1(1,1).YEndPoints+stdVal1,string(strcat("n=",num2str(numSelSom))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b1(1,2).XEndPoints,b1(1,2).YEndPoints+stdVal1,string(strcat("n=",num2str(numNonSelSom))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

b2 = bar(3,[numRespNonSom/numNonSom numNonRespNonSom/numNonSom],'grouped');
b2(1,1).FaceColor = [0 0.4470 0.7410];
b2(1,2).FaceColor = [0 0.4470 0.7410] .* 0.5;
stdVal2 = mannyPropSTD(numRespNonSom,numNonSom);
text(b2(1,1).XEndPoints,b2(1,1).YEndPoints+stdVal2,string(strcat("n=",num2str(numRespNonSom))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b2(1,2).XEndPoints,b2(1,2).YEndPoints+stdVal2,string(strcat("n=",num2str(numNonRespNonSom))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

b3 = bar(4,[numSelNonSom/numNonSom numNonSelNonSom/numNonSom],'grouped');
b3(1,1).FaceColor = [0.4660 0.6740 0.1880];
b3(1,2).FaceColor = [0.4660 0.6740 0.1880] .* 0.5;
stdVal3 = mannyPropSTD(numSelNonSom,numNonSom);
text(b3(1,1).XEndPoints,b3(1,1).YEndPoints+stdVal3,string(strcat("n=",num2str(numSelNonSom))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b3(1,2).XEndPoints,b3(1,2).YEndPoints+stdVal3,string(strcat("n=",num2str(numNonSelNonSom))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

errorbar(b0(1,1).XEndPoints,b0(1,1).YEndPoints,mannyPropSTD(numRespSom,numSom),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b0(1,2).XEndPoints,b0(1,2).YEndPoints,mannyPropSTD(numNonRespSom,numSom),...
    Color=[0,0,0],LineWidth=1.4);

errorbar(b1(1,1).XEndPoints,b1(1,1).YEndPoints,mannyPropSTD(numSelSom,numSom),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b1(1,2).XEndPoints,b1(1,2).YEndPoints,mannyPropSTD(numNonSelSom,numSom),...
    Color=[0,0,0],LineWidth=1.4);

errorbar(b2(1,1).XEndPoints,b2(1,1).YEndPoints,mannyPropSTD(numRespNonSom,numNonSom),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b2(1,2).XEndPoints,b2(1,2).YEndPoints,mannyPropSTD(numNonRespNonSom,numNonSom),...
    Color=[0,0,0],LineWidth=1.4);

errorbar(b3(1,1).XEndPoints,b3(1,1).YEndPoints,mannyPropSTD(numSelNonSom,numNonSom),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b3(1,2).XEndPoints,b3(1,2).YEndPoints,mannyPropSTD(numNonSelNonSom,numNonSom),...
    Color=[0,0,0],LineWidth=1.4);

ax0 = gca;
ax0.XTick = [1.5 3.4];
ax0.XTickLabel = {'SOM';'Non-SOM'};
ax0.YLim = [0 1];
ylabel('Percentage of Neurons')
title('Catagorization of Neurons ');
legend({'Responsive','Non-Responsive','Selective','Non-Selective'});


%% Frequency Tuning Distribution
f2 = figure;
somPrefFreq = [bigOut{isSom,10}];
somPrefFreq(somPrefFreq == 0) = [];
subplot(1,2,1);
bar(100:200:1100,histcounts(somPrefFreq,6)./numRespSom);
title(strcat("Pref. Freq. for ",num2str(numRespSom)," Som Neurons"));
xlabel('Frequency (Hz)');
ylabel('Percentage of Neurons');

nonsomPrefFreq = [bigOut{~isSom,10}];
nonsomPrefFreq(nonsomPrefFreq == 0) = [];
subplot(1,2,2);
bar(100:200:1100,histcounts(nonsomPrefFreq,6)./numRespNonSom);
title(strcat("Pref. Freq. for ",num2str(numRespNonSom)," Non-Som Neurons"));
xlabel('Frequency (Hz)');
ylabel('Number of Neurons');














