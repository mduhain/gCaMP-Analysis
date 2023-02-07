%meta analysis of all data

%% First Stage (individual im2p files)
meta810_11_02_exp1 = struct();
for n=1:a.somaticRoiCounter
    name = strcat("ROI",num2str(n));
    meta810_11_02_exp1.(name) = a.getTraces(a.somaticFDT(n,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,1);
    %close all
end

for n=1:a.redSomaticRoiCounter
    name = strcat("redROI",num2str(n));
    meta810_11_02_exp1.(name) = a.getTraces(a.redSomaticFDT(n,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,1);
    %close all
end

load("Pass3.mat");
varNames = fieldnames(meta);
numFreqs = 6;
freqs = [100,300,500,700,900,1100];
numExps = size(varNames,1);

%Dimensions: EXP# / ROI# / Freq.
allZTuning = nan(numExps,100,numFreqs);
allMeanTuning = nan(numExps,100,numFreqs);
allMedTuning = nan(numExps,100,numFreqs);
allResponsivity = nan(numExps,100,1);
allSelectivity = nan(numExps,100,1);
totalROIs = 0;

%% SECOND STAGE
for i = 1 : size(varNames,1)
    %Vars
    expName = varNames{i};
    numROIs = size(fieldnames(meta.(expName)),1);
    flag = 0;
    flagVal = 0;
    %Extraction loop
    for n=1:numROIs
        logIndx = contains(fieldnames(meta.(expName)),'red');
        if logIndx(n) == 1
            if flag == 0
                flagVal = n-1;
                flag = 1;
            end
            fName = strcat("redROI",num2str(n-flagVal));
        else
            %fName = strcat("ROI",num2str(n));
            continue;
        end
        allZTuning(i,n,:) = meta.(expName).(fName).respFreqZP(:,1)';
        try
            allResponsivity(i,n) = meta.(expName).(fName).responsivityZ;
        catch
            %not important anyway
        end
        allSelectivity(i,n) = meta.(expName).(fName).selectivity_bestMinWorst;
        allMeanTuning(i,n,:) = meta.(expName).(fName).avgRespPost;
        allMedTuning(i,n,:) = meta.(expName).(fName).medRespPost;
    end
    totalROIs = totalROIs + numROIs;
end

%% PLOTTING
%all tuning curves
totTune = allZTuning(~isnan(allZTuning));
totTune = reshape(totTune,[size(totTune,1)/numFreqs,numFreqs]);
[maxZ, maxZindx] = max(totTune,[],2);
[minZ, minZindx] = min(totTune,[],2);
numPosZUnits = size(maxZ(maxZ>=1.96),1);
numNegZUnits = size(minZ(minZ<=-1.96),1);
otherZUnits = size(totTune,1) - (numNegZUnits + numPosZUnits);

%Pie chart of significantly modulated units
figure;
pie([numPosZUnits numNegZUnits otherZUnits])
legend({'Excited','Inhibited','Unresponsive'});
title(strcat("Neurons (n=",num2str(size(totTune,1)),") with significant response to at least one frequency"));

sigMax = maxZ >= 1.96;
sigMin = minZ <= -1.96;
[~, maxZindxs] = max(totTune(sigMax,:),[],2);
[~, minZindxs] = min(totTune(sigMin,:),[],2);

%histogram of optimal frequency positive mod
figure;
C = categorical(maxZindxs,[1 2 3 4 5 6],{'100','300','500','700','900','1100'})
histogram(C);
title("Optimal frequency for all neurons");
subtitle("From neurons with significant positive modulation")
xlabel("Frequency (Hz)");
ylabel("Number of ROIs");

%histogram of optimal frequency negative mod
figure;
C = categorical(minZindxs,[1 2 3 4 5 6],{'100','300','500','700','900','1100'})
histogram(C);
title("Optimal frequency for all neurons");
subtitle("From neurons with significant negative modulation")
xlabel("Frequency (Hz)");
ylabel("Number of ROIs");


disp(" ");
disp(strcat("Total ROIs analyzed: ",num2str(numROIs)));
numSelectiveElements = allSelectivity;
disp(strcat("Of ",num2str(numROIs)," ROIs, ",num2str(numPosZUnits), ""));


%SOME PLOTS
figure;
bar(sort(allSelectivity)); hold on;
title("Selectivity of all ROIs");
xlabel("ROI number");
ylabel("Selectivity index");


figure;
[M,I] = max(ZTuning');
C = categorical(I,[1 2 3 4 5 6],{'100','300','500','700','900','1100'})
histogram(C);
title("Optimal Freq. for all ROIs");
xlabel("Frequency (Hz)");
ylabel("Number of ROIs");

figure;
bar(sort(Responsivity)); hold on;
plot([0 numROIs],[1.96 1.96],'g--');
plot([0 numROIs],[-1.96 -1.96],'g--');
title("Responsivity of all ROIs");
xlabel("Number of ROIs");
ylabel("Z-Stat");







%VISUALIZATIN ONLY
for n=1:a.somaticRoiCounter
    a.getTraces(a.somaticF(n,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,1);
end




%FINDING
for n=1:length(fieldnames(meta810_11_02_exp1))
    name = strcat("ROI",num2str(n));
    disp(strcat("N:",num2str(n)," sel: ",num2str(meta810_11_02_exp1.(name).selectivity)));
end

figure;
shadedErrorBar(meta810_11_02_exp1.ROI10.xRange,meta810_11_02_exp1.ROI10.f900.avgTrace,...
    meta810_11_02_exp1.ROI10.f900.semTraces); hold on;
plot([0 0],[meta810_11_02_exp1.ROI10.yAxMinMax(1) meta810_11_02_exp1.ROI10.yAxMinMax(2)],'k--');
plot([1 1],[meta810_11_02_exp1.ROI10.yAxMinMax(1) meta810_11_02_exp1.ROI10.yAxMinMax(2)],'k--');

figure; plot(meta810_11_02_exp1.ROI10.f900.allTraces');

figure;
plot(0,meta810_11_02_exp1.ROI10.f900.preStim,'k.'); hold on;
plot(0,mean(meta810_11_02_exp1.ROI10.f900.preStim),'ro');
plot(1,meta810_11_02_exp1.ROI10.f900.avgRespPostStim,'k.');
plot(1,mean(meta810_11_02_exp1.ROI10.f900.avgRespPostStim),'ro');
ax = gca;
ax.XLim = [-1 2];


figure;
boxplot([meta810_11_02_exp1.ROI10.f900.preStim meta810_11_02_exp1.ROI10.f900.avgRespPostStim],'Notch','on', ...
        'Labels',{'pre-stim','post-stim'});




