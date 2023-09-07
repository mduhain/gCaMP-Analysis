% quickResponseMaps.m
%
% script from @mduhain, 2023-03-24

% Inputs
fSignal = double(gfpStack);
stimFrameTimes = a.stimOnFrameNum;
stimFreqs = a.stimFreqs;
audioTrials = a.audioTrials;
    
preStimT = 1; %sec before sitm on;
postStimT = 6; %sec after stim on;
frameDur = 0.1;
Hz = 10; %frameRate
allFreqs = unique(stimFreqs);
numTrials = length(stimFreqs);
keyMap = zeros(numTrials,length(allFreqs));
for n=1:length(allFreqs)
    keyMap(:,n) = stimFreqs == allFreqs(n);
end

keyMap = logical(keyMap);
xRange = [-1*preStimT:frameDur:postStimT-frameDur];
frameBin = [-1*Hz:1:(postStimT)*Hz];
frameBin(frameBin==0) = [];

%% 
tic;

figure;
hold on;

for i = 1:length(allFreqs)
    thisFreqOnTimes = floor(stimFrameTimes(keyMap(:,i))/3);
    big = NaN(size(fSignal,1),size(fSignal,2),size(frameBin,2),length(thisFreqOnTimes));
    
    for n=1:length(thisFreqOnTimes)
        if thisFreqOnTimes(n)+(postStimT*Hz) > size(fSignal,3)
            disp("-1");
            break;
        end
        if thisFreqOnTimes(n)-(preStimT*Hz) < 0
            disp("-1");
            continue;
        end
        raw = fSignal(:,:,thisFreqOnTimes(n)+frameBin);
        meanPre = mean(raw(:,:,1:preStimT*Hz),3);
        big(:,:,:,n) = raw - meanPre;
    end
    
    clean = squeeze(mean(big,4,'omitnan'));
    respMap = squeeze(mean(clean(:,:,11:30),3));
    locMax = max(max(max(respMap)));
    locMin = min(min(min(respMap)));
    
    respMapPos = respMap + (-1*locMin);
    normVal = 2/max(max(respMapPos));
    respMapNorm = respMapPos.*normVal-1;
    
    subplot(2,3,i);
    imagesc(respMapNorm);
    title(strcat("Norm. Resp. @ ",num2str(allFreqs(i))," Hz."));
    colorbar;
    
    clear raw clean
end 

toc;
