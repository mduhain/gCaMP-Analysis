% dataSurfer2.m
%
% mduhain 2023-04-26
% - Are PV neurons responsive and tuned to vibrotactile frequency in a
%   mouse with cre-dependent GCaMP7c?
%
%--------------------------------------------------------------------------

% Directory assignment
startingDir = 'V:\Haptics-2P\2P-Data\';
%rawImgDir = 'V:\Haptics-2P\2P-Data\';
rawImgDir = 'C:\dataParser\2P-Data\';

%hardcode values
freqList = [100;300;500;700;900;1100];

%BIG OUT
bigOut = cell(100,5);
boc = 1; %bigOut Counter;

cd(startingDir);
[mouseNum, ~] = analyzeDir();

for nm = 1 : length(mouseNum) %mouse number loop
    if ~strcmp(mouseNum(nm),'m792')
        disp(strcat("Skipping ",mouseNum(nm),"..."));
        continue
    end
    mouseDir = strcat(startingDir,'\',mouseNum(nm));
    cd(mouseDir);
    [dates, ~] = analyzeDir();
    for nd = 1 : length(dates) %date loop, nd = date number
        if strcmp(dates(nd),"ForepawMap")
            disp("Skipping ForepawMap...");
            continue;
        end 
        cd(strcat(mouseDir,'\',dates(nd)));
        [expNames,~] = analyzeDir();
        for ne = 1 : length(expNames) %experiment loop, ne = exp number
            cd(strcat(mouseDir,'\',dates(nd),'\',expNames(ne)));
            [depths, fileNames] = analyzeDir();
            if any(contains(fileNames,'im2p')) && any(contains(fileNames,'analyzed'))
                %Load in im2p file for Exp Metadata
                targetList = fileNames(contains(fileNames,'im2p-analyzed'));
                for nf = 1 : length(targetList) %im2p file loop, nf = file number
                    disp(strcat("loading ",targetList(nf)));
                    load(targetList(nf)); %loads in im2p
                    
                    % For each red ROI's, find overlap ROI in GREEN
                    for nr = 1 : a.redSomaticRoiCounter % red ROI loop, nr = red ROI number
                        thisR = a.redSomaticROIs{1,nr};
                        for ng = 1 : a.somaticRoiCounter %green ROI loop, ng = green ROI number
                            thisG = a.somaticROIs{1,ng};
                            if nnz(thisR & thisG) > 0 %ROI overlap
                                out = a.getTraces(a.somaticFDT(ng,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,1);
                                disp(strcat("PV Selectivity = ",num2str(out.selectivity_bestMinWorst)));
                                avgResps = NaN(floor(size(out.f100.avgRespPostStim,1)*1.1),size(out.respFreqZP,1)); %pad +10%
                                %find frequency names in fieldnames of struct 'out'
                                expKey = regexp(fieldnames(out),'f+\d*'); %regex for letter 'f' followed by any number
                                expKey =  ~cellfun(@isempty,expKey); %logical arrray for matches
                                freqs = fieldnames(out);
                                freqs = freqs(expKey);
                                for na = 1 : size(out.respFreqZP,1)
                                    avgResps(1:size(out.(freqs{na}).avgRespPostStim,1),na) = out.(freqs{na}).avgRespPostStim;
                                end
                                
                                %Gaussian Fit
                                [Fg,GOFg] = fit(freqList,out.avgRespPost,'gauss1');
%                                 figure; hold on
%                                 for np = 1 : 6
%                                 plot(freqList(np),avgResps(:,np),'k.');
%                                 plot(freqList(np),out.avgRespPost(np),'bo');
%                                 end
%                                 plot(Fg)
%                                 cf = gca;
%                                 cf.XLim = [0,1200];
%                                 xlabel('Frequency (Hz)');
%                                 ylabel('dF/F');
%                                 title(strcat("Gaussian fit to mean, R^2= ",num2str(GOFg.rsquare)));
                                sigma = 1/(Fg.a1*sqrt(pi*2));

                                %Exponential Fit
                                [Fe,GOFe] = fit(freqList,out.avgRespPost,'exp1');
%                                 figure; hold on
%                                 for np = 1 : 6
%                                 plot(freqList(np),avgResps(:,np),'k.');
%                                 plot(freqList(np),out.avgRespPost(np),'bo');
%                                 end
%                                 plot(Fe)
%                                 cf = gca;
%                                 cf.XLim = [0,1200];
%                                 xlabel('Frequency (Hz)');
%                                 ylabel('dF/F');
%                                 title(strcat("Exponential fit to mean, R^2= ",num2str(GOFe.rsquare)));
                                %sigma = 1/(Fe.a1*sqrt(pi*2));

                                %Exponential Fit
                                [Fp,GOFp] = fit(freqList,out.avgRespPost,'poly2');
%                                 figure; hold on
%                                 for np = 1 : 6
%                                 plot(freqList(np),avgResps(:,np),'k.');
%                                 plot(freqList(np),out.avgRespPost(np),'bo');
%                                 end
%                                 plot(Fp)
%                                 cf = gca;
%                                 cf.XLim = [0,1200];
%                                 xlabel('Frequency (Hz)');
%                                 ylabel('dF/F');
%                                 title(strcat("Exponential fit to mean, R^2= ",num2str(GOFp.rsquare)));
                                %sigma = 1/(Fp.a1*sqrt(pi*2));
                                

                                if GOFg.rsquare > GOFe.rsquare && GOFg.rsquare > GOFp.rsquare
                                    %Gaussian best fit
                                    bigOut{boc,1} = 'G';
                                    bigOut{boc,2} = sigma;
                                    bigOut{boc,3} = GOFg.rsquare;
                                    bigOut{boc,5} = Fg;
                                    bigOut{boc,6} = GOFg;
                                end

                                if GOFe.rsquare > GOFg.rsquare && GOFe.rsquare > GOFp.rsquare
                                    %Exponential best fit
                                    bigOut{boc,1} = 'E';
                                    bigOut{boc,2} = NaN;
                                    bigOut{boc,3} = GOFe.rsquare;
                                    bigOut{boc,5} = Fe;
                                    bigOut{boc,6} = GOFe;
                                end

                                if GOFp.rsquare > GOFg.rsquare && GOFp.rsquare > GOFe.rsquare
                                    %Polynomail best fit
                                    bigOut{boc,1} = 'P';
                                    bigOut{boc,2} = NaN;
                                    bigOut{boc,3} = GOFp.rsquare;
                                    bigOut{boc,5} = Fp;
                                    bigOut{boc,6} = GOFp;
                                end

                                [m,indx] = max(out.avgRespPost);
                                bigOut{boc,4} = freqList(indx);
                                bigOut{boc,7} = out;
                                boc = boc + 1;
                                

                                



                            end
                            
                        end
                        
                    end
 
                end    
                    
            end 
            
        end
        
    end
    
end


