
% Meta analysis of correlations between cell pairs

% 1: Tuning curve best fit (Gaussian, Polynomial, Exponential)
% 2: Gaussian Tuning Width
% 3: R^2 value for goodness of fit
% 4: Best Frequency for neuron
% 5: Cell ID (GFP-Only vs. Overlap)
% 6: cfit struct for Col. 1-3
% 7: Struct for fitting parameters
% 8: Response struct contaiting raw data, trial traces, etc 
% 9: Mouse number
% 10: Experiment Data
% 11: Experiment Depth
% 12: Number of freqs with significant response
% 13: [X,Y] location of somal center
% 14: Mouse Exp ID
% 15: Neuron responsive yes/no
% 16: Neuron selective yes/no
% 17: Response type (Pos,Neg,Mix)

% hardcode values
freqList = [100;300;500;700;900;1100];
load('PreTrainingData.mat');
sessNames = unique([bigOut{:,14}]);
corrMatrix = cell(length(sessNames),2);

for ns = 1 : length(sessNames) %number of sessions
    corrMatrix{ns,1} = sessNames(ns);
    BOid = find([bigOut{:,14}] == sessNames(ns));
    corrMat = zeros(length(BOid),length(BOid),2); %cell,cell,[R P]
    for nn = 1 : length(BOid) %number of first dim (corr matrix)
        for nj = nn+1 : length(BOid) % number of second dimension (corr matrix)
            i1 = BOid(nn); %index number
            i2 = BOid(nj); %index number
            %reconstruct complete df/f trace for all frequencies & reps
            rawTrace1 = [];
            rawTrace2 = [];
            for nf = 1 : length(freqList) %number frequencies
                fn1 = strcat('f',num2str(freqList(nf))); %field name here 'f000'
                rawTrace1 = [rawTrace1; bigOut{i1,8}.(fn1).allTraces(:)];
                rawTrace2 = [rawTrace2; bigOut{i2,8}.(fn1).allTraces(:)];
            end
            [RHO,PVAL] = corr([rawTrace1 rawTrace2]);
            corrMat(nj,nn,1) = RHO(1,2);
            corrMat(nj,nn,2) = PVAL(1,2);
        end
    end

    corrMatrix{ns,2} = corrMat;

    figure; hold on;
    imagesc(corrMat(:,:,1));
    xlabel("ROI N");
    ylabel("ROI M");
    colorbar;
    sc = strsplit(sessNames(ns),'_');
    title(strcat("Pearson correlation matrix for neuron pairs from ",sc(1),"-",sc(2)));
    hold off;

end


%% Compare all correlations of PV:PV, Non-PV:Non-PV, and PV:Non-PV

for ns = 1 : length(sessNames) %number of sessions
    BOid = find([bigOut{:,14}] == sessNames(ns));
    isPV = false(length(BOid),1);
    for nn = 1 : length(BOid) %number of first dim (corr matrix)
        if strcmp(bigOut{nn,5},'Overlap')
            isPV(nn) = true;
        end
    end



end






























