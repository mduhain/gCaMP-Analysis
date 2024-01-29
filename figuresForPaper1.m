% figuresForPaper1.m
% 
% Script @mduhain, July 17th, 2023
% Showing the identity %'s of each neuron type (PV vs. non-PV)
% also showing % responsive (PV/non-PV) both before and after training.
% also showing % selective (PV/non-PV) both before and after training.
%
%
% BigOut Column Titles (1 row = 1 neuron)
% -------------------------------------------------------------------------
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
% -------------------------------------------------------------------------

%Load in the curated data
cd('C:\Users\ramirezlab\Desktop');
preD = load('PreTrainingData.mat');
posD = load('PostTrainingData.mat');

%frequency list
freqList = [100,300,500,700,900,1100];

% PRE TRAINING - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
isPVpre = zeros(1,length(preD.bigOut));
is810pre = zeros(1,length(preD.bigOut));
is806pre = zeros(1,length(preD.bigOut));
is792pre = zeros(1,length(preD.bigOut));
for n = 1 : length(preD.bigOut)
    preD.bigOut{n,15} = 0;
    preD.bigOut{n,16} = 0;
    if any(preD.bigOut{n,8}.respFreqZP(:,2)<0.05)
        % Neuron responsive
        preD.bigOut{n,15} = 1;
        sigResps = preD.bigOut{n,8}.respFreqZP(:,2)<0.05;
        zVals = preD.bigOut{n,8}.respFreqZP(sigResps,1);
        if any(zVals < 0) && any(zVals > 0) %Mixed responsive neuron, skip
            preD.bigOut{n,17} = 'mix';
            [~,IndxMax] = max(preD.bigOut{n,8}.respFreqZP(:,1)); %max resp
            [~,IndxMin] = min(preD.bigOut{n,8}.respFreqZP(:,1)); %min resp
            maxResps = preD.bigOut{n,8}.(strcat('f',num2str(freqList(IndxMax)))).avgRespPostStim;
            minResps = preD.bigOut{n,8}.(strcat('f',num2str(freqList(IndxMin)))).avgRespPostStim;
            [P,H] = ranksum(maxResps,minResps);
            if P < 0.05
                preD.bigOut{n,16} = 1;
            else 
                preD.bigOut{n,16} = -1;
            end
        elseif all(zVals > 0) %Positive responding neuron
            preD.bigOut{n,17} = 'pos';
            [~,IndxMax] = max(preD.bigOut{n,8}.respFreqZP(:,1)); %max resp
            [~,IndxMin] = min(preD.bigOut{n,8}.respFreqZP(:,1)); %min resp
            maxResps = preD.bigOut{n,8}.(strcat('f',num2str(freqList(IndxMax)))).avgRespPostStim;
            minResps = preD.bigOut{n,8}.(strcat('f',num2str(freqList(IndxMin)))).avgRespPostStim;
            [P,H] = ranksum(maxResps,minResps);
            if P < 0.05
                preD.bigOut{n,16} = 1;
            else 
                preD.bigOut{n,16} = -1;
            end
        elseif all(zVals < 0) %negatively responding neuron, flip vals
            preD.bigOut{n,17} = 'neg';
            [~,IndxMax] = max(abs(preD.bigOut{n,8}.respFreqZP(:,1))); %max resp with abs() because it's neg
            [~,IndxMin] = min(abs(preD.bigOut{n,8}.respFreqZP(:,1))); %min resp with abs() because it's neg
            maxResps = preD.bigOut{n,8}.(strcat('f',num2str(freqList(IndxMax)))).avgRespPostStim;
            minResps = preD.bigOut{n,8}.(strcat('f',num2str(freqList(IndxMin)))).avgRespPostStim;
            [P,H] = ranksum(maxResps,minResps);
            if P < 0.05
                preD.bigOut{n,16} = 1;
            else 
                preD.bigOut{n,16} = -1;
            end
        end
    else
        %neuron non-responsive
        preD.bigOut{n,15} = 0;
    end
    %check ROI Overlap (PV+ / SOM+)
    if strcmp(preD.bigOut{n,5},'Overlap')
        isPVpre(n) = 1;
    end
    %assign mouse numbers
    if strcmp(preD.bigOut{n,9},'m810')
        is810pre(n) = 1;
    elseif strcmp(preD.bigOut{n,9},'m806')
        is806pre(n) = 1;
    elseif strcmp(preD.bigOut{n,9},'m792')
        is792pre(n) = 1;
    else
        error(strcat("Mouse number, on row ",num2str(n)," not recognized!"));
    end
end
isResponsivePreTrain = logical([preD.bigOut{:,15}]);
isSelectivePreTrain = logical([preD.bigOut{:,16}]);
is810pre = logical(is810pre);
is806pre = logical(is806pre);
is792pre = logical(is792pre);
isPVpre = logical(isPVpre);
isPVpre = logical(isPVpre);

% POST TRAINING - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
isPVpost = zeros(1,length(posD.bigOut));
is810post = zeros(1,length(posD.bigOut));
is806post = zeros(1,length(posD.bigOut));
is792post = zeros(1,length(posD.bigOut));
for n = 1 : length(posD.bigOut)
    posD.bigOut{n,15} = 0;
    posD.bigOut{n,16} = 0;
    if any(posD.bigOut{n,8}.respFreqZP(:,2)<0.05)
        % Neuron responsive
        posD.bigOut{n,15} = 1;
        sigResps = posD.bigOut{n,8}.respFreqZP(:,2)<0.05;
        zVals = posD.bigOut{n,8}.respFreqZP(sigResps,1);
        if any(zVals < 0) && any(zVals > 0) %Mixed responsive neuron, skip
            posD.bigOut{n,17} = 'mix';
            [~,IndxMax] = max(posD.bigOut{n,8}.respFreqZP(:,1)); %max resp
            [~,IndxMin] = min(posD.bigOut{n,8}.respFreqZP(:,1)); %min resp
            maxResps = posD.bigOut{n,8}.(strcat('f',num2str(freqList(IndxMax)))).avgRespPostStim;
            minResps = posD.bigOut{n,8}.(strcat('f',num2str(freqList(IndxMin)))).avgRespPostStim;
            [P,H] = ranksum(maxResps,minResps);
            if P < 0.05
                posD.bigOut{n,16} = 1;
            else 
                posD.bigOut{n,16} = -1;
            end
        elseif all(zVals > 0) %Positive responding neuron
            posD.bigOut{n,17} = 'pos';
            [~,IndxMax] = max(posD.bigOut{n,8}.respFreqZP(:,1)); %max resp
            [~,IndxMin] = min(posD.bigOut{n,8}.respFreqZP(:,1)); %min resp
            maxResps = posD.bigOut{n,8}.(strcat('f',num2str(freqList(IndxMax)))).avgRespPostStim;
            minResps = posD.bigOut{n,8}.(strcat('f',num2str(freqList(IndxMin)))).avgRespPostStim;
            [P,H] = ranksum(maxResps,minResps);
            if P < 0.05
                posD.bigOut{n,16} = 1;
            else 
                posD.bigOut{n,16} = -1;
            end
        elseif all(zVals < 0) %negatively responding neuron, flip vals
            posD.bigOut{n,17} = 'neg';
            [~,IndxMax] = max(abs(posD.bigOut{n,8}.respFreqZP(:,1))); %max resp with abs() because its negative
            [~,IndxMin] = min(abs(posD.bigOut{n,8}.respFreqZP(:,1))); %min resp with abs() because its negative
            maxResps = posD.bigOut{n,8}.(strcat('f',num2str(freqList(IndxMax)))).avgRespPostStim;
            minResps = posD.bigOut{n,8}.(strcat('f',num2str(freqList(IndxMin)))).avgRespPostStim;
            [P,H] = ranksum(maxResps,minResps);
            if P < 0.05
                posD.bigOut{n,16} = 1;
            else 
                posD.bigOut{n,16} = -1;
            end
        end
    else
        %neuron non-responsive
        posD.bigOut{n,15} = 0;
    end
    %check for ROI overlap (PV+ / SOM+)
    if strcmp(posD.bigOut{n,5},'Overlap')
        isPVpost(n) = 1;
    end
    %assign mouse numbers
    if strcmp(posD.bigOut{n,9},'m810')
        is810post(n) = 1;
    elseif strcmp(posD.bigOut{n,9},'m806')
        is806post(n) = 1;
    elseif strcmp(posD.bigOut{n,9},'m792')
        is792post(n) = 1;
    else
        error(strcat("Mouse number, on row ",num2str(n)," not recognized!"));
    end
end

isResponsivePostTrain = logical([posD.bigOut{:,15}]);
isSelectivePostTrain = logical([posD.bigOut{:,16}]);
is810post = logical(is810post);
is806post = logical(is806post);
is792post = logical(is792post);
isPVpost = logical(isPVpost);
isPVpost = logical(isPVpost);

%% NEW PEARSON CORRELATION SECTION (ADDED 2023-08-16)

% Load and calculate from Pre-training Data
sessNamesPre = unique([preD.bigOut{:,14}]);
% preallocate cell array of correlation matricies for each seesion
% seesion (col.1), r value (col.2) and p val (col.3)
corrMatrixPre = cell(length(sessNamesPre),2);
for ns = 1 : length(sessNamesPre) %number of unique sessions
    % preallocate cell array of correlation matricies for each seesion
    corrMatrixPre{ns,1} = sessNamesPre(ns); 
    BOid = find([preD.bigOut{:,14}] == sessNamesPre(ns)); %Big Out ID
    rhoMat = NaN(length(BOid),length(BOid)); % pearson r value matrix [ROI1 , ROI2]
    pvalMat = NaN(length(BOid),length(BOid)); % pearson p value matrix [ROI1 , ROI2]
    for nn = 1 : length(BOid) %number of first dim (corr matrix)
        for nj = 1 : nn % number of second dimension (corr matrix)
            i1 = BOid(nn); %index number
            i2 = BOid(nj); %index number
            response1 = [];
            response2 = [];
            for nf = 1 : length(freqList) % loop through number frequencies
                fn1 = strcat('f',num2str(freqList(nf))); %field name here 'f000'
                response1 = [response1; preD.bigOut{i1,8}.(fn1).avgRespPostStim];
                response2 = [response2; preD.bigOut{i2,8}.(fn1).avgRespPostStim];
            end
            [rhoMat(nn,nj),pvalMat(nn,nj)] = corr(response1,response2);
        end
    end
    corrMatrixPre{ns,2} = rhoMat;
    corrMatrixPre{ns,3} = pvalMat;
end

%--------------------------------------------------------------------------
% Load and calculate from Post-training Data
sessNamesPost = unique([posD.bigOut{:,14}]);
% preallocate cell array of correlation matricies for each seesion
% seesion (col.1), r value (col.2) and p val (col.3)
corrMatrixPost = cell(length(sessNamesPost),2);

for ns = 1 : length(sessNamesPost) %number of unique sessions
    % preallocate cell array of correlation matricies for each seesion
    corrMatrixPost{ns,1} = sessNamesPost(ns); 
    BOid = find([posD.bigOut{:,14}] == sessNamesPost(ns)); %Big Out ID
    rhoMat = NaN(length(BOid),length(BOid),length(freqList)); % pearson r value matrix [ROI1 , ROI2 , FREQ]
    pvalMat = NaN(length(BOid),length(BOid),length(freqList)); % pearson p value matrix [ROI1 , ROI2 , FREQ]
    for nn = 1 : length(BOid) %number of first dim (corr matrix)
        for nj = 1 : nn % number of second dimension (corr matrix)
            i1 = BOid(nn); %index number
            i2 = BOid(nj); %index number
            response1 = [];
            response2 = [];
            for nf = 1 : length(freqList) % loop through number frequencies
                fn1 = strcat('f',num2str(freqList(nf))); %field name here 'f000'
                response1 = [response1; posD.bigOut{i1,8}.(fn1).avgRespPostStim];
                response2 = [response2; posD.bigOut{i2,8}.(fn1).avgRespPostStim];
            end
            [rhoMat(nn,nj),pvalMat(nn,nj)] = corr(response1,response2);
        end
    end
    corrMatrixPost{ns,2} = rhoMat;
    corrMatrixPost{ns,3} = pvalMat;
end



%% CORRELATION SECTION (ADDED 2023-08-01)

% -----------------------------------------------------------------------------------------
% Load and calculate from Pre-training Data
sessNamesPre = unique([preD.bigOut{:,14}]);
corrMatrixPre = cell(length(sessNamesPre),2);
for ns = 1 : length(sessNamesPre) %number of sessions
    corrMatrixPre{ns,1} = sessNamesPre(ns);
    BOid = find([preD.bigOut{:,14}] == sessNamesPre(ns)); %Big Out ID
    corrMat = NaN(length(BOid),length(BOid),2); %cell,cell,[R P]
    for nn = 1 : length(BOid) %number of first dim (corr matrix)
        for nj = nn : length(BOid) % number of second dimension (corr matrix)
            i1 = BOid(nn); %index number
            i2 = BOid(nj); %index number
            %reconstruct complete df/f trace for all frequencies & reps
            rawTrace1 = [];
            rawTrace2 = [];
            for nf = 1 : length(freqList) %number frequencies
                fn1 = strcat('f',num2str(freqList(nf))); %field name here 'f000'
                rawTrace1 = [rawTrace1; preD.bigOut{i1,8}.(fn1).allTraces(:)];
                rawTrace2 = [rawTrace2; preD.bigOut{i2,8}.(fn1).allTraces(:)];
            end
            [RHO,PVAL] = corr([rawTrace1 rawTrace2]);
            corrMat(nj,nn,1) = RHO(1,2);
            corrMat(nj,nn,2) = PVAL(1,2);
        end
    end
    corrMatrixPre{ns,2} = corrMat;
end

% -----------------------------------------------------------------------------------------
% Load and calculate from Post-training Data
sessNamesPost = unique([posD.bigOut{:,14}]);
corrMatrixPost = cell(length(sessNamesPost),2);
for ns = 1 : length(sessNamesPost) %number of sessions
    corrMatrixPost{ns,1} = sessNamesPost(ns);
    BOid = find([posD.bigOut{:,14}] == sessNamesPost(ns));
    corrMat = NaN(length(BOid),length(BOid),2); %cell,cell,[R P]
    for nn = 1 : length(BOid) %number of first dim (corr matrix)
        for nj = nn : length(BOid) % number of second dimension (corr matrix)
            i1 = BOid(nn); %index number
            i2 = BOid(nj); %index number
            %reconstruct complete df/f trace for all frequencies & reps
            rawTrace1 = [];
            rawTrace2 = [];
            for nf = 1 : length(freqList) %number frequencies
                fn1 = strcat('f',num2str(freqList(nf))); %field name here 'f000'
                rawTrace1 = [rawTrace1; posD.bigOut{i1,8}.(fn1).allTraces(:)];
                rawTrace2 = [rawTrace2; posD.bigOut{i2,8}.(fn1).allTraces(:)];
            end
            [RHO,PVAL] = corr([rawTrace1 rawTrace2]);
            corrMat(nj,nn,1) = RHO(1,2);
            corrMat(nj,nn,2) = PVAL(1,2);
        end
    end
    corrMatrixPost{ns,2} = corrMat;
end

%% Compare all correlations of PV:PV, Non-PV:Non-PV, and PV:Non-PV

R_PV_PV = [];
R_NONPV_NONPV = [];
R_PV_NONPV = [];
for ns = 1 : length(sessNamesPre) %number of sessions
    isCurrSess = [preD.bigOut{:,14}] == sessNamesPre(ns);
    pvROIs = find(isCurrSess & isPVpre); %curr session and PV
    nonpvROIs = find(isCurrSess & ~isPVpre); %curr session and nonPV
    floorVal = min([pvROIs nonpvROIs]) - 1; %correction
    pvROIs = pvROIs - floorVal;
    nonpvROIs = nonpvROIs - floorVal;
    % PV : PV
    pvpvRVals1 = squeeze(mean(corrMatrixPre{ns,2}(pvROIs,pvROIs,:),3));
    pvpvRVals1(logical(eye(size(pvpvRVals1)))) = NaN; %remove diagonal correlations "1"
    R_PV_PV = [R_PV_PV pvpvRVals1(~isnan(pvpvRVals1))'];
    % NONPV : NONPV
    nonpvnonpvRVals1 = squeeze(mean(corrMatrixPre{ns,2}(nonpvROIs,nonpvROIs,1),3));
    nonpvnonpvRVals1(logical(eye(size(nonpvnonpvRVals1)))) = NaN; %remove diagonal correlations "1"
    R_NONPV_NONPV = [R_NONPV_NONPV nonpvnonpvRVals1(~isnan(nonpvnonpvRVals1))'];
    % PV : NONPV
    array1 = squeeze(mean(corrMatrixPre{ns,2}(pvROIs,nonpvROIs,1),3));
    array1 = array1(~isnan(array1))'; %remove NaNs
    if size(array1,1) > size(array1,2)
        array1 = array1';
    end
    R_PV_NONPV = [R_PV_NONPV array1];
    array2 = squeeze(mean(corrMatrixPre{ns,2}(nonpvROIs,pvROIs,1),3));
    array2 = array2(~isnan(array2))';
    if size(array2,1) > size(array2,2)
        array2 = array2';
    end
    R_PV_NONPV = [R_PV_NONPV array2];
end 

forplotting = {R_PV_PV',R_NONPV_NONPV',R_PV_NONPV'};
figure("Color","w");
violin(forplotting);
ax = gca;
ax.XTick = [1 2 3];
ax.XTickLabel = {'PV:PV', 'Non-PV:Non-PV', 'PV:Non-PV'};
ylabel('Pearsons r');
title('Pearson Correlations between Cell Identity Pre-Training');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

R_PV_PV = [];
R_NONPV_NONPV = [];
R_PV_NONPV = [];
for ns = 1 : length(sessNamesPost) %number of sessions
    isCurrSess = [posD.bigOut{:,14}] == sessNamesPost(ns);
    pvROIs = find(isCurrSess & isPVpost); %curr session and PV
    nonpvROIs = find(isCurrSess & ~isPVpost); %curr session and nonPV
    floorVal = min([pvROIs nonpvROIs]) - 1; %correction
    pvROIs = pvROIs - floorVal;
    nonpvROIs = nonpvROIs - floorVal;
    % PV : PV
    pvpvRVals1 = squeeze(mean(corrMatrixPost{ns,2}(pvROIs,pvROIs,1),3));
    pvpvRVals1(logical(eye(size(pvpvRVals1)))) = NaN; %remove diagonal correlations "1"
    R_PV_PV = [R_PV_PV pvpvRVals1(~isnan(pvpvRVals1))'];
    % NONPV : NONPV
    nonpvnonpvRVals1 = squeeze(mean(corrMatrixPost{ns,2}(nonpvROIs,nonpvROIs,1),3));
    nonpvnonpvRVals1(logical(eye(size(nonpvnonpvRVals1)))) = NaN; %remove diagonal correlations "1"
    R_NONPV_NONPV = [R_NONPV_NONPV nonpvnonpvRVals1(~isnan(nonpvnonpvRVals1))'];
    % PV : NONPV
    array1 = squeeze(mean(corrMatrixPost{ns,2}(pvROIs,nonpvROIs,1),3));
    array1 = array1(~isnan(array1))'; %remove NaNs
    if size(array1,1) > size(array1,2)
        array1 = array1';
    end
    R_PV_NONPV = [R_PV_NONPV array1];
    array2 = squeeze(mean(corrMatrixPost{ns,2}(nonpvROIs,pvROIs,1),3));
    array2 = array2(~isnan(array2))';
    if size(array2,1) > size(array2,2)
        array2 = array2';
    end
    R_PV_NONPV = [R_PV_NONPV array2];
end 

forplotting = {R_PV_PV',R_NONPV_NONPV',R_PV_NONPV'};
figure("Color","w");
violin(forplotting);
ax = gca;
ax.XTick = [1 2 3];
ax.XTickLabel = {'PV:PV', 'Non-PV:Non-PV', 'PV:Non-PV'};
ylabel('Pearsons r');
title('Pearson Correlations between Cell Identity Post-Training');


%% Compare correlations of selective PV:PV, Non-PV:Non-PV, and PV:Non-PV (pre training)

R_PV_PV1 = [];
R_PV_PV2 = [];
R_NONPV_NONPV1 = [];
R_NONPV_NONPV2 = [];
R_PV_NONPV1 = [];
R_PV_NONPV2 = [];
for ns = 1 : length(sessNamesPre) %number of sessions
    isCurrSess = [preD.bigOut{:,14}] == sessNamesPre(ns);

    pvROIs1 = find(isCurrSess & isPVpre & ~isSelectivePreTrain); %curr session and PV, non selective
    nonpvROIs1 = find(isCurrSess & ~isPVpre & ~isSelectivePreTrain); %curr session and nonPV, non selective
    pvROIs2 = find(isCurrSess & isPVpre & isSelectivePreTrain); %curr session and PV, selective
    nonpvROIs2 = find(isCurrSess & ~isPVpre & isSelectivePreTrain); %curr session and nonPV, selective

    floorVal = min([pvROIs1 nonpvROIs1 pvROIs2 nonpvROIs2]) - 1; %correction
    pvROIs1 = pvROIs1 - floorVal;
    nonpvROIs1 = nonpvROIs1 - floorVal;
    pvROIs2 = pvROIs2 - floorVal;
    nonpvROIs2 = nonpvROIs2 - floorVal;

    % PV : PV
    pvpvRVals1 = corrMatrixPre{ns,2}(pvROIs1,pvROIs1);
    pvpvRVals1(logical(eye(size(pvpvRVals1)))) = NaN; %remove diagonal correlations "1"
    R_PV_PV1 = [R_PV_PV1 pvpvRVals1(~isnan(pvpvRVals1))'];

    pvpvRVals2 = corrMatrixPre{ns,2}(pvROIs2,pvROIs2);
    pvpvRVals2(logical(eye(size(pvpvRVals2)))) = NaN; %remove diagonal correlations "1"
    R_PV_PV2 = [R_PV_PV2 pvpvRVals2(~isnan(pvpvRVals2))'];

    % NONPV : NONPV
    nonpvnonpvRVals1 = corrMatrixPre{ns,2}(nonpvROIs1,nonpvROIs1);
    nonpvnonpvRVals1(logical(eye(size(nonpvnonpvRVals1)))) = NaN; %remove diagonal correlations "1"
    R_NONPV_NONPV1 = [R_NONPV_NONPV1 nonpvnonpvRVals1(~isnan(nonpvnonpvRVals1))'];

    nonpvnonpvRVals2 = corrMatrixPre{ns,2}(nonpvROIs2,nonpvROIs2);
    nonpvnonpvRVals2(logical(eye(size(nonpvnonpvRVals2)))) = NaN; %remove diagonal correlations "1"
    R_NONPV_NONPV2 = [R_NONPV_NONPV2 nonpvnonpvRVals2(~isnan(nonpvnonpvRVals2))'];

    % PV : NONPV
    array1 = corrMatrixPre{ns,2}(pvROIs1,nonpvROIs1);
    array1 = array1(~isnan(array1))'; %remove NaNs
    if size(array1,1) > size(array1,2)
        array1 = array1';
    end
    R_PV_NONPV1 = [R_PV_NONPV1 array1];
    array2 = corrMatrixPre{ns,2}(nonpvROIs1,pvROIs1);
    array2 = array2(~isnan(array2))';
    if size(array2,1) > size(array2,2)
        array2 = array2';
    end
    R_PV_NONPV1 = [R_PV_NONPV1 array2];

    array3 = corrMatrixPre{ns,2}(pvROIs2,nonpvROIs2);
    array3 = array3(~isnan(array3))'; %remove NaNs
    if size(array3,1) > size(array3,2)
        array3 = array3';
    end
    R_PV_NONPV2 = [R_PV_NONPV2 array3];
    array4 = corrMatrixPre{ns,2}(nonpvROIs2,pvROIs2);
    array4 = array4(~isnan(array4))';
    if size(array4,1) > size(array4,2)
        array4 = array4';
    end
    R_PV_NONPV2 = [R_PV_NONPV2 array4];
end 

forplotting = {R_PV_PV1',R_PV_PV2',  R_NONPV_NONPV1',R_NONPV_NONPV2',  R_PV_NONPV1',R_PV_NONPV2'};
figure("Color","w");
colorMat = [1 0.5 0; 0.5 1 0; 1 0.5 0; 0.5 1 0; 1 0.5 0; 0.5 1 0];
violin(forplotting,'facecolor',colorMat);
ax = gca;
ax.XTick = [1.5 3.5 5.5];
ax.XTickLabel = {'PV:PV', 'Non-PV:Non-PV','PV:Non-PV'};
ylabel('Pearsons r');
title('Pearson Correlations between Cell Identity & Selectivity (Pre-Training)');
legend({'non-selective','mean','median','selective'});

%%  Compare correlations of selective PV:PV, Non-PV:Non-PV, and PV:Non-PV (post training)

R_PV_PV1 = [];
R_PV_PV2 = [];
R_NONPV_NONPV1 = [];
R_NONPV_NONPV2 = [];
R_PV_NONPV1 = [];
R_PV_NONPV2 = [];
for ns = 1 : length(sessNamesPost) %number of sessions
    isCurrSess = [posD.bigOut{:,14}] == sessNamesPost(ns);

    pvROIs1 = find(isCurrSess & isPVpost & ~isSelectivePostTrain); %curr session and PV, non selective
    nonpvROIs1 = find(isCurrSess & ~isPVpost & ~isSelectivePostTrain); %curr session and nonPV, non selective
    pvROIs2 = find(isCurrSess & isPVpost & isSelectivePostTrain); %curr session and PV, selective
    nonpvROIs2 = find(isCurrSess & ~isPVpost & isSelectivePostTrain); %curr session and nonPV, selective

    floorVal = min([pvROIs1 nonpvROIs1 pvROIs2 nonpvROIs2]) - 1; %correction
    pvROIs1 = pvROIs1 - floorVal;
    nonpvROIs1 = nonpvROIs1 - floorVal;
    pvROIs2 = pvROIs2 - floorVal;
    nonpvROIs2 = nonpvROIs2 - floorVal;

    % PV : PV
    pvpvRVals1 = corrMatrixPost{ns,2}(pvROIs1,pvROIs1);
    pvpvRVals1(logical(eye(size(pvpvRVals1)))) = NaN; %remove diagonal correlations "1"
    R_PV_PV1 = [R_PV_PV1 pvpvRVals1(~isnan(pvpvRVals1))'];

    pvpvRVals2 = corrMatrixPost{ns,2}(pvROIs2,pvROIs2);
    pvpvRVals2(logical(eye(size(pvpvRVals2)))) = NaN; %remove diagonal correlations "1"
    R_PV_PV2 = [R_PV_PV2 pvpvRVals2(~isnan(pvpvRVals2))'];

    % NONPV : NONPV
    nonpvnonpvRVals1 = corrMatrixPost{ns,2}(nonpvROIs1,nonpvROIs1);
    nonpvnonpvRVals1(logical(eye(size(nonpvnonpvRVals1)))) = NaN; %remove diagonal correlations "1"
    R_NONPV_NONPV1 = [R_NONPV_NONPV1 nonpvnonpvRVals1(~isnan(nonpvnonpvRVals1))'];

    nonpvnonpvRVals2 = corrMatrixPost{ns,2}(nonpvROIs2,nonpvROIs2);
    nonpvnonpvRVals2(logical(eye(size(nonpvnonpvRVals2)))) = NaN; %remove diagonal correlations "1"
    R_NONPV_NONPV2 = [R_NONPV_NONPV2 nonpvnonpvRVals2(~isnan(nonpvnonpvRVals2))'];

    % PV : NONPV
    array1 = corrMatrixPost{ns,2}(pvROIs1,nonpvROIs1);
    array1 = array1(~isnan(array1))'; %remove NaNs
    if size(array1,1) > size(array1,2)
        array1 = array1';
    end
    R_PV_NONPV1 = [R_PV_NONPV1 array1];
    array2 = corrMatrixPost{ns,2}(nonpvROIs1,pvROIs1);
    array2 = array2(~isnan(array2))';
    if size(array2,1) > size(array2,2)
        array2 = array2';
    end
    R_PV_NONPV1 = [R_PV_NONPV1 array2];

    array3 = corrMatrixPost{ns,2}(pvROIs2,nonpvROIs2);
    array3 = array3(~isnan(array3))'; %remove NaNs
    if size(array3,1) > size(array3,2)
        array3 = array3';
    end
    R_PV_NONPV2 = [R_PV_NONPV2 array3];
    array4 = corrMatrixPost{ns,2}(nonpvROIs2,pvROIs2);
    array4 = array4(~isnan(array4))';
    if size(array4,1) > size(array4,2)
        array4 = array4';
    end
    R_PV_NONPV2 = [R_PV_NONPV2 array4];
end 

forplotting = {R_PV_PV1',R_PV_PV2',  R_NONPV_NONPV1',R_NONPV_NONPV2',  R_PV_NONPV1',R_PV_NONPV2'};
figure("Color","w");
colorMat = [1 0.5 0; 0.5 1 0; 1 0.5 0; 0.5 1 0; 1 0.5 0; 0.5 1 0];
violin(forplotting,'facecolor',colorMat);
ax = gca;
ax.XTick = [1.5 3.5 5.5];
ax.XTickLabel = {'PV:PV', 'Non-PV:Non-PV','PV:Non-PV'};
ylabel('Pearsons r');
title('Pearson Correlations between Cell Identity & Selectivity (Post-Training)');
legend({'non-selective','mean','median','selective'});

%% PLOTTING SECTION ALL NEURONS FIG 1.0

nTotPre = length(isResponsivePreTrain);
nRespPre = sum(isResponsivePreTrain);
nNonRespPre = length(isResponsivePreTrain) - nRespPre;
nSelPre = sum(isSelectivePreTrain == 1);
nNonSelPre = sum(isSelectivePreTrain == -1);

nTotPost = length(isResponsivePostTrain);
nRespPost = sum(isResponsivePostTrain);
nNonRespPost = length(isResponsivePostTrain) - nRespPost;
nSelPost = sum(isSelectivePostTrain == 1);
nNonSelPost = sum(isSelectivePostTrain == -1);


f0 = figure;  hold on;
% x0 = [1:4]
% y0 = [nRespPre/nTotPre nNonRespPre/nTotPre; nSelPre/nTotPre nNonSelPre/nTotPre;
%     nRespPost/nTotPost nNonRespPost/nTotPost; nSelPost/nTotPost nNonSelPost/nTotPost]';
% b = bar(x0,y0);

b0 = bar(1,[nRespPre/nTotPre nNonRespPre/nTotPre],'grouped');
b0(1,1).FaceColor = [0 0.4470 0.7410];
b0(1,2).FaceColor = [0 0.4470 0.7410] .* 0.5;
stdVal0 = mannyPropSTD(nRespPre,nTotPre);
text(b0(1,1).XEndPoints,b0(1,1).YEndPoints+stdVal0,string(strcat("n=",num2str(nRespPre))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b0(1,2).XEndPoints,b0(1,2).YEndPoints+stdVal0,string(strcat("n=",num2str(nNonRespPre))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

b1 = bar(2,[nSelPre/nTotPre nNonSelPre/nTotPre],'grouped');
b1(1,1).FaceColor = [0.4660 0.6740 0.1880];
b1(1,2).FaceColor = [0.4660 0.6740 0.1880] .* 0.5;
stdVal1 = mannyPropSTD(nSelPre,nTotPre);
text(b1(1,1).XEndPoints,b1(1,1).YEndPoints+stdVal1,string(strcat("n=",num2str(nSelPre))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b1(1,2).XEndPoints,b1(1,2).YEndPoints+stdVal1,string(strcat("n=",num2str(nNonSelPre))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

b2 = bar(3,[nRespPost/nTotPost nNonRespPost/nTotPost],'grouped');
b2(1,1).FaceColor = [0 0.4470 0.7410];
b2(1,2).FaceColor = [0 0.4470 0.7410] .* 0.5;
stdVal2 = mannyPropSTD(nRespPost,nTotPost);
text(b2(1,1).XEndPoints,b2(1,1).YEndPoints+stdVal2,string(strcat("n=",num2str(nRespPost))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b2(1,2).XEndPoints,b2(1,2).YEndPoints+stdVal2,string(strcat("n=",num2str(nNonRespPost))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

b3 = bar(4,[nSelPost/nTotPost nNonSelPost/nTotPost],'grouped');
b3(1,1).FaceColor = [0.4660 0.6740 0.1880];
b3(1,2).FaceColor = [0.4660 0.6740 0.1880] .* 0.5;
stdVal3 = mannyPropSTD(nSelPost,nTotPost);
text(b3(1,1).XEndPoints,b3(1,1).YEndPoints+stdVal3,string(strcat("n=",num2str(nSelPost))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b3(1,2).XEndPoints,b3(1,2).YEndPoints+stdVal3,string(strcat("n=",num2str(nNonSelPost))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

errorbar(b0(1,1).XEndPoints,b0(1,1).YEndPoints,mannyPropSTD(nRespPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b0(1,2).XEndPoints,b0(1,2).YEndPoints,mannyPropSTD(nNonRespPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b1(1,1).XEndPoints,b1(1,1).YEndPoints,mannyPropSTD(nSelPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b1(1,2).XEndPoints,b1(1,2).YEndPoints,mannyPropSTD(nNonSelPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b2(1,1).XEndPoints,b2(1,1).YEndPoints,mannyPropSTD(nRespPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b2(1,2).XEndPoints,b2(1,2).YEndPoints,mannyPropSTD(nNonRespPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b3(1,1).XEndPoints,b3(1,1).YEndPoints,mannyPropSTD(nSelPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b3(1,2).XEndPoints,b3(1,2).YEndPoints,mannyPropSTD(nNonSelPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);

ax0 = gca;
ax0.XTick = [1.5 3.4];
ax0.XTickLabel = {'Pre-training';'Post-training'};
ax0.YLim = [0 1];
ylabel('Percentage of Neurons')
title('Catagorization of Neurons ');
legend({'Responsive','Non-Responsive','Selective','Non-Selective'});


%% PLOTTING SECTION , POOL ALL PRE- and POST- Data, split by ID

nTotPV = sum([isPVpre isPVpost]);
nRespPV = sum(isPVpre & isResponsivePreTrain) + sum(isPVpost & isResponsivePostTrain);
nNonRespPV = nTotPV - nRespPV;
iSELpre = isSelectivePreTrain == 1;
iSELpost = isSelectivePostTrain == 1;
nSelPV = sum(iSELpre & isPVpre) + sum(iSELpost & isPVpost);
nNonSelPV = nRespPV - nSelPV;

nTotNPV = sum([~isPVpre ~isPVpost]);
nRespNPV = sum(~isPVpre & isResponsivePreTrain) + sum(~isPVpost & isResponsivePostTrain);
nNonRespNPV = nTotNPV - nRespNPV;
nSelNPV = sum(iSELpre & ~isPVpre) + sum(iSELpost & ~isPVpost);
nNonSelNPV = nRespNPV - nSelNPV;


f0 = figure("Color","w");  hold on;
% x0 = [1:4]
% y0 = [nRespPV/nTotPV nNonRespPV/nTotPV; nSelPV/nTotPV nNonSelPV/nTotPV;
%     nRespNPV/nTotNPV nNonRespNPV/nTotNPV; nSelNPV/nTotNPV nNonSelNPV/nTotNPV]';
% b = bar(x0,y0);

b0 = bar(1,[nRespPV/nTotPV nNonRespPV/nTotPV],'grouped');
b0(1,1).FaceColor = [0 0.4470 0.7410];
b0(1,2).FaceColor = [0 0.4470 0.7410] .* 0.5;
stdVal0 = mannyPropSTD(nRespPV,nTotPV);
text(b0(1,1).XEndPoints,b0(1,1).YEndPoints+stdVal0,string(strcat("n=",num2str(nRespPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b0(1,2).XEndPoints,b0(1,2).YEndPoints+stdVal0,string(strcat("n=",num2str(nNonRespPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

b1 = bar(2,[nSelPV/nTotPV nNonSelPV/nTotPV],'grouped');
b1(1,1).FaceColor = [0.4660 0.6740 0.1880];
b1(1,2).FaceColor = [0.4660 0.6740 0.1880] .* 0.5;
stdVal1 = mannyPropSTD(nSelPV,nTotPV);
text(b1(1,1).XEndPoints,b1(1,1).YEndPoints+stdVal1,string(strcat("n=",num2str(nSelPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b1(1,2).XEndPoints,b1(1,2).YEndPoints+stdVal1,string(strcat("n=",num2str(nNonSelPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

b2 = bar(3,[nRespNPV/nTotNPV nNonRespNPV/nTotNPV],'grouped');
b2(1,1).FaceColor = [0 0.4470 0.7410];
b2(1,2).FaceColor = [0 0.4470 0.7410] .* 0.5;
stdVal2 = mannyPropSTD(nRespNPV,nTotNPV);
text(b2(1,1).XEndPoints,b2(1,1).YEndPoints+stdVal2,string(strcat("n=",num2str(nRespNPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b2(1,2).XEndPoints,b2(1,2).YEndPoints+stdVal2,string(strcat("n=",num2str(nNonRespNPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

b3 = bar(4,[nSelNPV/nTotNPV nNonSelNPV/nTotNPV],'grouped');
b3(1,1).FaceColor = [0.4660 0.6740 0.1880];
b3(1,2).FaceColor = [0.4660 0.6740 0.1880] .* 0.5;
stdVal3 = mannyPropSTD(nSelNPV,nTotNPV);
text(b3(1,1).XEndPoints,b3(1,1).YEndPoints+stdVal3,string(strcat("n=",num2str(nSelNPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b3(1,2).XEndPoints,b3(1,2).YEndPoints+stdVal3,string(strcat("n=",num2str(nNonSelNPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

errorbar(b0(1,1).XEndPoints,b0(1,1).YEndPoints,mannyPropSTD(nRespPV,nTotPV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b0(1,2).XEndPoints,b0(1,2).YEndPoints,mannyPropSTD(nNonRespPV,nTotPV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b1(1,1).XEndPoints,b1(1,1).YEndPoints,mannyPropSTD(nSelPV,nTotPV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b1(1,2).XEndPoints,b1(1,2).YEndPoints,mannyPropSTD(nNonSelPV,nTotPV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b2(1,1).XEndPoints,b2(1,1).YEndPoints,mannyPropSTD(nRespPV,nTotPV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b2(1,2).XEndPoints,b2(1,2).YEndPoints,mannyPropSTD(nNonRespPV,nTotPV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b3(1,1).XEndPoints,b3(1,1).YEndPoints,mannyPropSTD(nSelPV,nTotPV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b3(1,2).XEndPoints,b3(1,2).YEndPoints,mannyPropSTD(nNonSelPV,nTotPV),...
    Color=[0,0,0],LineWidth=1.4);

ax0 = gca;
ax0.XTick = [1.5 3.4];
ax0.XTickLabel = {'PV';'Non-PV'};
ax0.YLim = [0 1];
ylabel('Percentage of Neurons')
title('Catagorization of Neural Responses by Cell Type ');
legend({'Responsive','Non-Responsive','Selective','Non-Selective'});

%% Comparrison between cell identity (subplotted) and training stage vs. responsivity and selectivity
% PLOTTING SECTION of PV-ONLY FIG 1.1

nTotPre = sum(isPVpre);
nRespPre = sum(isResponsivePreTrain & isPVpre);
nNonRespPre = nTotPre - nRespPre;
nSelPre = sum(isSelectivePreTrain == 1 & isPVpre);
nNonSelPre = sum(isSelectivePreTrain == -1 & isPVpre);

nTotPost = sum(isPVpost);
nRespPost = sum(isResponsivePostTrain & isPVpost);
nNonRespPost = nTotPost - nRespPost;
nSelPost = sum(isSelectivePostTrain == 1 & isPVpost);
nNonSelPost = sum(isSelectivePostTrain == -1 & isPVpost);

f1 = figure("Color","w"); hold on
subplot(1,2,1); hold on;
% x0 = [1:4]
% y0 = [nRespPre/nTotPre nNonRespPre/nTotPre; nSelPre/nTotPre nNonSelPre/nTotPre;
%     nRespPost/nTotPost nNonRespPost/nTotPost; nSelPost/nTotPost nNonSelPost/nTotPost]';
% b = bar(x0,y0);

b0 = bar(1,[nRespPre/nTotPre nNonRespPre/nTotPre],'grouped');
b0(1,1).FaceColor = [0 0.4470 0.7410];
b0(1,2).FaceColor = [0 0.4470 0.7410] .* 0.6;
stdVal0 = mannyPropSTD(nRespPre,nTotPre);
text(b0(1,1).XEndPoints,b0(1,1).YEndPoints+stdVal0,string(strcat("n=",num2str(nRespPre))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b0(1,2).XEndPoints,b0(1,2).YEndPoints+stdVal0,string(strcat("n=",num2str(nNonRespPre))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

b1 = bar(2,[nSelPre/nRespPre nNonSelPre/nRespPre],'grouped');
b1(1,1).FaceColor = [0.4660 0.6740 0.1880];
b1(1,2).FaceColor = [0.4660 0.6740 0.1880] .* 0.6;
stdVal1 = mannyPropSTD(nSelPre,nRespPre);
text(b1(1,1).XEndPoints,b1(1,1).YEndPoints+stdVal1,string(strcat("n=",num2str(nSelPre))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b1(1,2).XEndPoints,b1(1,2).YEndPoints+stdVal1,string(strcat("n=",num2str(nNonSelPre))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

b2 = bar(3,[nRespPost/nTotPost nNonRespPost/nTotPost],'grouped');
b2(1,1).FaceColor = [0 0.4470 0.7410];
b2(1,2).FaceColor = [0 0.4470 0.7410] .* 0.6;
stdVal2 = mannyPropSTD(nRespPost,nTotPost);
text(b2(1,1).XEndPoints,b2(1,1).YEndPoints+stdVal2,string(strcat("n=",num2str(nRespPost))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b2(1,2).XEndPoints,b2(1,2).YEndPoints+stdVal2,string(strcat("n=",num2str(nNonRespPost))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

b3 = bar(4,[nSelPost/nRespPost nNonSelPost/nRespPost],'grouped');
b3(1,1).FaceColor = [0.4660 0.6740 0.1880];
b3(1,2).FaceColor = [0.4660 0.6740 0.1880] .* 0.6;
stdVal3 = mannyPropSTD(nSelPost,nRespPost);
text(b3(1,1).XEndPoints,b3(1,1).YEndPoints+stdVal3,string(strcat("n=",num2str(nSelPost))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b3(1,2).XEndPoints,b3(1,2).YEndPoints+stdVal3,string(strcat("n=",num2str(nNonSelPost))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

errorbar(b0(1,1).XEndPoints,b0(1,1).YEndPoints,mannyPropSTD(nRespPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b0(1,2).XEndPoints,b0(1,2).YEndPoints,mannyPropSTD(nNonRespPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b1(1,1).XEndPoints,b1(1,1).YEndPoints,mannyPropSTD(nSelPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b1(1,2).XEndPoints,b1(1,2).YEndPoints,mannyPropSTD(nNonSelPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b2(1,1).XEndPoints,b2(1,1).YEndPoints,mannyPropSTD(nRespPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b2(1,2).XEndPoints,b2(1,2).YEndPoints,mannyPropSTD(nNonRespPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b3(1,1).XEndPoints,b3(1,1).YEndPoints,mannyPropSTD(nSelPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b3(1,2).XEndPoints,b3(1,2).YEndPoints,mannyPropSTD(nNonSelPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);

ax0 = gca;
ax0.XTick = [1.5 3.4];
ax0.XTickLabel = {'Pre-training';'Post-training'};
ax0.YLim = [0 1];
ylabel('Percentage of Neurons')
title('Categorization of PV+ Neurons ');
legend({'Responsive','Non-Responsive','Selective','Non-Selective'});



% PLOTTING SECTION of NON PV-ONLY FIG 1.2

nTotPre = sum(~isPVpre);
nRespPre = sum(isResponsivePreTrain & ~isPVpre);
nNonRespPre = nTotPre - nRespPre;
nSelPre = sum(isSelectivePreTrain == 1 & ~isPVpre);
nNonSelPre = sum(isSelectivePreTrain == -1 & ~isPVpre);

nTotPost = sum(~isPVpost);
nRespPost = sum(isResponsivePostTrain & ~isPVpost);
nNonRespPost = nTotPost - nRespPost;
nSelPost = sum(isSelectivePostTrain == 1 & ~isPVpost);
nNonSelPost = sum(isSelectivePostTrain == -1 & ~isPVpost);


subplot(1,2,2); hold on;
% x0 = [1:4]
% y0 = [nRespPre/nTotPre nNonRespPre/nTotPre; nSelPre/nTotPre nNonSelPre/nTotPre;
%     nRespPost/nTotPost nNonRespPost/nTotPost; nSelPost/nTotPost nNonSelPost/nTotPost]';
% b = bar(x0,y0);

b0 = bar(1,[nRespPre/nTotPre nNonRespPre/nTotPre],'grouped');
b0(1,1).FaceColor = [0 0.4470 0.7410];
b0(1,2).FaceColor = [0 0.4470 0.7410] .* 0.6;
stdVal0 = mannyPropSTD(nRespPre,nTotPre);
text(b0(1,1).XEndPoints,b0(1,1).YEndPoints+stdVal0,string(strcat("n=",num2str(nRespPre))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b0(1,2).XEndPoints,b0(1,2).YEndPoints+stdVal0,string(strcat("n=",num2str(nNonRespPre))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

b1 = bar(2,[nSelPre/nRespPre nNonSelPre/nRespPre],'grouped');
b1(1,1).FaceColor = [0.4660 0.6740 0.1880];
b1(1,2).FaceColor = [0.4660 0.6740 0.1880] .* 0.6;
stdVal1 = mannyPropSTD(nSelPre,nRespPre);
text(b1(1,1).XEndPoints,b1(1,1).YEndPoints+stdVal1,string(strcat("n=",num2str(nSelPre))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b1(1,2).XEndPoints,b1(1,2).YEndPoints+stdVal1,string(strcat("n=",num2str(nNonSelPre))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

b2 = bar(3,[nRespPost/nTotPost nNonRespPost/nTotPost],'grouped');
b2(1,1).FaceColor = [0 0.4470 0.7410];
b2(1,2).FaceColor = [0 0.4470 0.7410] .* 0.6;
stdVal2 = mannyPropSTD(nRespPost,nTotPost);
text(b2(1,1).XEndPoints,b2(1,1).YEndPoints+stdVal2,string(strcat("n=",num2str(nRespPost))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b2(1,2).XEndPoints,b2(1,2).YEndPoints+stdVal2,string(strcat("n=",num2str(nNonRespPost))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

b3 = bar(4,[nSelPost/nRespPost nNonSelPost/nRespPost],'grouped');
b3(1,1).FaceColor = [0.4660 0.6740 0.1880];
b3(1,2).FaceColor = [0.4660 0.6740 0.1880] .* 0.6;
stdVal3 = mannyPropSTD(nSelPost,nRespPost);
text(b3(1,1).XEndPoints,b3(1,1).YEndPoints+stdVal3,string(strcat("n=",num2str(nSelPost))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b3(1,2).XEndPoints,b3(1,2).YEndPoints+stdVal3,string(strcat("n=",num2str(nNonSelPost))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

errorbar(b0(1,1).XEndPoints,b0(1,1).YEndPoints,mannyPropSTD(nRespPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b0(1,2).XEndPoints,b0(1,2).YEndPoints,mannyPropSTD(nNonRespPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b1(1,1).XEndPoints,b1(1,1).YEndPoints,mannyPropSTD(nSelPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b1(1,2).XEndPoints,b1(1,2).YEndPoints,mannyPropSTD(nNonSelPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b2(1,1).XEndPoints,b2(1,1).YEndPoints,mannyPropSTD(nRespPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b2(1,2).XEndPoints,b2(1,2).YEndPoints,mannyPropSTD(nNonRespPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b3(1,1).XEndPoints,b3(1,1).YEndPoints,mannyPropSTD(nSelPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b3(1,2).XEndPoints,b3(1,2).YEndPoints,mannyPropSTD(nNonSelPre,nTotPre),...
    Color=[0,0,0],LineWidth=1.4);

ax0 = gca;
ax0.XTick = [1.5 3.4];
ax0.XTickLabel = {'Pre-training';'Post-training'};
ax0.YLim = [0 1];
ylabel('Percentage of Neurons')
title('Categorization of Non-PV Neurons ');
legend({'Responsive','Non-Responsive','Selective','Non-Selective'});


%% Comparrison between training stage (subplotted) and cell identity vs. responsivity and selectivity (2023/08/15)

% DEFINE LOGICAL INDEXES TO SPLIT DATA
nTotPrePV = sum(isPVpre);
nRespPrePV = sum(isResponsivePreTrain & isPVpre);
nNonRespPrePV = nTotPrePV - nRespPrePV;
nSelPrePV = sum(isSelectivePreTrain == 1 & isPVpre);
nNonSelPrePV = sum(isSelectivePreTrain == -1 & isPVpre);
nTotPostPV = sum(isPVpost);
nRespPostPV = sum(isResponsivePostTrain & isPVpost);
nNonRespPostPV = nTotPostPV - nRespPostPV;
nSelPostPV = sum(isSelectivePostTrain == 1 & isPVpost);
nNonSelPostPV = sum(isSelectivePostTrain == -1 & isPVpost);
nTotPreNPV = sum(~isPVpre);
nRespPreNPV = sum(isResponsivePreTrain & ~isPVpre);
nNonRespPreNPV = nTotPreNPV - nRespPreNPV;
nSelPreNPV = sum(isSelectivePreTrain == 1 & ~isPVpre);
nNonSelPreNPV = sum(isSelectivePreTrain == -1 & ~isPVpre);
nTotPostNPV = sum(~isPVpost);
nRespPostNPV = sum(isResponsivePostTrain & ~isPVpost);
nNonRespPostNPV = nTotPostNPV - nRespPostNPV;
nSelPostNPV = sum(isSelectivePostTrain == 1 & ~isPVpost);
nNonSelPostNPV = sum(isSelectivePostTrain == -1 & ~isPVpost);

% ---------------------------------------------------------------------------------------------------
% PLOTTING SECTION of PRE TRAINING ONLY
f1 = figure("Color","w"); subplot(1,2,1); hold on;


% PV | PRE-T | RESP_NONRESP
b0 = bar(1,[nRespPrePV/nTotPrePV nNonRespPrePV/nTotPrePV],'grouped');
b0(1,1).FaceColor = [0 0.4470 0.7410];
b0(1,2).FaceColor = [0 0.4470 0.7410] .* 0.6;
stdVal0 = mannyPropSTD(nRespPrePV,nTotPrePV);
text(b0(1,1).XEndPoints,b0(1,1).YEndPoints+stdVal0,string(strcat("n=",num2str(nRespPrePV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b0(1,2).XEndPoints,b0(1,2).YEndPoints+stdVal0,string(strcat("n=",num2str(nNonRespPrePV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

% PV | PRE-T | SEL_NONSEL
b1 = bar(2,[nSelPrePV/nRespPrePV nNonSelPrePV/nRespPrePV],'grouped');
b1(1,1).FaceColor = [0.4660 0.6740 0.1880];
b1(1,2).FaceColor = [0.4660 0.6740 0.1880] .* 0.6;
stdVal1 = mannyPropSTD(nSelPrePV,nRespPrePV);
text(b1(1,1).XEndPoints,b1(1,1).YEndPoints+stdVal1,string(strcat("n=",num2str(nSelPrePV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b1(1,2).XEndPoints,b1(1,2).YEndPoints+stdVal1,string(strcat("n=",num2str(nNonSelPrePV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

% PV | PRE-T | ERRORBARS (formatted this way for legend accuracy)
errorbar(b0(1,1).XEndPoints,b0(1,1).YEndPoints,mannyPropSTD(nRespPrePV,nTotPrePV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b0(1,2).XEndPoints,b0(1,2).YEndPoints,mannyPropSTD(nNonRespPrePV,nTotPrePV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b1(1,1).XEndPoints,b1(1,1).YEndPoints,mannyPropSTD(nSelPrePV,nTotPrePV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b1(1,2).XEndPoints,b1(1,2).YEndPoints,mannyPropSTD(nNonSelPrePV,nTotPrePV),...
    Color=[0,0,0],LineWidth=1.4);

% NON-PV | PRE-T | RESP_NONRESP
b2 = bar(3,[nRespPreNPV/nTotPreNPV nNonRespPreNPV/nTotPreNPV],'grouped');
b2(1,1).FaceColor = [0 0.4470 0.7410];
b2(1,2).FaceColor = [0 0.4470 0.7410] .* 0.6;
stdVal0 = mannyPropSTD(nRespPreNPV,nTotPreNPV);
text(b2(1,1).XEndPoints,b2(1,1).YEndPoints+stdVal0,string(strcat("n=",num2str(nRespPreNPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b2(1,2).XEndPoints,b2(1,2).YEndPoints+stdVal0,string(strcat("n=",num2str(nNonRespPreNPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
errorbar(b2(1,1).XEndPoints,b2(1,1).YEndPoints,mannyPropSTD(nRespPreNPV,nTotPreNPV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b2(1,2).XEndPoints,b2(1,2).YEndPoints,mannyPropSTD(nNonRespPreNPV,nTotPreNPV),...
    Color=[0,0,0],LineWidth=1.4);

% NON-PV | PRE-T | SEL_NONSEL
b3 = bar(4,[nSelPreNPV/nRespPreNPV nNonSelPreNPV/nRespPreNPV],'grouped');
b3(1,1).FaceColor = [0.4660 0.6740 0.1880];
b3(1,2).FaceColor = [0.4660 0.6740 0.1880] .* 0.6;
stdVal1 = mannyPropSTD(nSelPreNPV,nRespPreNPV);
text(b3(1,1).XEndPoints,b3(1,1).YEndPoints+stdVal1,string(strcat("n=",num2str(nSelPreNPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b3(1,2).XEndPoints,b3(1,2).YEndPoints+stdVal1,string(strcat("n=",num2str(nNonSelPreNPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
errorbar(b3(1,1).XEndPoints,b3(1,1).YEndPoints,mannyPropSTD(nSelPreNPV,nTotPreNPV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b3(1,2).XEndPoints,b3(1,2).YEndPoints,mannyPropSTD(nNonSelPreNPV,nTotPreNPV),...
    Color=[0,0,0],LineWidth=1.4);

% TITLE / LABELS / LEGEND
ax0 = gca;
ax0.XTick = [1.5 3.4];
ax0.XTickLabel = {'PV';'Non-PV'};
ax0.YLim = [0 1];
ylabel('Percentage of Neurons')
title('Categorization of Neurons Pre-Training');
legend({'Responsive','Non-Responsive','Selective','Non-Selective'});

% ---------------------------------------------------------------------------------------------------
% PLOTTING SECTION of POST TRAINING ONLY
subplot(1,2,2); hold on;

% PV | POST-T | RESP_NONRESP
b0 = bar(1,[nRespPostPV/nTotPostPV nNonRespPostPV/nTotPostPV],'grouped');
b0(1,1).FaceColor = [0 0.4470 0.7410];
b0(1,2).FaceColor = [0 0.4470 0.7410] .* 0.6;
stdVal2 = mannyPropSTD(nRespPostPV,nTotPostPV);
text(b0(1,1).XEndPoints,b0(1,1).YEndPoints+stdVal2,string(strcat("n=",num2str(nRespPostPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b0(1,2).XEndPoints,b0(1,2).YEndPoints+stdVal2,string(strcat("n=",num2str(nNonRespPostPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

% PV | POST-T | SEL_NONSEL
b1 = bar(2,[nSelPostPV/nRespPostPV nNonSelPostPV/nRespPostPV],'grouped');
b1(1,1).FaceColor = [0.4660 0.6740 0.1880];
b1(1,2).FaceColor = [0.4660 0.6740 0.1880] .* 0.6;
stdVal3 = mannyPropSTD(nSelPostPV,nRespPostPV);
text(b1(1,1).XEndPoints,b1(1,1).YEndPoints+stdVal3,string(strcat("n=",num2str(nSelPostPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b1(1,2).XEndPoints,b1(1,2).YEndPoints+stdVal3,string(strcat("n=",num2str(nNonSelPostPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');

% PV | POST-T | ERRORBARS (formated for legend accuracy)
errorbar(b1(1,1).XEndPoints,b1(1,1).YEndPoints,mannyPropSTD(nSelPrePV,nTotPrePV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b1(1,2).XEndPoints,b1(1,2).YEndPoints,mannyPropSTD(nNonSelPrePV,nTotPrePV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b0(1,1).XEndPoints,b0(1,1).YEndPoints,mannyPropSTD(nRespPrePV,nTotPrePV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b0(1,2).XEndPoints,b0(1,2).YEndPoints,mannyPropSTD(nNonRespPrePV,nTotPrePV),...
    Color=[0,0,0],LineWidth=1.4);

% NON-PV | POST-T | RESP_NONRESP
b2 = bar(3,[nRespPostNPV/nTotPostNPV nNonRespPostNPV/nTotPostNPV],'grouped');
b2(1,1).FaceColor = [0 0.4470 0.7410];
b2(1,2).FaceColor = [0 0.4470 0.7410] .* 0.6;
stdVal2 = mannyPropSTD(nRespPostNPV,nTotPostNPV);
text(b2(1,1).XEndPoints,b2(1,1).YEndPoints+stdVal2,string(strcat("n=",num2str(nRespPostNPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b2(1,2).XEndPoints,b2(1,2).YEndPoints+stdVal2,string(strcat("n=",num2str(nNonRespPostNPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
errorbar(b2(1,1).XEndPoints,b2(1,1).YEndPoints,mannyPropSTD(nRespPreNPV,nTotPreNPV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b2(1,2).XEndPoints,b2(1,2).YEndPoints,mannyPropSTD(nNonRespPreNPV,nTotPreNPV),...
    Color=[0,0,0],LineWidth=1.4);

% NON-PV | POST-T | SEL_NONSEL
b3 = bar(4,[nSelPostNPV/nRespPostNPV nNonSelPostNPV/nRespPostNPV],'grouped');
b3(1,1).FaceColor = [0.4660 0.6740 0.1880];
b3(1,2).FaceColor = [0.4660 0.6740 0.1880] .* 0.6;
stdVal3 = mannyPropSTD(nSelPostNPV,nRespPostNPV);
text(b3(1,1).XEndPoints,b3(1,1).YEndPoints+stdVal3,string(strcat("n=",num2str(nSelPostNPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b3(1,2).XEndPoints,b3(1,2).YEndPoints+stdVal3,string(strcat("n=",num2str(nNonSelPostNPV))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
errorbar(b3(1,1).XEndPoints,b3(1,1).YEndPoints,mannyPropSTD(nSelPreNPV,nTotPreNPV),...
    Color=[0,0,0],LineWidth=1.4);
errorbar(b3(1,2).XEndPoints,b3(1,2).YEndPoints,mannyPropSTD(nNonSelPreNPV,nTotPreNPV),...
    Color=[0,0,0],LineWidth=1.4);

% TITLE / LABELS / LEGEND
ax0 = gca;
ax0.XTick = [1.5 3.4];
ax0.XTickLabel = {'PV';'Non-PV'};
ax0.YLim = [0 1];
ylabel('Percentage of Neurons')
title('Categorization of Neurons Post-Training');
legend({'Responsive','Non-Responsive','Selective','Non-Selective'});



%% FIGURE 2: PREFFERED FREQUENCY FOR ALL NEURONS PROPORTIONAL PV/NON_PV


freqCountsPre = zeros(1,6); %bin array for counting combinations of selectivity and pref. frequency
for n = 1 : size(preD.bigOut,1)
    if preD.bigOut{n,12} > 0 & strcmp(preD.bigOut{n,5},'Overlap')
        switch preD.bigOut{n,4}
            case 100
                freqCountsPre(1) = freqCountsPre(1) + 1;
            case 300
                freqCountsPre(2) = freqCountsPre(2) + 1;
            case 500
                freqCountsPre(3) = freqCountsPre(3) + 1;
            case 700
                freqCountsPre(4) = freqCountsPre(4) + 1;
            case 900
                freqCountsPre(5) = freqCountsPre(5) + 1;
            case 1100
                freqCountsPre(6) = freqCountsPre(6) + 1;
        end
    end
end

freqCountsPost = zeros(1,6); %bin array for counting combinations of selectivity and Pref. frequency
for n = 1 : size(posD.bigOut,1)
    if posD.bigOut{n,12} > 0 & strcmp(posD.bigOut{n,5},'Overlap')
        switch posD.bigOut{n,4}
            case 100
                freqCountsPost(1) = freqCountsPost(1) + 1;
            case 300
                freqCountsPost(2) = freqCountsPost(2) + 1;
            case 500
                freqCountsPost(3) = freqCountsPost(3) + 1;
            case 700
                freqCountsPost(4) = freqCountsPost(4) + 1;
            case 900
                freqCountsPost(5) = freqCountsPost(5) + 1;
            case 1100
                freqCountsPost(6) = freqCountsPost(6) + 1;
        end
    end
end

freqCountsPre2 = zeros(1,6); %bin array for counting combinations of selectivity and pref. frequency
for n = 1 : size(preD.bigOut,1)
    if preD.bigOut{n,12} > 0 & strcmp(preD.bigOut{n,5},'GFP-only')
        switch preD.bigOut{n,4}
            case 100
                freqCountsPre2(1) = freqCountsPre2(1) + 1;
            case 300
                freqCountsPre2(2) = freqCountsPre2(2) + 1;
            case 500
                freqCountsPre2(3) = freqCountsPre2(3) + 1;
            case 700
                freqCountsPre2(4) = freqCountsPre2(4) + 1;
            case 900
                freqCountsPre2(5) = freqCountsPre2(5) + 1;
            case 1100
                freqCountsPre2(6) = freqCountsPre2(6) + 1;
        end
    end
end

freqCountsPost2 = zeros(1,6); %bin array for counting combinations of selectivity and Pref. frequency
for n = 1 : size(posD.bigOut,1)
    if posD.bigOut{n,12} > 0 & strcmp(posD.bigOut{n,5},'GFP-only')
        switch posD.bigOut{n,4}
            case 100
                freqCountsPost2(1) = freqCountsPost2(1) + 1;
            case 300
                freqCountsPost2(2) = freqCountsPost2(2) + 1;
            case 500
                freqCountsPost2(3) = freqCountsPost2(3) + 1;
            case 700
                freqCountsPost2(4) = freqCountsPost2(4) + 1;
            case 900
                freqCountsPost2(5) = freqCountsPost2(5) + 1;
            case 1100
                freqCountsPost2(6) = freqCountsPost2(6) + 1;
        end
    end
end


f = figure; hold on;
x1 = [100 300 500 700 900 1100];
y1 = [(freqCountsPre(1) + freqCountsPost(1)) / nTotPV, ...
    (freqCountsPre2(1) + freqCountsPost2(1)) / nTotNPV; 
(freqCountsPre(2) + freqCountsPost(2)) / nTotPV, ...
    (freqCountsPre2(2) + freqCountsPost2(2)) / nTotNPV; 
(freqCountsPre(3) + freqCountsPost(3)) / nTotPV, ...
    (freqCountsPre2(3) + freqCountsPost2(3)) / nTotNPV;
(freqCountsPre(4) + freqCountsPost(4)) / nTotPV, ...
    (freqCountsPre2(4) + freqCountsPost2(4)) / nTotNPV;
(freqCountsPre(5) + freqCountsPost(5)) / nTotPV, ...
    (freqCountsPre2(5) + freqCountsPost2(5)) / nTotNPV;
(freqCountsPre(6) + freqCountsPost(6)) / nTotPV, ...
    (freqCountsPre2(6) + freqCountsPost2(6)) / nTotNPV; ]
b = bar(x1,y1,'grouped');
pvGroup = freqCountsPre + freqCountsPost;
npvGroup = freqCountsPre2 + freqCountsPost2;
for nf = 1 : 6
text(b(1,1).XEndPoints(nf)+25,b(1,1).YEndPoints(nf)+0.02,string(strcat("n=",num2str(pvGroup(nf)))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom','Rotation',60);
text(b(1,2).XEndPoints(nf)+25,b(1,2).YEndPoints(nf)+0.02,string(strcat("n=",num2str(npvGroup(nf)))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom','Rotation',60);
end

title("Preferred frequency for all responsive neurons");
ylabel('% of Neurons');
xlabel('Frequency (Hz)');
ax = gca;
ax.XTick = x1;
ax.XTickLabel = {'100','300','500','700','900','1100'};
legend('PV','Non-PV');


%% FIGURE 2: PREFFERED FREQUENCY FOR ALL NEURONS PROPORTIONAL PV/NON_PV

nTotPre = length(isResponsivePreTrain);
freqCountsPre = zeros(1,6); %bin array for counting combinations of selectivity and pref. frequency
for n = 1 : size(preD.bigOut,1)
    if preD.bigOut{n,12} > 0
        switch preD.bigOut{n,4}
            case 100
                freqCountsPre(1) = freqCountsPre(1) + 1;
            case 300
                freqCountsPre(2) = freqCountsPre(2) + 1;
            case 500
                freqCountsPre(3) = freqCountsPre(3) + 1;
            case 700
                freqCountsPre(4) = freqCountsPre(4) + 1;
            case 900
                freqCountsPre(5) = freqCountsPre(5) + 1;
            case 1100
                freqCountsPre(6) = freqCountsPre(6) + 1;
        end
    end
end

nTotPost = length(isResponsivePostTrain);
freqCountsPost = zeros(1,6); %bin array for counting combinations of selectivity and Pref. frequency
for n = 1 : size(posD.bigOut,1)
    if posD.bigOut{n,12} > 0
        switch posD.bigOut{n,4}
            case 100
                freqCountsPost(1) = freqCountsPost(1) + 1;
            case 300
                freqCountsPost(2) = freqCountsPost(2) + 1;
            case 500
                freqCountsPost(3) = freqCountsPost(3) + 1;
            case 700
                freqCountsPost(4) = freqCountsPost(4) + 1;
            case 900
                freqCountsPost(5) = freqCountsPost(5) + 1;
            case 1100
                freqCountsPost(6) = freqCountsPost(6) + 1;
        end
    end
end


f = figure; hold on;
x1 = [100 300 500 700 900 1100];
y1 = [freqCountsPre(1) / nTotPre, freqCountsPost(1) / nTotPost; 
freqCountsPre(2) / nTotPre, freqCountsPost(2) / nTotPost; 
freqCountsPre(3) / nTotPre, freqCountsPost(3) / nTotPost; 
freqCountsPre(4) / nTotPre, freqCountsPost(4) / nTotPost; 
freqCountsPre(5) / nTotPre, freqCountsPost(5) / nTotPost; 
freqCountsPre(6) / nTotPre, freqCountsPost(6) / nTotPost; ]
b = bar(x1,y1,'grouped');
for nf = 1 : 6
text(b(1,1).XEndPoints(nf),b(1,1).YEndPoints(nf),string(strcat("n=",num2str(freqCountsPre(nf)))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b(1,2).XEndPoints(nf),b(1,2).YEndPoints(nf),string(strcat("n=",num2str(freqCountsPost(nf)))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
end

title("Preferred frequency for all responsive neurons");
ylabel('% of Neurons');
xlabel('Frequency (Hz)');
ax = gca;
ax.XTick = x1;
ax.XTickLabel = {'100','300','500','700','900','1100'};
legend('Pre-training','Post-training');



%% FIGURE 2: PREFFERED FREQUENCY FOR ALL NEURONS PROPORTIONAL 
% FOR EACH MOUSE

% Data Extraction m810
nTotPre = sum(is810pre);
freqCountsPre = zeros(1,6); 
for n = 1 : size(preD.bigOut,1)
    if preD.bigOut{n,12} > 0 & logical(is810pre(n))
        switch preD.bigOut{n,4}
            case 100
                freqCountsPre(1) = freqCountsPre(1) + 1;
            case 300
                freqCountsPre(2) = freqCountsPre(2) + 1;
            case 500
                freqCountsPre(3) = freqCountsPre(3) + 1;
            case 700
                freqCountsPre(4) = freqCountsPre(4) + 1;
            case 900
                freqCountsPre(5) = freqCountsPre(5) + 1;
            case 1100
                freqCountsPre(6) = freqCountsPre(6) + 1;
        end
    end
end
nTotPost = sum(is810post);
freqCountsPost = zeros(1,6);
for n = 1 : size(posD.bigOut,1)
    if posD.bigOut{n,12} > 0 & logical(is810post(n))
        switch posD.bigOut{n,4}
            case 100
                freqCountsPost(1) = freqCountsPost(1) + 1;
            case 300
                freqCountsPost(2) = freqCountsPost(2) + 1;
            case 500
                freqCountsPost(3) = freqCountsPost(3) + 1;
            case 700
                freqCountsPost(4) = freqCountsPost(4) + 1;
            case 900
                freqCountsPost(5) = freqCountsPost(5) + 1;
            case 1100
                freqCountsPost(6) = freqCountsPost(6) + 1;
        end
    end
end

% plotting section m810
f = figure; hold on;
x1 = [100 300 500 700 900 1100];
y1 = [freqCountsPre(1) / nTotPre, freqCountsPost(1) / nTotPost 
freqCountsPre(2) / nTotPre, freqCountsPost(2) / nTotPost; 
freqCountsPre(3) / nTotPre, freqCountsPost(3) / nTotPost; 
freqCountsPre(4) / nTotPre, freqCountsPost(4) / nTotPost; 
freqCountsPre(5) / nTotPre, freqCountsPost(5) / nTotPost; 
freqCountsPre(6) / nTotPre, freqCountsPost(6) / nTotPost; ];
b = bar(x1,y1,'grouped');
for nf = 1 : 6
text(b(1,1).XEndPoints(nf),b(1,1).YEndPoints(nf),string(strcat("n=",num2str(freqCountsPre(nf)))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b(1,2).XEndPoints(nf),b(1,2).YEndPoints(nf),string(strcat("n=",num2str(freqCountsPost(nf)))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
end
title("Preferred frequency for all responsive neurons (m810)");
ylabel('% of Neurons');
xlabel('Frequency Hz, (R) = rewarded');
ax = gca;
ax.XTick = x1;
ax.XTickLabel = {'100','300','500','700','900(R)','1100(R)'};
legend('Pre-training','Post-training');

%--------------------------------------------------------------------------
% Data Extraction m806
nTotPre = sum(is806pre);
freqCountsPre = zeros(1,6); 
for n = 1 : size(preD.bigOut,1)
    if preD.bigOut{n,12} > 0 & logical(is806pre(n))
        switch preD.bigOut{n,4}
            case 100
                freqCountsPre(1) = freqCountsPre(1) + 1;
            case 300
                freqCountsPre(2) = freqCountsPre(2) + 1;
            case 500
                freqCountsPre(3) = freqCountsPre(3) + 1;
            case 700
                freqCountsPre(4) = freqCountsPre(4) + 1;
            case 900
                freqCountsPre(5) = freqCountsPre(5) + 1;
            case 1100
                freqCountsPre(6) = freqCountsPre(6) + 1;
        end
    end
end
nTotPost = sum(is806post);
freqCountsPost = zeros(1,6); 
for n = 1 : size(posD.bigOut,1)
    if posD.bigOut{n,12} > 0 & logical(is806post(n))
        switch posD.bigOut{n,4}
            case 100
                freqCountsPost(1) = freqCountsPost(1) + 1;
            case 300
                freqCountsPost(2) = freqCountsPost(2) + 1;
            case 500
                freqCountsPost(3) = freqCountsPost(3) + 1;
            case 700
                freqCountsPost(4) = freqCountsPost(4) + 1;
            case 900
                freqCountsPost(5) = freqCountsPost(5) + 1;
            case 1100
                freqCountsPost(6) = freqCountsPost(6) + 1;
        end
    end
end

% plotting section m806
f = figure; hold on;
x1 = [100 300 500 700 900 1100];
y1 = [freqCountsPre(1) / nTotPre, freqCountsPost(1) / nTotPost 
freqCountsPre(2) / nTotPre, freqCountsPost(2) / nTotPost; 
freqCountsPre(3) / nTotPre, freqCountsPost(3) / nTotPost; 
freqCountsPre(4) / nTotPre, freqCountsPost(4) / nTotPost; 
freqCountsPre(5) / nTotPre, freqCountsPost(5) / nTotPost; 
freqCountsPre(6) / nTotPre, freqCountsPost(6) / nTotPost; ];
b = bar(x1,y1,'grouped');
for nf = 1 : 6
text(b(1,1).XEndPoints(nf),b(1,1).YEndPoints(nf),string(strcat("n=",num2str(freqCountsPre(nf)))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b(1,2).XEndPoints(nf),b(1,2).YEndPoints(nf),string(strcat("n=",num2str(freqCountsPost(nf)))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
end
title("Preferred frequency for all responsive neurons (m806)");
ylabel('% of Neurons');
xlabel('Frequency Hz, (R) = rewarded');
ax = gca;
ax.XTick = x1;
ax.XTickLabel = {'100(R)','300(R)','500','700','900','1100'};
legend('Pre-training','Post-training');


%--------------------------------------------------------------------------
% Data Extraction m792
nTotPre = sum(is792pre);
freqCountsPre = zeros(1,6); 
for n = 1 : size(preD.bigOut,1)
    if preD.bigOut{n,12} > 0 & logical(is792pre(n))
        switch preD.bigOut{n,4}
            case 100
                freqCountsPre(1) = freqCountsPre(1) + 1;
            case 300
                freqCountsPre(2) = freqCountsPre(2) + 1;
            case 500
                freqCountsPre(3) = freqCountsPre(3) + 1;
            case 700
                freqCountsPre(4) = freqCountsPre(4) + 1;
            case 900
                freqCountsPre(5) = freqCountsPre(5) + 1;
            case 1100
                freqCountsPre(6) = freqCountsPre(6) + 1;
        end
    end
end
nTotPost = sum(is792post);
freqCountsPost = zeros(1,6); 
for n = 1 : size(posD.bigOut,1)
    if posD.bigOut{n,12} > 0 & logical(is792post(n))
        switch posD.bigOut{n,4}
            case 100
                freqCountsPost(1) = freqCountsPost(1) + 1;
            case 300
                freqCountsPost(2) = freqCountsPost(2) + 1;
            case 500
                freqCountsPost(3) = freqCountsPost(3) + 1;
            case 700
                freqCountsPost(4) = freqCountsPost(4) + 1;
            case 900
                freqCountsPost(5) = freqCountsPost(5) + 1;
            case 1100
                freqCountsPost(6) = freqCountsPost(6) + 1;
        end
    end
end

% plotting section m792
f = figure; hold on;
x1 = [100 300 500 700 900 1100];
y1 = [freqCountsPre(1) / nTotPre, freqCountsPost(1) / nTotPost 
freqCountsPre(2) / nTotPre, freqCountsPost(2) / nTotPost; 
freqCountsPre(3) / nTotPre, freqCountsPost(3) / nTotPost; 
freqCountsPre(4) / nTotPre, freqCountsPost(4) / nTotPost; 
freqCountsPre(5) / nTotPre, freqCountsPost(5) / nTotPost; 
freqCountsPre(6) / nTotPre, freqCountsPost(6) / nTotPost; ];
b = bar(x1,y1,'grouped');
for nf = 1 : 6
text(b(1,1).XEndPoints(nf),b(1,1).YEndPoints(nf),string(strcat("n=",num2str(freqCountsPre(nf)))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
text(b(1,2).XEndPoints(nf),b(1,2).YEndPoints(nf),string(strcat("n=",num2str(freqCountsPost(nf)))),...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
end
title("Preferred frequency for all responsive neurons (m792)");
ylabel('% of Neurons');
xlabel('Frequency Hz, (R) = rewarded');
ax = gca;
ax.XTick = x1;
ax.XTickLabel = {'100','300','500','700','900(R)','1100(R)'};
legend('Pre-training','Post-training');








