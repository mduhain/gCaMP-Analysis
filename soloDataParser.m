%soloDataParser.m
%
% mduhain 2023/12/14:
% Designed to be called from inside a 2P experiment folder
% i.e. /Haptics-2P/2P-Data/m000/2023-12-12/Exp000
% transverses individual depth folders (e.g. 'Z000/') and extract gfp
% images into an im2p placed one dir prior


%% Begin Dir

% id directory
startingDir = pwd;

%find im2p files and Z000 depth folders
[depths, fileNames] = analyzeDir();

%% Match variables from file path
bDir = strsplit(startingDir,'\');
for n=1:length(bDir)
    if ~isempty(regexp(bDir{n},'[m]\d{3}', 'once')) %'m000'
        mouseNum = bDir{n};
    end
    if ~isempty(regexp(bDir{n},'(\d{4})-(\d{2})-(\d{2})', 'once')) %'yyyy-mm-dd'
       ddate = bDir{n};
    end
    if ~isempty(regexp(bDir{n},'[exp]\d{3}', 'once')) %'exp000'
        expName = bDir{n};
    end
end

%% Check to make sure no im2p files exist already
if any(contains(fileNames,'im2p'))
    disp(strcat("IM2P found, skipping ",mouseDir,'\',dates(nd),'\',expNames(ne)));
    return
end     

%% Parse Data in Depth Folders
for k = 1 : length(depths) %loop through different depths
    %change to Z000 directory (depth)
    cd(strcat(startingDir,'\',depths(k)));

    %get info number of frames
    numFiles = length(dir('*.tif'));

    %extract GFP frames (downsample 256x256px @ 10fps)
    %extract whole frame avg trace from RFP frames
    if numFiles == 1
        % One Large TIF stack
        fileN = dir('*.tif');
        [gfpStack, rfpTrace, ~] = loadExp(fileN.name);
    elseif numFiles > 1
        % Many Individual TIFs
        gfpStack = gfpFromDir();
        rfpTrace = rfpTraceFromAllFiles();
    elseif numFiles == 0
        error("Check Dir, no files found...");
    end

    %save parameters into im2p
    a = im2p;

    % back out 1 directory to 'Exp000'
    cd(startingDir);

    % add other params. rfp imgs, vasculature 
    a.gfp3FrameAvg = gfpStack;
    a.rfpExpTrace = rfpTrace;

    % save updated im2p
    fileOutName = strcat(mouseNum,'-',ddate,'-',expName,'-',depths(k),'-im2p.mat');
    save(fileOutName,'a','-v7.3');

    %clear workspace
    clear a gfpStack rfpTrace

end

