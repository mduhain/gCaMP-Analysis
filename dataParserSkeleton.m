% dataParserSkeleton.m

startingDir = 'Y:\Michael\Haptics-2P\2P-Data';

cd(startingDir);
[mouseNum, ~] = analyzeDir();

for n = 1 : length(mouseNum) %mouse number loop
    mouseDir = strcat(startingDir,'\',mouseNum(n));
    cd(mouseDir);
    [dates, ~] = analyzeDir();
    for nd = 1 : length(dates) %date loop
        cd(strcat(mouseDir,'\',dates(nd)));
        [expNames,fileNames0] = analyzeDir();
        if sum(contains(fileNames0,'im2p')) > 0
            disp(strcat(mouseNum(n),'-',dates(nd),'-',expNames(ne)," im2p found in wrong place."));
        end
        for ne = 1 : length(expNames) %experiment loop
            cd(strcat(mouseDir,'\',dates(nd),'\',expNames(ne)));
            [dirNames, fileNames] = analyzeDir();
            %check for previous im2p files
            if sum(contains(fileNames,'im2p')) > 0
                disp(strcat(mouseNum(n),'-',dates(nd),'-',expNames(ne)," im2p found."));
                break;
            end  
            if isempty(dirNames)
                %new style Dec2023 directory, no depth folders
            else
                soloDataParser();
            end
        end
    end
end

disp('Data Parser complete.');
    
    



%% Local Functions
%--------------------------------------------------------------------------













