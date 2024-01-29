% dataParser2
%
% mduhain 2023/12/14
% Improvements of original idea for dataParser. A function to be called at 
% end of day; searches through for new data on server to be converted --> im2p
%
%
startingDir = 'C:\Users\ramirezlab\Documents\Haptics-2P\2P-Data';
targetMice = ["m630","m627"];
cd(startingDir);
[mouseNum, ~] = analyzeDir();
disp(strcat(string(datetime("now"))," starting dataParser2.m"));
for n = 1 : length(mouseNum) %mouse number loop
    if ~any(contains(targetMice, mouseNum(n),'IgnoreCase',true))
        continue;
    end
    mouseDir = strcat(startingDir,'\',mouseNum(n));
    cd(mouseDir);
    [dates, ~] = analyzeDir();
    for nd = 1 : length(dates) %date loop
        %check for non-date folder, skip if true
        currentDate = char(dates(nd));
        if currentDate(1) ~= '2'
            continue;
        end
        cd(strcat(mouseDir,'\',dates(nd)));
        [expNames,fileNames0] = analyzeDir();
        if sum(contains(fileNames0,'im2p')) > 0
            disp(strcat(mouseNum(n),'-',dates(nd),'-',expNames(ne)," im2p found in wrong place."));
        end
        for ne = 1 : length(expNames) %experiment loop
            disp(strcat("Starting: ",mouseNum(n),'-',dates(nd),'-',expNames(ne),'...'));
            cd(strcat(mouseDir,'\',dates(nd),'\',expNames(ne)));
            [dirNames, fileNames] = analyzeDir();
            %check for previous im2p files
            if sum(contains(fileNames,'im2p')) > 0
                disp("im2p found proceeding to next folder...");
                break;
            end  
            if isempty(dirNames) %no depth data
                %new style Dec2023 directory, no depth folders
                for nf = 1 : length(fileNames)
                    if contains(fileNames(nf),'RFP','IgnoreCase',true) 
                        continue;
                    elseif contains(fileNames(nf),'vasc','IgnoreCase',true)
                        continue;
                    elseif contains(fileNames(nf),'tactilestim','IgnoreCase',true)
                        %local data parser
                        disp(string(datetime("now")));
                        [gfpStack, rfpTrace, ~] = loadExp(fileNames(nf));
                        a = im2p;
                        a.gfp3FrameAvg = gfpStack;
                        a.rfpExpTrace = rfpTrace;
                        disp('Exporting im2p...');
                        fNameHere = char(fileNames(nf));
                        if any(fNameHere=='Z')
                            indx = find(fNameHere=='Z');
                            depthChar = fNameHere(indx:indx+3);
                        elseif any(fNameHere=='z')
                            indx = find(fNameHere=='z');
                            depthChar = fNameHere(indx:indx+3);
                        else 
                            depthChar = 'Z000';
                        end
                        depthChar = fNameHere(indx:indx+3);
                        fileOutName = strcat(mouseNum(n),'-',dates(nd),'-',expNames(ne),'-',depthChar,'-im2p.mat');
                        save(fileOutName,'a','-v7.3');
                    end
                end
            else
                soloDataParser();
            end
        end
    end
end

disp('Data Parser complete.');
    
    



%% Local Functions
%--------------------------------------------------------------------------













