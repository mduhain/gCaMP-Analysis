% struct2Im2p.m
% a helpful script for manually assigning from 1x1 Struct --> im2p

function a = struct2Im2p(stringIn)

    load(stringIn);
    D=im2p;
    
    %RFP COMPONENTS    
    D.rfpPVImgPre = a.rfpPVImgPre;
    D.rfpExpTrace = a.rfpExpTrace;
    
    %GFP COMPONENTS
    D.gfp3FrameAvg = a.gfp3FrameAvg;
    D.gfpMotionCorrected = a.gfpMotionCorrected;
    D.templateImg4MotionCorrection = a.templateImg4MotionCorrection;
    D.gfpXCorrImg = a.gfpXCorrImg;
    D.somaticF = a.somaticF;
    D.somaticFDT = a.somaticFDT;
    
    %META DATA COMPONENTS
    D.expMetaData = a.expMetaData;
    D.numFrames = a.numFrames;
    D.stimOnTimes = a.stimOnTimes;
    D.stimOnFrameNum = a.stimOnFrameNum;
    D.stimFreqs = a.stimFreqs;
    D.stimAmps = a.stimAmps;
    D.stimDurs = a.stimDurs;
    D.audioTrials = a.audioTrials;
    D.imgDepth = a.imgDepth;
    D.graspTimes = a.graspTimes;
    D.graspDurs = a.graspDurs;
    
    %ROI IDENTIFICATION
    D.somaticROI_PixelLists = a.somaticROI_PixelLists;
    D.somaticROIBoundaries = a.somaticROIBoundaries;
    D.somaticROICenters = a.somaticROICenters;
    D.somaticRoiCounter = a.somaticRoiCounter;
    D.somaticROIs = a.somaticROIs;
    
    clear a;
    a = D;
    clear D;

end