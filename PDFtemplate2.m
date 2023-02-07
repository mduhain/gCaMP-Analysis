%% Template for exporting PDFs
% J. Fritzinger, updated 8/30/22
%
% M. DuHain, use for 2P data analysis early November

%% Initialize report
import mlreportgen.dom.*
import mlreportgen.report.*

% Experiment information 
mouse = "810";
date = "2023-01-24";
expNum = "001";

%create titles (.pdf and page title)
report_name = strcat("m",mouse,"-",date,"-Exp",expNum,"-report");
headingText = strcat("Tactile Stim: m",mouse," on ",date," Exp",expNum);

% Initialize report
images = {}; %hold all plots as images, need to delete when finished
rpt = Document(report_name, 'pdf'); %initialize report document as pdf
open(rpt);
pm = rpt.CurrentPageLayout;

% Set page header dimensions
pm.PageMargins.Top = '0.1in';
pm.PageMargins.Header = '0.1in';
pm.PageMargins.Bottom = '0.1in';
pm.PageMargins.Footer = '0.1in';
pm.PageMargins.Left = '0.2in';
pm.PageMargins.Right = '0.2in';

% Generate page heading
h = Heading(2, headingText);
b = Border();
b.BottomStyle = 'single';
b.BottomColor = 'LightGray';
b.BottomWidth = '1pt';
h.Style = [h.Style {Color('Black'), b}, {PageBreakBefore()}];
append(rpt,h);

% figB = figure;
% title = 'TextInfo';
% figure;
% ax = gca
% ax.XTickLabel = [];
% ax.YTickLabel = [];
% text(0.02,0.95,strcat(num2str(a.redSomaticRoiCounter+a.somaticRoiCounter)," total ROIs."));
% text(0.02,0.85,strcat(num2str(a.redSomaticRoiCounter),"(rfp) and ",num2str(a.somaticRoiCounter),"(gfp)"));


for n=1:size(a.somaticF,1)
    %Responsivity
    fig1 = figure;
    title = strcat('Plot',num2str(4*n-3));
    %[tuningMat,tuningMatSTD,yRange] = a.psth2(a.somaticFDT(n,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,n,4);
    out = a.getResponsivity(a.somaticFDT(n,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,1);
    [plt1, images] = addToPDF(images, fig1, title, [8.5, 3]);
    append(rpt, plt1);

    %PSTH by Freq.
    fig2 = figure;
    title = strcat('Plot',num2str(4*n-2));
    %[tuningMat,tuningMatSTD,yRange] = a.psth2(a.somaticFDT(n,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,n,4);
    out = a.getTraces(a.somaticFDT(n,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,1);
    [plt2, images] = addToPDF(images, fig2, title, [8.5, 4]);
    append(rpt, plt2);

    %AUDIO-STIM RESPONSE
    fig4 = figure;
    title = strcat('Plot',num2str(4*n-1));
    a.psthAudioTrials(a.somaticFDT(n,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,n,out.yAxMinMax);
    [plt4, images] = addToPDF(images, fig4, title, [3.5, 2.5]);

    %ROI LOCATION
    fig3 = figure;
    title = strcat('Plot',num2str(4*n));
    imagesc(a.gfpXCorrImg); hold on;
    plot(a.somaticROIBoundaries{1,n}{1,1}(:,2),a.somaticROIBoundaries{1,n}{1,1}(:,1),'LineWidth',2,'Color',[1 0 0]);
    [plt3, images] = addToPDF(images, fig3, title, [3.5, 2.5]);

    %JOIN TUNING CURVE AND AUDIO-STIM RESPONSE
    table = Table({plt3,plt4});
    table.Style = {Width('100%'), ResizeToFitContents(false)};
    table.BorderColor = 'White';
    append(rpt, table);
end

np = size(a.somaticF,1);

% % NOW FOR RED
% for n=1:size(a.redSomaticF,1)
%     %Responsivity
%     fig1 = figure;
%     title = strcat('Plot',num2str(4*n-3+np));
%     %[tuningMat,tuningMatSTD,yRange] = a.psth2(a.somaticFDT(n,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,n,4);
%     out = a.getResponsivity(a.redSomaticFDT(n,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,1);
%     [plt1, images] = addToPDF(images, fig1, title, [8.5, 3]);
%     append(rpt, plt1);
% 
%     %PSTH by Freq.
%     fig2 = figure;
%     title = strcat('Plot',num2str(4*n-2+np));
%     %[tuningMat,tuningMatSTD,yRange] = a.psth2(a.somaticFDT(n,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,n,4);
%     out = a.getTraces(a.redSomaticFDT(n,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,1);
%     [plt2, images] = addToPDF(images, fig2, title, [8.5, 4]);
%     append(rpt, plt2);
% 
%     %AUDIO-STIM RESPONSE
%     fig4 = figure;
%     title = strcat('Plot',num2str(4*n-1+np));
%     a.psthAudioTrials(a.redSomaticFDT(n,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,n,out.yAxMinMax);
%     [plt4, images] = addToPDF(images, fig4, title, [3.5, 2.5]);
% 
%     %ROI LOCATION
%     fig3 = figure;
%     title = strcat('Plot',num2str(4*n+np));
%     imagesc(a.gfpXCorrImg); hold on;
%     plot(a.redSomaticROIBoundaries{1,n}{1,1}(:,2),a.redSomaticROIBoundaries{1,n}{1,1}(:,1),'LineWidth',2,'Color',[1 0 0]);
%     [plt3, images] = addToPDF(images, fig3, title, [3.5, 2.5]);
% 
%     %JOIN TUNING CURVE AND AUDIO-STIM RESPONSE
%     table = Table({plt3,plt4});
%     table.Style = {Width('100%'), ResizeToFitContents(false)};
%     table.BorderColor = 'White';
%     append(rpt, table);
% end





close(rpt);
for i = 1:length(images)
    delete(images{1,i}.Path);
end
% 
% rptview(rpt)

%% 
function [img, images] = addToPDF(images, fig, title, values)
import mlreportgen.dom.*

% Set figure size, recommended
%values = [4.5, 3.5];
fig.PaperSize = values;
fig.PaperPosition = [0 0 values];
fig.Units = 'inches';
fig.Position(3:4) = values;

% Add the plot to the document
name = sprintf('%s.svg', title);
print(fig, name, '-dsvg');
img = Image(name);
delete(fig) %delete plot figure window
images = [images {img}];

end
