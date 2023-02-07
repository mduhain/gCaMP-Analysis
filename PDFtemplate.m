%% Template for exporting PDFs
% J. Fritzinger, updated 8/30/22
%
% M. DuHain, use for 2P data analysis early October

%% Initialize report
import mlreportgen.dom.*
import mlreportgen.report.*

% Title information 
report_name = 'm748-2022-10-07-Exp002-tuningSEM';
headingText = 'Tactile Stim: m748 on 2022-10-07, Exp002';

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

% Example plot 


for n=1:size(a.somaticF,1)
    %PSTH by Freq.
    fig1 = figure;
    title = strcat('Plot',num2str(4*n-3));
    [tuningMat,tuningMatSTD,yRange] = a.psth2(a.somaticFDT(n,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,n,4);
    [plt1, images] = addToPDF(images, fig1, title, [7.5, 4]);
    append(rpt, plt1);
    
    %TUNING CURVE
    fig2 = figure;
    title = strcat('Plot',num2str(4*n-2));
    a.simpleTuningCurve(a.stimFreqs,tuningMat,tuningMatSTD,n);
    [plt2, images] = addToPDF(images, fig2, title, [4.2, 3]);
    %append(rpt, plt2);

    %AUDIO-STIM RESPONSE
    fig4 = figure;
    title = strcat('Plot',num2str(4*n-1));
    a.psthAudioTrials(a.somaticFDT(n,:),a.stimOnFrameNum,a.stimFreqs,a.audioTrials,n,yRange);
    [plt4, images] = addToPDF(images, fig4, title, [3.5, 2.5]);

    %JOIN TUNING CURVE AND AUDIO-STIM RESPONSE
    table = Table({plt2,plt4});
    table.Style = {Width('100%'), ResizeToFitContents(false)};
    table.BorderColor = 'White';
    append(rpt, table);

    %ROI LOCATION
    fig3 = figure;
    title = strcat('Plot',num2str(4*n));
    imagesc(a.gfpXCorrImg); hold on;
    plot(a.somaticROIBoundaries{1,n}{1,1}(:,2),a.somaticROIBoundaries{1,n}{1,1}(:,1),'LineWidth',2,'Color',[1 0 0]);
    [plt3, images] = addToPDF(images, fig3, title, [3.5, 2.5]);
    append(rpt, plt3);
end



% table = Table({plt2,plt3});
% table.Style = {Width('100%'), ResizeToFitContents(false)};
% table.BorderColor = 'White';
% append(rpt, table);

% Closes and opens PDF to view
close(rpt);
for i = 1:length(images)
    delete(images{1,i}.Path);
end
rptview(rpt)

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
