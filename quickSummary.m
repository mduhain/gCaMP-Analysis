%summary resp struct

function [f] = quickSummary(respStruct)

    if isstruct(respStruct)
    
        f = figure;
        subplot(2,5,[4,5,9,10]); hold on;
        title("All responses (1 sec post stim) per freq.");
        xlabel('frequency(Hz)');
        ylabel('dF/F');
        plot(100,respStruct.f100.avgRespPostStim,'k.'); hold on;
        plot(300,respStruct.f300.avgRespPostStim,'k.');
        plot(500,respStruct.f500.avgRespPostStim,'k.');
        plot(700,respStruct.f700.avgRespPostStim,'k.');
        plot(900,respStruct.f900.avgRespPostStim,'k.');
        plot(1100,respStruct.f1100.avgRespPostStim,'k.');
        hold off;
    
        % 100 HZ
        subplot(2,5,1); hold on;
        if respStruct.respFreqZP(1,2) < 0.05
            shadedErrorBar(-1:0.1:4.9,respStruct.f100.avgTrace,respStruct.f100.semTraces,'lineProps','g-');
        else
            shadedErrorBar(-1:0.1:4.9,respStruct.f100.avgTrace,respStruct.f100.semTraces,'lineProps','r-');
        end
        ax = gca;
        plot([0 0],[ax.YLim(1),ax.YLim(2)],'k-');
        plot([1 1],[ax.YLim(1),ax.YLim(2)],'k-');
        title("PSTH 100 HZ");
        xlabel('time (s) post-stim');
        ylabel('dF/F'); hold off;
    
        % 300 HZ
        subplot(2,5,2); hold on;
        if respStruct.respFreqZP(2,2) < 0.05
            shadedErrorBar(-1:0.1:4.9,respStruct.f300.avgTrace,respStruct.f300.semTraces,'lineProps','g-');
        else
            shadedErrorBar(-1:0.1:4.9,respStruct.f300.avgTrace,respStruct.f300.semTraces,'lineProps','r-');
        end
        ax = gca;
        plot([0 0],[ax.YLim(1),ax.YLim(2)],'k-');
        plot([1 1],[ax.YLim(1),ax.YLim(2)],'k-');
        title("PSTH 300 HZ");
        xlabel('time (s) post-stim');
        ylabel('dF/F'); hold off;
    
        % 500 HZ
        subplot(2,5,3); hold on;
        if respStruct.respFreqZP(3,2) < 0.05
            shadedErrorBar(-1:0.1:4.9,respStruct.f500.avgTrace,respStruct.f500.semTraces,'lineProps','g-');
        else
            shadedErrorBar(-1:0.1:4.9,respStruct.f500.avgTrace,respStruct.f500.semTraces,'lineProps','r-');
        end
        ax = gca;
        plot([0 0],[ax.YLim(1),ax.YLim(2)],'k-');
        plot([1 1],[ax.YLim(1),ax.YLim(2)],'k-');
        title("PSTH 500 HZ");
        xlabel('time (s) post-stim');
        ylabel('dF/F'); hold off;
    
        % 700 HZ
        subplot(2,5,6); hold on;
        if respStruct.respFreqZP(4,2) < 0.05
            shadedErrorBar(-1:0.1:4.9,respStruct.f700.avgTrace,respStruct.f700.semTraces,'lineProps','g-');
        else
            shadedErrorBar(-1:0.1:4.9,respStruct.f700.avgTrace,respStruct.f700.semTraces,'lineProps','r-');
        end
        ax = gca;
        plot([0 0],[ax.YLim(1),ax.YLim(2)],'k-');
        plot([1 1],[ax.YLim(1),ax.YLim(2)],'k-');
        title("PSTH 700 HZ");
        xlabel('time (s) post-stim');
        ylabel('dF/F'); hold off;
    
        % 900 HZ
        subplot(2,5,7); hold on;
        if respStruct.respFreqZP(5,2) < 0.05
            shadedErrorBar(-1:0.1:4.9,respStruct.f900.avgTrace,respStruct.f900.semTraces,'lineProps','g-');
        else
            shadedErrorBar(-1:0.1:4.9,respStruct.f900.avgTrace,respStruct.f900.semTraces,'lineProps','r-');
        end
        ax = gca;
        plot([0 0],[ax.YLim(1),ax.YLim(2)],'k-');
        plot([1 1],[ax.YLim(1),ax.YLim(2)],'k-');
        title("PSTH 900 HZ");
        xlabel('time (s) post-stim');
        ylabel('dF/F'); hold off;
    
        % 1100 HZ
        subplot(2,5,8); hold on;
        if respStruct.respFreqZP(6,2) < 0.05
            shadedErrorBar(-1:0.1:4.9,respStruct.f1100.avgTrace,respStruct.f1100.semTraces,'lineProps','g-');
        else
            shadedErrorBar(-1:0.1:4.9,respStruct.f1100.avgTrace,respStruct.f1100.semTraces,'lineProps','r-');
        end
        ax = gca;
        plot([0 0],[ax.YLim(1),ax.YLim(2)],'k-');
        plot([1 1],[ax.YLim(1),ax.YLim(2)],'k-');
        title("PSTH 1100 HZ");
        xlabel('time (s) post-stim');
        ylabel('dF/F'); hold off;

        f.Position = [200, 200, 1400, 600];

    elseif iscell(respStruct) %multiple entries, loop through them
        for n = 1 : size(respStruct,1)
            respStructLoc = respStruct{n};
            f = figure;
            subplot(2,5,[4,5,9,10]); hold on;
            title("All responses (1 sec post stim) per freq.");
            xlabel('frequency(Hz)');
            ylabel('dF/F');
            plot(100,respStructLoc.f100.avgRespPostStim,'k.'); hold on;
            plot(300,respStructLoc.f300.avgRespPostStim,'k.');
            plot(500,respStructLoc.f500.avgRespPostStim,'k.');
            plot(700,respStructLoc.f700.avgRespPostStim,'k.');
            plot(900,respStructLoc.f900.avgRespPostStim,'k.');
            plot(1100,respStructLoc.f1100.avgRespPostStim,'k.');
            hold off;
        
            % 100 HZ
            subplot(2,5,1); hold on;
            if respStructLoc.respFreqZP(1,2) < 0.05
                shadedErrorBar(-1:0.1:4.9,respStructLoc.f100.avgTrace,respStructLoc.f100.semTraces,'lineProps','g-');
            else
                shadedErrorBar(-1:0.1:4.9,respStructLoc.f100.avgTrace,respStructLoc.f100.semTraces,'lineProps','r-');
            end
            ax = gca;
            plot([0 0],[ax.YLim(1),ax.YLim(2)],'k-');
            plot([1 1],[ax.YLim(1),ax.YLim(2)],'k-');
            title("PSTH 100 HZ");
            xlabel('time (s) post-stim');
            ylabel('dF/F'); hold off;
        
            % 300 HZ
            subplot(2,5,2); hold on;
            if respStructLoc.respFreqZP(2,2) < 0.05
                shadedErrorBar(-1:0.1:4.9,respStructLoc.f300.avgTrace,respStructLoc.f300.semTraces,'lineProps','g-');
            else
                shadedErrorBar(-1:0.1:4.9,respStructLoc.f300.avgTrace,respStructLoc.f300.semTraces,'lineProps','r-');
            end
            ax = gca;
            plot([0 0],[ax.YLim(1),ax.YLim(2)],'k-');
            plot([1 1],[ax.YLim(1),ax.YLim(2)],'k-');
            title("PSTH 300 HZ");
            xlabel('time (s) post-stim');
            ylabel('dF/F'); hold off;
        
            % 500 HZ
            subplot(2,5,3); hold on;
            if respStructLoc.respFreqZP(3,2) < 0.05
                shadedErrorBar(-1:0.1:4.9,respStructLoc.f500.avgTrace,respStructLoc.f500.semTraces,'lineProps','g-');
            else
                shadedErrorBar(-1:0.1:4.9,respStructLoc.f500.avgTrace,respStructLoc.f500.semTraces,'lineProps','r-');
            end
            ax = gca;
            plot([0 0],[ax.YLim(1),ax.YLim(2)],'k-');
            plot([1 1],[ax.YLim(1),ax.YLim(2)],'k-');
            title("PSTH 500 HZ");
            xlabel('time (s) post-stim');
            ylabel('dF/F'); hold off;
        
            % 700 HZ
            subplot(2,5,6); hold on;
            if respStructLoc.respFreqZP(4,2) < 0.05
                shadedErrorBar(-1:0.1:4.9,respStructLoc.f700.avgTrace,respStructLoc.f700.semTraces,'lineProps','g-');
            else
                shadedErrorBar(-1:0.1:4.9,respStructLoc.f700.avgTrace,respStructLoc.f700.semTraces,'lineProps','r-');
            end
            ax = gca;
            plot([0 0],[ax.YLim(1),ax.YLim(2)],'k-');
            plot([1 1],[ax.YLim(1),ax.YLim(2)],'k-');
            title("PSTH 700 HZ");
            xlabel('time (s) post-stim');
            ylabel('dF/F'); hold off;
        
            % 900 HZ
            subplot(2,5,7); hold on;
            if respStructLoc.respFreqZP(5,2) < 0.05
                shadedErrorBar(-1:0.1:4.9,respStructLoc.f900.avgTrace,respStructLoc.f900.semTraces,'lineProps','g-');
            else
                shadedErrorBar(-1:0.1:4.9,respStructLoc.f900.avgTrace,respStructLoc.f900.semTraces,'lineProps','r-');
            end
            ax = gca;
            plot([0 0],[ax.YLim(1),ax.YLim(2)],'k-');
            plot([1 1],[ax.YLim(1),ax.YLim(2)],'k-');
            title("PSTH 900 HZ");
            xlabel('time (s) post-stim');
            ylabel('dF/F'); hold off;
        
            % 1100 HZ
            subplot(2,5,8); hold on;
            if respStructLoc.respFreqZP(6,2) < 0.05
                shadedErrorBar(-1:0.1:4.9,respStructLoc.f1100.avgTrace,respStructLoc.f1100.semTraces,'lineProps','g-');
            else
                shadedErrorBar(-1:0.1:4.9,respStructLoc.f1100.avgTrace,respStructLoc.f1100.semTraces,'lineProps','r-');
            end
            ax = gca;
            plot([0 0],[ax.YLim(1),ax.YLim(2)],'k-');
            plot([1 1],[ax.YLim(1),ax.YLim(2)],'k-');
            title("PSTH 1100 HZ");
            xlabel('time (s) post-stim');
            ylabel('dF/F'); hold off;

            f.Position = [200, 200, 1400, 600];
            xi = input('');
            close(f)

        end
    end

end