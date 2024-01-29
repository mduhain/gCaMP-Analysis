%mkGif.m
% - makes a gif


filename = 'gfpMotionCorrected.gif';

%draw first frame
imagesc(rawFrames(:,:,1));
colormap gray
axis square
ax = gca;
ax.XTickLabel = [];
ax.YTickLabel = [];

gif(filename,'DelayTime',1/30);

for n = 2 : 100
    imagesc(rawFrames(:,:,n));
    colormap gray
    axis square
    ax = gca;
    ax.XTickLabel = [];
    ax.YTickLabel = [];
    gif;
end
