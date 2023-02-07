%mkGif.m
% - makes a gif

figure;
h1 = gcf;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated1.gif';
for n = 1:1:1000
    % Draw plot for y = x.^n
    %x = 0:0.01:1;
    %y = x.^n;
    %plot(x,y) 
    %drawnow 
      % Capture the plot as an image 
      imagesc(rawFrames(:,:,1));
      view(n,35);
      frame = getframe(h1); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
  end