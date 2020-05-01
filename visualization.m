clc;
clear;

dataString = fileread('resultdata.txt');
fullData = sscanf(dataString(10:end), '%d;%f;%f;%f\n', [4 Inf]);

trackedBodies = max(fullData(1, :)) + 1;
datapointCount = size(fullData, 2)/trackedBodies;

fullData = reshape(fullData, 4, 36, datapointCount);

myVideo = VideoWriter('visualization');
myVideo.FrameRate = 24;
open(myVideo)
oservedBody = 4;
windowSize = 1e9;

for i = 1:datapointCount
    scatter3(fullData(2, end-1:end, i), fullData(3, end-1:end, i), fullData(4, end-1:end, i), 'filled', 'g');
    hold on;
    scatter3(fullData(2, 2:end-2, i), fullData(3, 2:end-2, i), fullData(4, 2:end-2, i), 'filled', 'b');
    scatter3(fullData(2, 1, i), fullData(3, 1, i), fullData(4, 1, i), 'filled', 'y');
    hold off;
    xlim([fullData(2, oservedBody, i) - windowSize, fullData(2, oservedBody, i) + windowSize]);
    ylim([fullData(3, oservedBody, i) - windowSize, fullData(3, oservedBody, i) + windowSize]);
    zlim([-5e11, 5e11]);
    view([0 0 1]);
    
    pause(0.01);
    
    frame = getframe(gcf);
    writeVideo(myVideo, frame);
end

close(myVideo)
