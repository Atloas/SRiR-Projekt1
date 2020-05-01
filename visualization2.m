clc;
clear;

dataString = fileread('resultdata.txt');
fullData = sscanf(dataString(10:end), '%d;%f;%f;%f\n', [4 Inf]);

trackedBodies = max(fullData(1, :)) + 1;
dataCount = size(fullData, 2)/trackedBodies;

myVideo = VideoWriter('visualization');
myVideo.FrameRate = 10;
open(myVideo)

for i = 1:trackedBodies:dataCount-trackedBodies
    scatter3(fullData(2, i:i+trackedBodies), fullData(3, i:i+trackedBodies), fullData(4, i:i+trackedBodies), 'filled');
    xlim([-5e8, 5e8]);
    ylim([-5e8, 5e8]);
    zlim([-5e8, 5e8]);
    view([0 0 1]);
    
    pause(0.01);
    
    frame = getframe(gcf);
    writeVideo(myVideo, frame);
end

close(myVideo)
