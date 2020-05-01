clear;
clc;

fid = fopen('resultdata.txt');
tline = fgetl(fid);

pos = [];
count = 0;
index = 1;
indexy = 1;
indexx = 1;

while ischar(tline)
    if mod(count,2)==1
        values = strsplit(tline,';');
        pos(indexy,indexx) = str2double(values{1});
        indexy = indexy + 1;
        pos(indexy,indexx) = str2double(values{2});
        indexy = indexy + 1;
        pos(indexy,indexx) = str2double(values{3});
        indexy = index;
        indexx = indexx + 1;
    end
    if count==23
        count = 0;
        indexx = 1;
        index = index + 3;
        indexy = index;
        while(ischar(tline) && (tline~="0:x;y;z"))
            tline = fgetl(fid);
        end
    end
    count = count + 1;
    tline = fgetl(fid);
end

X=[];
Y=[];
Z=[];
for i=1:3:size(pos,1)
    for j=1:1:size(pos,2)
        X(j)=pos(i,j);
        Y(j)=pos(i+1,j);
        Z(j)=pos(i+2,j);
    end
    pause(0.01);
    plt=scatter3(X,Y,Z);
    axis([-5e12 5e12 -5e12 5e12 -5e12 5e12]);
end

fclose(fid);