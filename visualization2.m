clc;
clear;

dataString = fileread('resultdata.txt');
fullData = sscanf(dataString, '%d:%f;%f;%f;', [4 Inf]);

