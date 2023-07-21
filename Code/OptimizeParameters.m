%% Created By Dhruv Khatri , IISER Pune. 
%{
The script can be used for optimizing semgentation parameters for a given
time series. The parameters can then be appended or exported to a csv file
that can be read by the KnotResolver function. 
%} 

%% Load the input times series 
[inputFile, pathFile]= uigetfile("*.csv*"); 
fullPath = fullfile(pathFile,inputFile); 
imInfo = imfinfo(fullPath);
% Display first frame 
I = imread(fullPath); 
figure(1) , imshow(I, []);
%% Segmentation optimizer

