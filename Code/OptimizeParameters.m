%% Created By Dhruv Khatri, IISER Pune.
%{
The script optimizes segmentation parameters for a given time series. The
parameters can be appended or exported to a CSV file for use with the
KnotResolver function.
%}

%% Input image selection
[inputFile, inputPath] = uigetfile('*.tif*', 'Select MT time series');
filepath = fullfile(inputPath, inputFile);
nFrames = imfinfo(filepath);
%% Set parameters for segmentation
% Define and initialize various segmentation parameters
visFrame = 1; % Select the frame for tuning segmentation parameters
filenameOut = replace(inputFile, ".tiff", ".mat");
threshold = 0.30; % Threshold value for binarization
numIteration = 50; % Number of iterations for active contour
contB = 0.25; % Contraction bias for activecontour function
smoothF = 0.5; % Smooth factor for activecontour function
resolveFile = "";
requiresManual = false; % Flag indicating manual intervention
manualcheck = []; % manual checking on these frames data
IgnoreFrames = []; % Frames to be ignored
useLogP = false; % Use Laplacian of Gaussian filter
logPk = 7; % filter parameters fize
logPs = 0.5; % filter parameter sigma
invertSeries = false; % Invert series flag
InvertManual = false; % Manual inversion flag
tuneParam = true; % Flag to indicate parameter tuning
%% Adjusting parameters and image processing
if tuneParam
    % Read and preprocess the image
    I = imread(filepath, visFrame);
    I = imadjust(I);
    medI = medfilt2(I);
    
    if useLogP
        % Apply Laplacian of Gaussian filter and binarize
        fil = -fspecial("log", logPk, logPs);
        medF = imfilter(I, fil);
        BW1 = imbinarize(medF);
    else
        % Threshold and binarize the image
        BW1 = imbinarize(medI, threshold);
    end
    
    % Perform active contour segmentation
    actCont = activecontour(medI, BW1, numIteration, 'Chan-vese', 'ContractionBias', contB, 'SmoothFactor', smoothF);
    
    % Visualize the result
    overImage1 = imoverlay(medI, bwperim(actCont));
    skelImage = bwskel(actCont);
    overImage2 = imoverlay(overImage1, skelImage, 'red');
    figure(2), imshow(overImage2);
end

%% Output Final segmentation parameters
% Define headers for the table
headers = {'Filename', 'FolderPath', 'segFile', 'segThresh', 'segIteration', ...
    'segContract', 'segSmooth', 'resFile', 'requiresManual', 'manualCheck', ...
    'IgnoreFrames', 'useLogP', 'logPk', 'logPs', 'invertSeries', 'InvertManual'};

% Convert some variables to cell arrays for table construction
filenameOut = cellstr(filenameOut);
inputPath = cellstr(inputPath);
filenameOut = cellstr(filenameOut);
resolveFile = cellstr(resolveFile);

% % Create a table with headers
% parameterTable = table(filenameOut, inputPath, filenameOut, threshold, ...
%     numIteration, contB, smoothF, ...
%     resolveFile, requiresManual, manualcheck, IgnoreFrames, useLogP, logPk, logPs, ...
%     invertSeries, InvertManual, 'VariableNames', headers);

outputString = sprintf('filenameOut: %s\ninputPath: %s\nthreshold: %.2f\nnumIteration: %d\ncontB: %.2f\nsmoothF: %.2f\nresolveFile: %s\nrequiresManual: %d\nmanualcheck: %s\nIgnoreFrames: %s\nuseLogP: %d\nlogPk: %d\nlogPs: %.2f\ninvertSeries: %d\nInvertManual: %d\n', ...
    filenameOut{1}, inputPath{1}, threshold, numIteration, contB, smoothF, resolveFile{1}, requiresManual, manualcheck, IgnoreFrames, useLogP, logPk, logPs, invertSeries, InvertManual);

% Display the formatted string
disp(outputString);


