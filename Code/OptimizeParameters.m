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
threshold = 0.46; % Threshold value for binarization
numIteration = 10; % Number of iterations for active contour
contB = 0.4; % Contraction bias for activecontour function
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

outputPath =fullfile(inputPath,'DemoOutput');
[status, msg, msgID] = mkdir(outputPath); 
segfileName = replace(inputFile,".tiff", "segOverlay.tiff"); 
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
%% Save the initial output movie and check for segmentation errors. 
selectedFilament = false; 
for f = 1:length(nFrames)
    I = imread(filepath,f);
    I = imadjust(I);
    medI = medfilt2(I); 
     if strcmp(useLogP, "TRUE") 
        fil = -fspecial("log", logPk,logPs); 
        medF = imfilter(medI, fil); 
        BW1 = imbinarize(medF);
        %figure(2), imshow(BW1); 
    else
        % Threshold and contour optimization 
        BW1 = imbinarize(medI, threshold);
        %figure(2), imshow(BW1); 
     end 

    actCont = activecontour(medI, BW1, numIteration, 'Chan-vese', ...
        'ContractionBias',contB, 'SmoothFactor',smoothF);
    actCont = bwareafilt(actCont, 2); 
    connectedComp = bwconncomp(actCont,8);
    filaprops = regionprops(connectedComp, 'PixelIdxList'); 
    if ~selectedFilament
        overImage1 = imoverlay(medI, bwperim(actCont)); 
        figure(2), imshow(overImage1); title('Select one filamnet to track'); 
        [x, y] = ginput(1);
        indxSelect = sub2ind(size(I),ceil(y), ceil(x)); 
        for m = 1:length(filaprops)
            if ismember(indxSelect, filaprops(m).PixelIdxList)
                selectedFilament = filaprops(m).PixelIdxList;
            end 
        end
    end  
    nullImage = zeros(size(I)); 
    for l = 1:length(filaprops)
        if any(ismember(filaprops(l).PixelIdxList , selectedFilament))
            nullImage(filaprops(l).PixelIdxList)=1;
            selectedFilament = filaprops(l).PixelIdxList; 
        end 
    end 

    skelImage = bwskel(logical(nullImage));
    bp = find(bwmorph(skelImage, 'branchpoints')); 

    % Set the fixed tip (Feb/08/2023) 
    E1 = find(bwmorph(skelImage, 'endpoints')); 
    [ye, xe] = ind2sub(size(skelImage), E1); 
    
    dist_matrix = pdist2([ye, xe], [y, x]); 
    [~, indxm] = min(dist_matrix, [],1);
    referenceEnd = E1(indxm); 
    data_struct(f).FixedEnd = referenceEnd; 
    % Set the fixed tip 


    data_struct(f).Contour = find(skelImage); % save only skeleton indices
    data_struct(f).FrameNumber = f;
    data_struct(f).ImageSize = size(I); 
    if isempty(bp) 
        data_struct(f).FilamentType = 'Unbranched'; 
    else
        data_struct(f).FilamentType = 'Branched';
    end 
    figure(1), imshow(imoverlay(medI, skelImage, 'red'),'InitialMagnification', 200,'Border','loose'); title(num2str(f));
    title(num2str(f))
    
    imgData = getframe(gcf).cdata; 
    imgData = imresize(imgData,size(I)); 
    imwrite(imgData, fullfile(outputPath, segfileName),'Writemode', 'append');

    %patch([xx' nan], [yy' nan], [yy' nan], [yy' nan], 'edgecolor', 'interp'); 
    %pause(0.25)
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


