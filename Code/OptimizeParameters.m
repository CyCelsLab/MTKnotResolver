%% Created By Dhruv Khatri, IISER Pune.
%{
The script optimizes segmentation parameters for a given time series. The
parameters can be appended or exported to a CSV file for use with the
KnotResolver function.
%}

%% Input csv file selection
[inputFile, inputPath] = uigetfile('*.csv', 'Select MT time series');
readCsv = readtable(fullfile(inputPath, inputFile));
%% Process the files 
for ff = 1:height(readCsv)
    filepath = fullfile(readCsv.FolderPath{ff}, readCsv.Filename{ff});
    nFrames = imfinfo(filepath);

    visFrame = 1; % Select the frame for tuning segmentation parameters
    filenameOut = strrep(inputFile, ".csv", ".mat");
    threshold = readCsv.segThresh(ff);
    numIteration = readCsv.segIteration(ff);
    contB = readCsv.segContract(ff);
    smoothF = readCsv.segSmooth(ff);
    resolveFile = '';
    requiresManual = false;
    manualcheck = [];
    IgnoreFrames = [];
    useLogP = false;
    logPk = 7;
    logPs = 0.5;
    invertSeries = false;
    InvertManual = false;
    tuneParam = true;

    outputPath = fullfile(readCsv.FolderPath{ff}, 'DemoOutput');
    [status, msg, msgID] = mkdir(outputPath);
    segfileName = strrep(readCsv.Filename{ff}, ".tif", "segOverlay.tiff");

    while tuneParam
        I = imread(filepath, visFrame);
        I = imadjust(I);
        medI = medfilt2(I);

        figure('Name', 'Segmentation Tuning', 'Position', [100, 100, 1600, 800]);
        originalAxes = axes('Position', [0.05, 0.2, 0.4, 0.7]);
        imshow(medI);
        title('Original Image');

        segmentedAxes = axes('Position', [0.55, 0.2, 0.4, 0.7]);
        title('Segmented Image');

       
        segmented = imbinarize(medI, threshold);
        actCont = activecontour(medI, segmented, numIteration, 'Chan-Vese', 'ContractionBias', contB, 'SmoothFactor', smoothF);
        overImage = imoverlay(medI, bwskel(actCont), 'red'); 
        overImage = imoverlay(overImage , bwperim(actCont), 'yellow'); 
        imshow(overImage, 'Parent', segmentedAxes);
        title(segmentedAxes, ['Segmented Image (Threshold = ' num2str(threshold) ')']);
        scrollbar = uicontrol('Style', 'slider', 'Min', 0, 'Max', 1, 'Value', threshold, ...
             'Position', [100, 10, 600, 20], 'Callback', @(src,~) updateSegmentation(src, medI, numIteration, contB, smoothF,segmentedAxes));
        
        frameScrollbar = uicontrol('Style', 'slider', 'String', 'Frame','Min', 1, 'Max', length(nFrames), 'Value', visFrame, ...
        'Position', [100, 40, 600, 20], 'Callback', @(src, ~) updateFrame(src, filepath,scrollbar.Value, numIteration, contB, smoothF,segmentedAxes));

       continueButton = uicontrol('Style', 'pushbutton', 'String', 'Continue', ...
        'Position', [100, 60, 100, 30], 'Callback', @continueScript);
        
        % Pause the execution of the script until the "Continue" button is clicked
        uiwait;
        threshold = scrollbar.Value; 
        readCsv.segThresh(ff) = threshold;
        tuneParam = false; 

    end
    close all 
    selectedFilament = false;

    fileToDelete = fullfile(outputPath, segfileName);

    if exist(fileToDelete, 'file') == 2
        delete(fileToDelete);
        disp(['Deleted ' fileToDelete]);
    else
        disp(['File ' fileToDelete ' does not exist.']);
    end

    data_struct = struct([]);
    for f = 1:length(nFrames)
        I = imread(filepath, f);
        I = imadjust(I);
        medI = medfilt2(I);

        if useLogP
            fil = -fspecial('log', logPk, logPs);
            medF = imfilter(medI, fil);
            BW1 = imbinarize(medF);
        else
            BW1 = imbinarize(medI, threshold);
        end

        actCont = activecontour(medI, BW1, numIteration, 'Chan-Vese', 'ContractionBias', contB, 'SmoothFactor', smoothF);
        actCont = imclearborder(actCont); 
        actCont = bwareafilt(actCont, 3);
        connectedComp = bwconncomp(actCont, 8);
        filaprops = regionprops(connectedComp, 'PixelIdxList');

        if ~selectedFilament
            overImage1 = imoverlay(medI, bwperim(actCont));
            figure(2), imshow(overImage1, 'InitialMagnification',600);
            title('Select one filament to track');
            [x, y] = ginput(1);
            indxSelect = sub2ind(size(I), ceil(y), ceil(x));
            for m = 1:length(filaprops)
                if ismember(indxSelect, filaprops(m).PixelIdxList)
                    selectedFilament = filaprops(m).PixelIdxList;
                end
            end
        end

        overlap = 0; 
        distance = Inf; 
        for l = 1:length(filaprops)
            if any(ismember(filaprops(l).PixelIdxList, selectedFilament))% ismember(indxSelect, filaprops(m).PixelIdxList)
                val = sum(ismember(filaprops(l).PixelIdxList, selectedFilament)); 
%                 indicCur = filaprops(l).PixelIdxList; 
%                 preVCur = selectedFilament; 
%                 distancMet = calculateDistanceBetweenIndices(indicCur,   preVCur, I); 
                if overlap <  val
                    overlap = val ; 
                    distance  = distancMet; 
                    nullImage = zeros(size(I)); 
                    nullImage(filaprops(l).PixelIdxList) = 1;
                    %selectedFilament = filaprops(l).PixelIdxList;
                end 
            end
        end

        skelImage = bwskel(logical(nullImage));
        bp = find(bwmorph(skelImage, 'branchpoints'));

        E1 = find(bwmorph(skelImage, 'endpoints'));
        [ye, xe] = ind2sub(size(skelImage), E1);

        dist_matrix = pdist2([ye, xe], [y, x]);
        [~, indxm] = min(dist_matrix, [], 1);
        referenceEnd = E1(indxm);
        data_struct(f).FixedEnd = referenceEnd;

        data_struct(f).Contour = find(skelImage);
        data_struct(f).FrameNumber = f;
        data_struct(f).ImageSize = size(I);
        if isempty(bp)
            data_struct(f).FilamentType = 'Unbranched';
        else
            data_struct(f).FilamentType = 'Branched';
        end

        figure(1), imshow(imoverlay(medI, skelImage, 'red'), 'InitialMagnification', 200, 'Border', 'loose');
        title(num2str(f));

        imgData = getframe(gcf).cdata;
        imgData = imresize(imgData, size(I));
        imwrite(imgData, fullfile(outputPath, segfileName), 'WriteMode', 'append');
    end

    headers = {'Filename', 'FolderPath', 'segFile', 'segThresh', 'segIteration', ...
        'segContract', 'segSmooth', 'resFile', 'requiresManual', 'manualCheck', ...
        'IgnoreFrames', 'useLogP', 'logPk', 'logPs', 'invertSeries', 'InvertManual'};

    outputString = sprintf('filenameOut: %s\ninputPath: %s\nthreshold: %.2f\nnumIteration: %d\ncontB: %.2f\nsmoothF: %.2f\nresolveFile: %s\nrequiresManual: %d\nmanualcheck: %s\nIgnoreFrames: %s\nuseLogP: %d\nlogPk: %d\nlogPs: %.2f\ninvertSeries: %d\nInvertManual: %d\n', ...
        filenameOut, inputPath, threshold, numIteration, contB, smoothF, resolveFile, requiresManual, manualcheck, IgnoreFrames, useLogP, logPk, logPs, invertSeries, InvertManual);

    disp(outputString);
end

function updateSegmentation(scrollbar, medI, numIteration, contB, smoothF,segmentedAxes)
    threshold = get(scrollbar, 'Value');
    segmented = imbinarize(medI, threshold);
    actCont = activecontour(medI, segmented, numIteration, 'Chan-Vese', 'ContractionBias', contB, 'SmoothFactor', smoothF);
    overImage = imoverlay(medI, bwskel(actCont), 'red'); 
    overImage = imoverlay(overImage , bwperim(actCont), 'yellow'); 

    imshow(overImage, 'Parent', segmentedAxes);
    title(segmentedAxes, ['Segmented Image (Threshold = ' num2str(threshold) ')']);
end

function updateFrame(scrollbar, filepath,threshold, numIteration, contB, smoothF,segmentedAxes)
    frame = get(scrollbar, 'Value'); 
    I = imread(filepath, ceil(frame));
    I = imadjust(I);
    medI = medfilt2(I);
    segmented = imbinarize(medI, threshold);
    actCont = activecontour(medI, segmented, numIteration, 'Chan-Vese', 'ContractionBias', contB, 'SmoothFactor', smoothF);
    overImage = imoverlay(medI, bwskel(actCont), 'red'); 
    overImage = imoverlay(overImage , bwperim(actCont), 'yellow'); 

    imshow(overImage, 'Parent', segmentedAxes);
    title(segmentedAxes, ['Segmented Image (Threshold = ' num2str(threshold) ')']);
end


function continueScript(~, ~)
    % Close the current figure to resume script execution
    uiresume;
end

function totalDistance = calculateDistanceBetweenIndices(indicCu, preVCur, I)
    % Convert indices to pixel coordinates
    [yCu, xCu] = ind2sub(size(I), indicCu);
    [yPreV, xPreV] = ind2sub(size(I), preVCur);
    
    % Calculate the Euclidean distance between the two points
    distance = sqrt((xCu - xPreV).^2 + (yCu - yPreV).^2);
    
    totalDistance = sum(distance);
end


