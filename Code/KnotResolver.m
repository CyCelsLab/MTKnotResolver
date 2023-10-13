%% Created by Dhruv Khatri , IISER pune.
%KnotRsolver.m script implements contour segmentation and branch resolution
% pipeline for a give file. The inputs are read from a csv file containing
% file path and input parameters. Input parameters can be optimized using
% OptimizeParameters.m file. The parameters can then be saved in a
% csv file and read directly with KnotResolver.m. Please refer to the
% sample cvs provided on github for naming the variables and the format. 
% Hosted at: https://github.com/CyCelsLab/MTLoopResolve
% Parameters/names coded in the script are on line 51, 126
clear all; close all 
%% Load the input csv data 
[inputFile, pathFile]= uigetfile('.csv'); 
outputPath =fullfile(pathFile,'DemoOutput');
[status, msg, msgID] = mkdir(outputPath); 
input_data =readtable(fullfile(pathFile,inputFile));
%% Load parameters 
for h = 1:height(input_data)
% Load the csv parameters 
threshold = input_data.segThresh(h); % threshold value for imbinarize 
numIteration = input_data.segIteration(h); % iterations for active contour
contB = input_data.segContract(h); % contraction bias for activecontour function 
smoothF =input_data.segSmooth(h); % smooth factor for activecontour function
inputFile = input_data.Filename{h}; 
filenameOut = replace(inputFile,'.tif', '.mat') ; 
segfileName = replace(inputFile, '.tif', 'seg.tif'); 
outfilename = fullfile(outputPath,filenameOut); 
useLogP = input_data.useLogP{h}; % Whether to use LOG to improve intital guess. 
folderPath = input_data.FolderPath{h}; 
logPk =input_data.logPk(h); logPs=input_data.logPs(h);
filepath = fullfile(folderPath ,inputFile); 
nFrames = length(imfinfo(filepath)); 


data_struct = struct([]); 
selectedFilament = false; 
for f = 1:nFrames
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
    actCont = imclearborder(actCont); 
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

%     figure(1), imshow(imoverlay(medI, skelImage, 'red'),'InitialMagnification', 200,'Border','loose'); title(num2str(f));
%     title(num2str(f))
%     imwrite(getframe(gcf).cdata, fullfile(outputPath, segfileName),'Writemode', 'append');
%     patch([xx' nan], [yy' nan], [yy' nan], [yy' nan], 'edgecolor', 'interp'); 
%     pause(0.25)
end 
% Can be used for analysis of binary contours separately. 
save(outfilename,'data_struct'); 
end 

%% Branch Resolver 
% 2. Branch resolution
segOutPath = outputPath;
for f = 1:height(input_data)
    segFile = input_data.segFile{f}; 
    inputFile = fullfile(segOutPath, segFile); 
    ignoreFrames = input_data.IgnoreFrames(f); 
    editFrames = input_data.manualCheck(f); 
    listedit = string(editFrames);
    frameEdit = str2num(listedit);  
    listFra = string(ignoreFrames); 
    frameIgnore =str2num(listFra); 
    frameInvert = str2num(string(input_data.InvertManual(f))); 
    resFileName = replace(segFile, '.mat', ''); 
    resfullPath  = fullfile(segOutPath, resFileName); 
    load(inputFile) % Load the data_struct 
    resolveCoordinates = SingleFilamentResolver(data_struct, frameIgnore, ...
        frameEdit, frameInvert,resfullPath); 
end 
%% Analysis Plots 
resolvedStruct = 'resolveCoordinates_v4.mat'; 
% Process each struct and plot graphs 
for d = 1:height(input_data)
    filePath = fullfile(segOutPath, replace(input_data.segFile{d}, '.mat', ''), resolvedStruct); 
    load(filePath) % Output folder is resolve Coordinates 
    dataArray = AnalysisFile(resolveCoordinates,filePath); 
    [figurePath,~,~] = fileparts(filePath); 
    %print(gcf, '-dpdf', fullfile(figurePath, sprintf("Over%d.pdf",d)),'-r600')
    %figure(2), plot(cat(1,resolveCoordinates.FrameNumber).*10,bendingarray, 'lineWidth', 3.0); 
    print(gcf, '-dpdf', fullfile(figurePath, sprintf("Al%d.pdf",d)),'-r600')
    close all
    ProcessCSV(dataArray, 106/1000, 0, 10); 
    print(gcf, '-dpdf', fullfile(figurePath, sprintf("Tip angle%d",d)), '-r600')
end 