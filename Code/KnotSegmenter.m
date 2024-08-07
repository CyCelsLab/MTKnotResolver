%% Created By Dhruv Khatri, IISER Pune.
%{
The script optimizes segmentation parameters for a given time series. The
parameters can be appended or exported to a CSV file for use with the
KnotResolver function.
%}
close all; clear; 
%% 1. Input csv file selection
[inputFile, inputPath] = uigetfile('*.csv', 'Select MT time series');
expectedHeaders = {'FileName', 'FolderPath'}; 

readCsv = readtable(fullfile(inputPath, inputFile), ...
    'ReadVariableNames',true, 'Delimiter',',');

actualHeaders = readCsv.Properties.VariableNames;

% Compare the actual headers with the expected headers
if isequal(actualHeaders, expectedHeaders)
    disp('Correct headers found.....');
else
    disp('Script stopped');
    return 
end 
 
% Specify the filename to check
filename = replace(inputFile, '.csv', '_optimized.csv');  % Replace with your actual CSV filename
optFilename = fullfile(inputPath, filename); 
% Check if the file exists in the current directory

if exist(optFilename, 'file') == 2
    % Read the first line of the CSV file to check headers
    disp(['Optimized file' , filename ,' exits in the directory'])
    optCsv = readtable(optFilename,  'ReadVariableNames',true, 'Delimiter',',');
    headerLine = optCsv.Properties.VariableNames;

    % Define the expected headers (replace with your actual headers)
    expectedHeaders = {'FileName', 'FolderPath', 'segThresh', 'segIteration', 'segContract', ...
    	'segSmooth', 'requiresManual','IgnoreFrames',...
    	'useLogP','logPk','logPs','InvertManual', 'manualCheck'} ;

    % Compare the actual headers with the expected headers
    actualHeaders = headerLine;
    if isequal(actualHeaders, expectedHeaders)
        userInput = input('OPT CSV file found. Do you want to re write the data (y/n): ', 's'); 
        if strcmpi(userInput, 'y')
            disp('File will be re written');
            optSeg = readCsv;
        elseif strcmpi(userInput, 'n')
            disp('Delete or rename the file')
            return 
        else
            return 
        end 
    else
        disp('OPT CSV seems incorrect, file will be over-written deleting the file') 
        optSeg = readCsv; 
    end
else
    disp(['Starting from the ' ...
        'first entry creating ', filename]);
    optSeg = readCsv;
end
%% 2. Default segmentation parameter structure
segD.segThresh =  0.46; % Threshold for imbinzrize 
segD.segIteration = 10; % Number of iterations for active contours
segD.segContract =  0.4; % Contraction bias for active contours
segD.segSmooth = 0.2; % Smooth factor for active contours 

% These parameters can be edited in the {inputFileName}_optimized.csv after checking the
% segmentation output 
segD.requiresManual = false; % Does any frame requires manual correction 
segD.manualCheck = []; % Frames that require manual correction 
segD.IgnoreFrames =  []; % Frame where segmentation should not be performed
%% Process the files 
for ff = 1:height(optSeg)
    filepath = fullfile(optSeg.FolderPath{ff}, optSeg.FileName{ff});
    nFrames = imfinfo(filepath);
    visFrame = 1; % Select the frame for tuning segmentation parameters
    filenameOut = strrep(inputFile, ".csv", ".mat");

    threshold = segD.segThresh;
    numIteration = segD.segIteration;
    contB = segD.segContract;
    smoothF = segD.segSmooth;

    resolveFile = '';
    requiresManual = segD.requiresManual;
    manualcheck = segD.manualCheck;
    IgnoreFrames = segD.IgnoreFrames; 

    tuneParam = true;
    while tuneParam

        I = imread(filepath, visFrame);
        %
        I = imadjust(I);
        medI = medfilt2(I); 

        % Create the main figure
        mainFig = figure('Position', [100, 100, 800, 800],"Name", "KnotSegmenter");
    
        % Create the segmented axis
        segmentedAxes = axes('Parent', mainFig, 'Position', [0.05, 0.35, 0.9, 0.55]);
        title('Segmented Image'); 

        segmented = imbinarize(medI, threshold);
        actCont = activecontour(medI, segmented, numIteration, 'Chan-Vese', 'ContractionBias', contB, 'SmoothFactor', smoothF);
        overImage = imoverlay(medI, bwskel(actCont), 'red'); 
        overImage = imoverlay(overImage , bwperim(actCont), 'yellow'); 
        imshow(overImage, 'Parent', segmentedAxes);
        title(segmentedAxes, ['Segmented Image (Threshold = ' num2str(threshold) ')'  ' Frame ' num2str(visFrame)]);
   

        % Default values
        defaultFrame = visFrame;
        defaultThreshold = threshold;
        defaultContB = contB;
        defaultSmoothFact = smoothF;
        defaultNumIteration = numIteration;
    
       % Create Frame slider
        frameSlider = createSlider('Frame', 1, length(nFrames), defaultFrame,[200 50 500 20], ...
            @(src, ~) updateSegmentationMain(src, filepath, frameSlider, thresholdSlider, ...
            contBEditBox, smoothFactEditBox, numIterationEditBox, segmentedAxes));
        
        % Create Threshold slider
        thresholdSlider = createSlider('Threshold', 0, 1, defaultThreshold, [200 20 500 20], ...
            @(src, ~) updateSegmentationMain(src, filepath, frameSlider, thresholdSlider, ...
            contBEditBox, smoothFactEditBox, numIterationEditBox, segmentedAxes));
        
        % Create Edit box for contB
        contBEditBox = createEditBox('Contraction Bias (-1, 1) ', defaultContB, [320, 200, 150, 50], ...
            @(src, ~)  updateSegmentationMain(src, filepath, frameSlider, thresholdSlider, ...
            contBEditBox, smoothFactEditBox, numIterationEditBox, segmentedAxes));
        
        % Create Edit box for smoothFact
        smoothFactEditBox = createEditBox('Smooth Factor (0 , 1)', defaultSmoothFact, [320, 160, 100, 30], ...
            @(src, ~)  updateSegmentationMain(src, filepath, frameSlider, thresholdSlider, ...
            contBEditBox, smoothFactEditBox, numIterationEditBox, segmentedAxes));
        
        % Create Edit box for numIteration
        numIterationEditBox = createEditBox('Iterations (5,50))', defaultNumIteration, [320, 110, 100, 30], ...
            @(src, ~) updateSegmentationMain(src, filepath, frameSlider, thresholdSlider, ...
            contBEditBox, smoothFactEditBox, numIterationEditBox, segmentedAxes));


         % Contine box for exit 
         continueButton = uicontrol('Style', 'pushbutton', 'String', 'Accept', ...
         'Position', [600, 80, 100, 30], 'Callback', @continueScript);

         updateSliderCallbacks(frameSlider, thresholdSlider, contBEditBox, smoothFactEditBox, numIterationEditBox, segmentedAxes,filepath);
        % Pause the execution of the script until the "Continue" button is clicked
        uiwait;
        optSeg.segThresh(ff) =  thresholdSlider.Value;
        optSeg.segIteration(ff) = str2double(numIterationEditBox.String);
        optSeg.segContract(ff) =  str2double(contBEditBox.String);
        optSeg.segSmooth(ff) = str2double(smoothFactEditBox.String);
        
        % These parameters can be edited in the {inputFileName}_optimized.csv after checking the
        % segmentation output 
        optSeg.requiresManual(ff) = false;
        optSeg.IgnoreFrames(ff) =  [0]; 
        optSeg.useLogP(ff) = "False";
        optSeg.logPk(ff) = "False"; 
        optSeg.logPs(ff) = "False";
        optSeg.InvertManual(ff) = "False"; 
        optSeg.manualCheck(ff) = "False"; 
        tuneParam = false; 

    end
    close all 
end

%% Save the OPTseg file 
outPutFile = fullfile(inputPath,filename);
writetable(optSeg, outPutFile); 
disp(["Ouput file with optimized parameters is written as: " , outPutFile])

%% RUN Knot RESOLVER
% Create a yes/no dialog box
choice = questdlg('Do you want to run KnotResolver ?', 'Confirmation', 'Yes', 'No', 'No');

switch choice
    case 'Yes'
        disp('Processing...');
        % Call your processing function or script here
        KnotResolver(filename, inputPath); 
        disp('Processing completed.');
    case 'No'
        disp('Operation canceled.');
end


%% Required functions 
function continueScript(~, ~)
    % Close the current figure to resume script execution
    uiresume;
end

function updateSegmentationMain(~,filepath,frameSlider, thresholdSlider, ...
            contBEditBox, smoothFactEditBox, numIterationEditBox,segmentedAxes)
    % Retrieve values from UI controls
    frame = ceil(frameSlider.Value);
    threshold = thresholdSlider.Value;
    contB = str2double(contBEditBox.String);
    smoothFact = str2double(smoothFactEditBox.String);
    numIteration = str2double(numIterationEditBox.String);
    I = imread(filepath, frame);
    I = imadjust(I);
    medI = medfilt2(I);
    segmented = imbinarize(medI, threshold);
    actCont = activecontour(medI, segmented, numIteration, 'Chan-Vese', 'ContractionBias', contB, 'SmoothFactor', smoothFact);
    overImage = imoverlay(medI, bwskel(actCont), 'red'); 
    overImage = imoverlay(overImage , bwperim(actCont), 'yellow'); 
    imshow(overImage, 'Parent', segmentedAxes);
    title(segmentedAxes, ['Segmented Image (Threshold = ' num2str(threshold) ')' ' Frame '  num2str(frame)]);

end

function updateSliderCallbacks(frameSlider, thresholdSlider, contBEditBox, smoothFactEditBox, numIterationEditBox, segmentedAxes,filepath)
    set(frameSlider, 'Callback', @(src, ~) updateSegmentationMain(src, filepath, frameSlider, thresholdSlider, ...
        contBEditBox, smoothFactEditBox, numIterationEditBox, segmentedAxes));
    set(thresholdSlider, 'Callback', @(src, ~) updateSegmentationMain(src, filepath, frameSlider, thresholdSlider, ...
        contBEditBox, smoothFactEditBox, numIterationEditBox, segmentedAxes));
    set(contBEditBox, 'Callback', @(src, ~) updateSegmentationMain(src, filepath, frameSlider, thresholdSlider, ...
        contBEditBox, smoothFactEditBox, numIterationEditBox, segmentedAxes));
    set(smoothFactEditBox, 'Callback', @(src, ~) updateSegmentationMain(src, filepath, frameSlider, thresholdSlider, ...
        contBEditBox, smoothFactEditBox, numIterationEditBox, segmentedAxes));
    set(numIterationEditBox, 'Callback', @(src, ~) updateSegmentationMain(src, filepath, frameSlider, thresholdSlider, ...
        contBEditBox, smoothFactEditBox, numIterationEditBox, segmentedAxes));
end

% Helper function to create sliders
function slider = createSlider(label, minValue, maxValue, defaultValue, Position,callback)
    boxPosition = Position;
    boxPosition(1) = boxPosition(1) -  100; 
    boxPosition(3) = boxPosition(1) + 10;  
    boxPosition(4) = boxPosition(4) + 1 ; 
    uicontrol('Style', 'text', 'Position', boxPosition, 'String', label, 'FontSize',12);
    slider = uicontrol('Style', 'slider', 'Min', minValue, 'Max', maxValue, 'Value', defaultValue, ...
        'Position', Position, 'Callback', callback);
end

% Helper function to create edit boxes
function editBox = createEditBox(label, defaultValue,Position, callback)
    boxPosition = Position;
    boxPosition(1) = 150; 
    boxPosition(3) = boxPosition(3) + 50; 
    uicontrol('Style', 'text', 'Position',boxPosition, 'String', label, 'FontSize',12);
    editBox = uicontrol('Style', 'edit', 'Position', Position, 'String', num2str(defaultValue), ...
        'Callback', callback);
end

