%% Written by Dhruv Khatri , Cycels IISER Pune, 2022. 
% Contributions from Prince Shah, Prof. Chaitanya Athale. 
function resolveCoordinates = SingleFilamentResolver(data_struct,frameIgnore, frameEdit, frameInvert,outFolder)
%SINGLEFILAMENTRESOLVER, created by Dhruv Khatri, Cycels, IISER Pune, Oct 3 2022
% Contributions from Prince Shah & Chaitanya A Athale 
% Script performs contour resolution of "branched" microtubules in a gliding
% assay. 
% Input : A struct (1 x N), where N is number of frames containing skeleton
% indices of segmeted microtubule and information on existence of branches.
% The field include "frameNumber","imageSize", "filamentType" and "Contour". 
% "filamentType": Should be either 'Unbranched' or 'Branched', please use
% the criteion mentioned in the paper for branchpoint identification.  
% "Contour": Binary image with a single connected object representing the
% filament contour. 
% Required Packages : MATLAB 2022a or newer, Image Processing Toolbox,
% Signal Processing Toolbox. 


%message = {'Enter the list of frames to ignore from analysis (sep= ",")';
%    'Manually edit the path on (sep = ",")'};
%titlebox = 'Input'; 
%dimsbox = [1 40];

%inputData = inputdlg(message, titlebox, dimsbox); 
%frameIgnore= frameIgnore; %str2double(inputData{1});
%frameEdit = frameEdit; % str2double(split(inputData{2}, ',')); 




% Check if all fields are provided in the given struct.  
if isfield(data_struct, 'FrameNumber') && isfield(data_struct, 'FilamentType') ...
        && isfield(data_struct, 'Contour') && isfield(data_struct, 'ImageSize')
    fprintf('Starting.....')
else 
    error(['One or more of the four fields required in the input struct are missing \n' ...
        'please check help SingleSilamentResolver for more information on the input'])
end

% Check if the first frame is branched or unbranched & if the fields are
% non empty 
% remove empty data_struct 

data_struct = data_struct(~cellfun(@isempty,{data_struct.Contour})); 
%data_struct(data_struct) = []
resolveCoordinates = data_struct;
resolveCoordinates(ismember(cat(1,resolveCoordinates.FrameNumber), frameIgnore')) = []; 
if ~isempty(resolveCoordinates)
    numFilaments = length(resolveCoordinates);
    fprintf('Number of frames provided %d', numFilaments); 
    if strcmp(resolveCoordinates(1).FilamentType, 'Branched')
        error(['CAUTION!, the first input contour is branched, either provide' ...
            'an unbranched contour in the first frame or flip the time series.']); 
    end 
else 
    error('Input data struct contains empty fields')  
end 



referenceEnd = []; 

f = waitbar(0, 'Processing....');
for f = 1:numFilaments
    waitbar(f/numFilaments);
    disp(f); 
    curContour = zeros(resolveCoordinates(f).ImageSize); 
    pxIndices = resolveCoordinates(f).Contour;
    curContour(pxIndices) = 1; 
    [IsolateImage, Offset] = ImIsolate(curContour,pxIndices); 
    if isempty(referenceEnd)
        referenceEnd = resolveCoordinates(1).FixedEnd; 
        [yr,xr] = ind2sub(size(curContour), referenceEnd); 
        referenceEnd  = sub2ind(size(IsolateImage), yr-Offset(1), xr-Offset(2)); 
        prevContour = size(IsolateImage); 
    else
        prevContour = resolveCoordinates(f-1).smallSize; 
    end
    curContour = logical(IsolateImage); 
    % need to adjust the referenceEnd after Isolation 
       

    if strcmp(resolveCoordinates(f).FilamentType, 'Unbranched')
        isBranched = 0;
    else
        isBranched = 1; 
    end 
    %imageSize = size(curContour); 
    if ~isBranched
        % Sort an unbranched structure 
        E1 = find(bwmorph(curContour, 'endpoints')); 
        if length(E1) > 2
            error(['A filament with more than two end points found please ' ...
                'remove any spurious branches beforehand, frame no: %d'], f);
        end 
        [sortedSkeleton, referenceEnd] = sortSkel(curContour,E1, ...
            referenceEnd, prevContour); 
        resolveCoordinates(f).Skeleton = sortedSkeleton;
        resolveCoordinates(f).Offset = Offset;
        resolveCoordinates(f).smallSize = size(curContour); 

    else 
        % sort a branched structure 
         E1 = find(bwmorph(curContour, 'endpoints'));
         %previousContour = resolveCoordinates(f-1).Contour; 
         if length(E1) > 2 
            sprintf(['More than two ends found for a branched filament ' ...
             'code will still try to resolve, but we recommed removing them, frame: %d'], f)

         end 
         % pick a reference contour from the last resolved/unbranched
         % contour, make sure that ignored frames are not picked 
         if ismember(resolveCoordinates(f).FrameNumber, frameEdit) 
             manualEdit = true;
             if ismember(resolveCoordinates(f).FrameNumber,frameInvert)
                 frameInvertTrue = true; 
             else
                 frameInvertTrue = false; 
             end 
         else 
             manualEdit = false;
             frameInvertTrue = false; 
         end 
         refSize  = resolveCoordinates(f-1).smallSize; 
         refOffset  = resolveCoordinates(f-1).Offset; 
         [sortedSkeleton, referenceEnd, ~, ~] = sortBranch_v3(curContour, E1, ...
             resolveCoordinates(f-1).Skeleton, referenceEnd, ...
             refSize,refOffset,Offset,manualEdit,frameInvertTrue); 
         resolveCoordinates(f).Skeleton = sortedSkeleton; 
         %resolveCoordinates(f).Possbilities = allpos;
         %resolveCoordinates(f).Scores = allScore; 
         resolveCoordinates(f).Offset = Offset;
         resolveCoordinates(f).smallSize = size(curContour); 

    end
end
close all 
%% Display the resolved structure 
if ~isfolder(outFolder)
    mkdir(outFolder); 
end 

if exist(fullfile(outFolder, 'segOverlay.tif'), 'file') 
    delete(fullfile(outFolder,'segOverlay.tif'));
end


if exist(fullfile(outFolder, 'Contour.tif'), 'file') 
    delete(fullfile(outFolder,'Contour.tif'));
end

for d = 1:length(resolveCoordinates)
    curC = resolveCoordinates(d).Contour;
    ims = zeros(resolveCoordinates(d).ImageSize);
    ims(curC) = 1; 
    indxC = resolveCoordinates(d).Skeleton; 
    [yy,xx] = ind2sub(resolveCoordinates(d).smallSize, indxC); 
    offss = resolveCoordinates(d).Offset; 

    xx = xx + offss(2); 
    yy = yy + offss(1); 
    ims = ims.*255; 
    ims = repmat(ims,[1 1 3]); 
    figure(1), imshow(ims, 'border','loose'); title(resolveCoordinates(d).FrameNumber); hold on 
    %figure(1), plot(xx, yy,'r-');hold off 
    try
        figure(1),patch([xx' nan], [yy' nan],[1:length(xx) nan], 'edgecolor', 'interp', 'LineWidth', 2.0); hold off
        imwrite(getframe(gcf).cdata, fullfile(outFolder, 'segOverlay_v4.tif'),'Writemode', 'append');
    catch 
        imwrite(getframe(gcf).cdata, fullfile(outFolder, 'segOverlay_v4.tif'),'Writemode', 'append');

    end 
    %imwrite(getframe(gcf).cdata, fullfile(outFolder, 'Contour.tif'), 'Writemode', 'append');
end 
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);
save(fullfile(outFolder, "resolveCoordinates_v4.mat" ),"resolveCoordinates"); 
end 

function  [IsolateImage, Offset]= ImIsolate(Im,PixelIndices)
%IMISOLATE, creates a binary image that bounds the object to be resolved
%INPUT > requires the original image, i.e, Im & pixel indices of the object
%      > to be resolved

%OUTPUT > Gives IsolateImage with reduced image size
%       > Offset will help us to correlate back to the the original indices after being
%       resolved 


% Decoy = zeros(size(Im)); 
% Decoy(PixelIndices) = 1; 
% figure ,imshow(Decoy)
% convert index to pixel indices
[yy, xx] = ind2sub(size(Im), PixelIndices);
% these values will decide the area bounded by the cropped image 
minX = (min(xx));
minY = min(yy);
maxX = max(xx);
maxY = max(yy);
Offset = [minY-2 minX-2];

% crop the image 
IsolateImage = Im(minY-1:maxY+1, minX-1:maxX+1); 
IsolateImage = IsolateImage*.0;
ss =size(IsolateImage);
row = xx-Offset(1,2);
col=  yy - Offset(1,1);
TruePixel = sub2ind(ss, col, row); 
IsolateImage(TruePixel) = 1 ; 
% figure , imshow(IsolateImage)
end



