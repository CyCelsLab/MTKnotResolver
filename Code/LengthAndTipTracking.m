%% Created by Dhruv Khatri ,IISER Pune 
%{
Function reads in a csv file and gets the length of the filament in the
selected frame. Also detects the tips and gets the tip coordinates of the
beating filament. 
%} 
%% Read in the Data 
% Specify the path to your CSV file
pathFile  = '/home/dhruv/Documents/Shivani/DifferentDensity/beating_newdyn_analysis/DataSheet.csv'; 
data =readtable(pathFile, "NumHeaderLines",0);
%% Iterate over each csv file 
% Iterate over each row in the CSV file
for i = 1:size(data, 1)
    % Extract file path, folder name, and file name
    filePath = char(data{i, 'directory'});
    folderName = char(data{i, 'folder'});
    fileName = char(data{i, 'filename'});
    lengthCalon = data{i,'frame'};
    % Create a outputFolder 
  
    % Combine the directory, folder, and file name to get the full file path
    fullFilePath = fullfile(filePath, folderName, fileName);
    filenameOut  = replace(fullFilePath, '.tif', '.mat') ; 
    createOutpulFolder =replace(fullFilePath, '.tif', '');
    [status,msg]=mkdir(createOutpulFolder); 
    % Load the file or perform any necessary preprocessing
    % Assuming you have a function to load the data or you can use imread
    %image = imread(fullFilePath);
    
    % Get the number of frames in the file
    numFrames = length(imfinfo(fullFilePath)); % Assuming the third dimension represents frames
    selectedFilament = false; 
    % Iterate over each frame and perform segmentation
    % data strucutre 
    dataFile = struct; 
    for frame = 1:numFrames
        frameData = imread(fullFilePath, frame); % Extract the current frame
        %frameData = imadjust(frameData); 
        %frameData =im2uint8(frameData); 
        % Perform segmentation using the MTseg function
        [segmentedImage, selectedFilament] = MTseg(frameData,selectedFilament); % Replace with your segmentation code  

        % get end points and fit the gaussianFitting function 
        tipCoords = findTips(segmentedImage,frameData); 
        if frame == lengthCalon 
            try
            getLength = filamentLength(segmentedImage); 
            catch
                getLength = 0 ; 
            end 
        end
        dataFile(frame).Tip = tipCoords; 
        
    end
    dataFile(1).Length = getLength*0.106;

    save(fullfile(createOutpulFolder, 'outstats') , "dataFile")
end

%% Defatult test 
% [X,Y] = meshgrid(-8.5:8.5);
% MdataSize = 17;
% xdata = zeros(size(X,1),size(Y,2),2);
% GTpara=  [0, 0,2.5,212,2.14]; 
% defaultImage = FilamentTipFunction(GTpara, xdata); 
% figure(1), subplot(1,4,1) , imshow(defaultImage,[]); title('Original Image') ; hold on 
% plot(GTpara(1)+8.5, GTpara(2)+8.5, 'ro'); hold off 
% 
% addNoise = defaultImage + rand(size(defaultImage,1),size(defaultImage,2))*50; 
% 
% 
% x0 = [-1, 1,1.5,150,1.14];
% lb = [-MdataSize/2,-MdataSize/2,0,0,0];
% ub = [MdataSize/2,(MdataSize/2),MdataSize/2,256,pi];
% figure(1), subplot(1,4,3),  imshow(FilamentTipFunction(x0,xdata), []); title('Initial Guess'); 
% figure(1),subplot(1,4,2),  imshow(addNoise, []); title('Noisy Image'); 
% hold on; plot(x0(1)+8.5, x0(2)+8.5, 'ro'); hold off 
% options = optimoptions('lsqcurvefit', 'FunctionTolerance',1e-7);
% [x,resnorm,residual,exitflag] = lsqcurvefit(@FilamentTipFunction,x0,xdata,addNoise,lb,ub,options);
% 
% figure(1),subplot(1,4,4),  imshow(FilamentTipFunction(x,xdata), []); title('Final Fit'); 
% hold on; plot(x(1)+8.5, x(2)+8.5, 'ro'); hold off 

%% Required Functions 
function [nullImage , selectedFilament] = MTseg(frame,selectedFilament)
% Performs mtsegmentation 
I = frame; 
medI = medfilt2(I);
% Threshold and contour optimization
BW1 = imbinarize(medI,'global');
actCont = activecontour(medI, BW1, 100, 'edge', ...
    'ContractionBias',-0.1, 'SmoothFactor',0.5);

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
end 


function  tipCoords=   findTips(segImage,frameIntensity)
segImage = logical(segImage); 
skelImage = bwskel(segImage);%, "skeleton",Inf); 
endpoints = find(bwmorph(skelImage, 'endpoints')); 
tipCoords = []; 
for e = 1:length(endpoints)
    [rr,cc] = ind2sub(size(segImage), endpoints(e)); 
    cropimage = imcrop(frameIntensity , [cc-8, rr-8, 17 ,17 ]); 
    if size(cropimage,1) ~= 18 || size(cropimage,2) ~= 18
        cropimage = padarray(cropimage, [18 - size(cropimage,1), abs(18 - size(cropimage,2))], 0,'pre'); 
    end 
    cropimage = double(cropimage); 
    [X,Y] = meshgrid(-8.5:8.5);
    xdata = zeros(size(X,1),size(Y,2),2);
    xdata(:,:,1) = X;
    xdata(:,:,2) = Y;
    x0 = [0, 0,5,212,2.14];
    lb = [-8.5,-8.5,0,0,-pi/4];
    ub = [8.5,8.5,8.5,256,pi];
    options = optimoptions('lsqcurvefit', 'FunctionTolerance',1e-7);
    [x,~,~,~] = lsqcurvefit(@FilamentTipFunction,x0,xdata,cropimage,lb,ub,options);
    
    coordx = x(1) + 9.5 + cc-9;  
    coordy = x(2) + 9.5 + rr-9; 

% 
%     figure(1), subplot(1,2,1), imshow(FilamentTipFunction(x,xdata), []); hold on
%     plot(x(1)+9.5, x(2)+9.5, 'ro'); hold off
%     figure(1), subplot(1,2,2), imshow(cropimage, []);
% 
% 
%     figure(1), imshow(frameIntensity, []); hold on 
%     plot(coordx, coordy, 'ro'); 
tipCoords = [tipCoords; coordx, coordy]; 
end 
end 


function  lengthEstimate = filamentLength(segImage) 
segImage = logical(segImage); 
skelImage = bwskel(segImage);%, "skeleton",Inf); 
endpoints = find(bwmorph(skelImage, 'endpoints'));
[rr,cc] = ind2sub(size(segImage), endpoints); 
endpoints = [cc,rr]; 
centroid = regionprops(bwconncomp(segImage), 'Centroid').Centroid; 
[max_point1, max_point2] = find_two_farthest(endpoints, centroid); 
D = bwdistgeodesic(skelImage,max_point1(1), max_point2(2), 'quasi-euclidean'); 
coordsFilament = find(skelImage); 
arrangeFil = D(coordsFilament); 
[~,I] = sort(arrangeFil); 
coordsFilament = coordsFilament(I); 
[rr,cc] = ind2sub(size(segImage), coordsFilament); 
try 
    rr = smooth(rr,0.2,'sgolay'); 
    cc = smooth(cc, 0.2, 'sgolay'); 
catch

end 
% figure(1), imshow(segImage); hold on ;
% plot(cc,rr, 'r--'); 
tip1Distance = sqrt(diff(rr).^2 + diff(cc).^2);

pathLength1 = sum(tip1Distance);
lengthEstimate = pathLength1; 
end 

function d = distance_squared(x1, y1, x2, y2)
    d = (x1 - x2)^2 + (y1 - y2)^2;
end

function [max_point1, max_point2] = find_two_farthest(points, centroid)
    max_heap = cell(size(points, 1), 1);

    for i = 1:size(points, 1)
        dist = distance_squared(points(i, 1), points(i, 2), centroid(1), centroid(2));
        max_heap{i} = {-dist, points(i, :)};
    end

    [~, idx] = max(cellfun(@(x) x{1}, max_heap));
    max_point1 = max_heap{idx}{2};
    
    % Remove the selected point from the heap
    max_heap(idx) = [];
    
    max_heap_points = cell2mat(cellfun(@(x) x{2}, max_heap, 'UniformOutput', false));
    max_point2 = find_farthest_points(max_heap_points, max_point1);
end

function max_point = find_farthest_points(points, reference_point)
    max_distance = -Inf;
    max_point = [];

    for i = 1:size(points, 1)
        dist = distance_squared(points(i, 1), points(i, 2), reference_point(1), reference_point(2));
        if dist > max_distance
            max_distance = dist;
            max_point = points(i, :);
        end
    end
end
% % Load your image (replace 'your_image.png' with your image file)
% image = imread('/home/dhruv/Documents/Dhruv/3.Analysis/Branch_resolving/FilResolver/MTLoopResolve/InputFiles/tipDetectionImage.png');
% image = double(image); 
% image = imresize(image,[16 16]); 
% figure(1), imshow(image,[], 'Border', 'tight')
% %% Plot the spot locations on the original image 
% Z = image; 
% MdataSize = 15; 
% [X,Y] = meshgrid(-MdataSize/2:MdataSize/2);
% xdata = zeros(size(X,1),size(Y,2),2);
% xdata(:,:,1) = X;
% xdata(:,:,2) = Y;
% x0 = [0, 0,5,212,2.14];
% lb = [-MdataSize/2,-MdataSize/2,0,0,-pi/4];
% ub = [MdataSize/2,(MdataSize/2),MdataSize/2,256,pi];
% options = optimoptions('lsqcurvefit', 'FunctionTolerance',1e-7);
% [x,resnorm,residual,exitflag] = lsqcurvefit(@FilamentTipFunction,x0,xdata,Z,lb,ub,options);
% %%
% figure(1), subplot(1,2,1), 
% FilamentTip(x(1),x(2),x(3),x(4),x(5))
% figure(1), subplot(1,2,2),
% imshow(Z, []); hold on 
% plot(x(1)+MdataSize/2+1,x(2)+MdataSize/2+1, 'ro'); hold off

%% Required Functions
function FilamentTip(x,y,sigma,H,theta)
[XG,YG] = meshgrid(-7.5:7.5);
f = -( XG - x ) * sin(theta-pi/2) + ( YG - y ) * cos(theta-pi/2) + 0.5;
f( f < 0 ) = 0;
f( f > 1 ) = 1;
A = 1/(2*sigma^2);
I = H*( f .*exp( -A*( ( XG - x ).^2 + ( YG - y ).^2 ) )+...
(1-f).*exp( -A*(-( XG - x )*sin(theta) + ( YG - y )*cos(theta)).^2) );
surface(XG,YG,I);
end
%%
function I = FilamentTipFunction(x0, xdata)
% X0: x1, y1, theta, H
sigma = x0(3); 
x = x0(1); 
y = x0(2); 
H= x0(4); 
theta = x0(5); 

XG =  xdata(:,:,1);
YG = xdata(:,:,2); 
f = -( XG - x ) * sin(theta-pi/2) + ( YG - y ) * cos(theta-pi/2) + 0.5;
f( f < 0 ) = 0;
f( f > 1 ) = 1;
A = 1/(2*sigma^2);
I = H*( f .*exp( -A*( ( XG - x ).^2 + ( YG - y ).^2 ) )+...
(1-f).*exp( -A*(-( XG - x )*sin(theta) + ( YG - y )*cos(theta)).^2) );
end 