function [sortedSkeleton, referenceEnd] = sortSkelManual(curContour,E1, ...
                referenceEnd)

%% Set a reference point
% this is the selected fixed point which will arrange the contour indices. 
sizeImage = size(curContour); 
if isempty(referenceEnd)
    % if no  reference is provided just select the first end point
    referenceEnd = E1(1); 
else 
    % Else from the current endpoints E1, select the closest reference
    % point, i.e, referenceStart, this is the input for the next iteration
    % and is passed as the second output 
    [ye, xe] = ind2sub(sizeImage, E1); 
    [yr, xr] = ind2sub(sizeImage, referenceEnd);
    dist_matrix = pdist2([ye, xe], [yr, xr]); 
    [~, I] = min(dist_matrix, [],1);
    referenceEnd = E1(I); 
    if length(E1) == 2
        referenceStart =  setdiff(E1, referenceEnd, 'stable'); 
    end 
end 
%%  Get the graph structure
branchImage = bwlookup(curContour, makelut(@(x) sum(x(:)) >=4 & x(5) ==1,3)); % this is the 
% branch criterion ref: Mathworks
rmBranch = curContour&~branchImage; 
[s,t,pxT, pxC, label_image,maxLabelImg] = get_graph(rmBranch, branchImage); 

% Create the graph table 
G0 = graph(s,t); 
table_graph = sort(G0.Edges.EndNodes,2); 
table_unique = unique(table_graph); 

% We now create a digraph where the start and end points can be reached
% only once, all intermediate segments can be traversed in both directions
startPoint = label_image(referenceEnd); 
group = unique([s', t'], 'rows'); 
group = [group; group(:,2), group(:,1)];
group((group(:,2) == startPoint),:) = []; 

if length(E1) ==2 
    endPoint = label_image(referenceStart); 
    group((group(:,1) == endPoint),:) = []; 
else
    endPoint = [];
end 

G = digraph(group(:,1), group(:,2)); 
%figure(1), imshow(label2rgb(label_image));
filamentTable = table(t',s',pxC',pxT','VariableNames',...
    ["From","To", "FromEnd", "ToStart"]); 
%figure(100), plot(G); 
%% Get loops
Trueloop = []; 
for l = 1:length(table_unique)
    lab = table_unique(l); 
    [indx_col, ~] = find(table_graph == lab); 
    if length(indx_col) > 1
        all_values = table_graph(indx_col, :); 
        all_values = setdiff(all_values(:), lab); 
        if length(unique(all_values)) == 1
            % But this has to avoid selecting the first indice
            Trueloop  = [Trueloop, lab]; 
        end 
    end 
end
%% Add end information to filament table
% Create coordinates indices of all possibilities and return the best match
% if we have a manual edit, get the correct path input manually; 
manualEdit = 1; 
if manualEdit
    set(0,'units','pixels')
    dim_desk = get(0,'ScreenSize');
    pts_or = readPoints(label2rgb(label_image), dim_desk); 
    pts = round(pts_or'); 
    pts_2_ind = sub2ind([512, 512], pts(:,2), pts(:,1)); 
    resizeLabel = imresize(label_image, [512 512], "nearest"); 
    labels_correct = resizeLabel(pts_2_ind); 
    %possiblePaths = {}; 
    possiblePaths = labels_correct';
    % Also return the selection in original image dimension (resize them) 
    original_size = size(label_image);
    scal_col = 512 / original_size(2); 
    scal_row = 512 / original_size(1); 
    ori_pts = [round(pts_or(1,:)/scal_row)'  round(pts_or(2,:)/scal_col)']; 
    ori_pts = sub2ind([original_size(2) original_size(1)], ori_pts(:,2), ori_pts(:,1));
end 
%% Arrange coordinates
% Arrange the coordinates for all paths according to the filament table
% How to avoid redundant calculations? We will address this in Future :)
szP = size(possiblePaths); 
arrangement1 = cell(szP(1),1); 
for l = 1:szP(1)
    curPath = nonzeros(possiblePaths(l,:)); 
    %curPath = nonzeros(possiblePaths{l}); 
    contourIndice  = sortCoordinates(curPath,label_image, ...
    Trueloop,filamentTable); 
    arrangement1{l} = contourIndice; 
end 
%% Format coordinates and compute score 
% Currently we will only compute two combinates, i.e, the solution performs
% poorly when more than one loop is present, [Future Work]
sortedSkeleton = arrangement1{l};

set_coord = ori_pts(labels_correct == Trueloop);
if length(set_coord) > 1
    error("More than One loops encountered in first frame, current" + ...
        " implementation does not handles more than one loop")
end
set_coord= set_coord(1); 
seq1 = sortedSkeleton(:,1);
findRep = find(~sortedSkeleton(:,2));
seq2 = sortedSkeleton(:,2);
loopSearch = find(sortedSkeleton(:,2));
seq2(findRep) = seq1(findRep);

% Only handles single loops [FUTURE WORK] 
[x1, y1] = ind2sub([original_size(2) original_size(1)], seq1); 
[x2, y2] = ind2sub([original_size(2) original_size(1)], seq2);
[xxx, yyy] = ind2sub([original_size(2) original_size(1)], set_coord);
distances1 = pdist2([x1, y1], [xxx, yyy]);
distances2 = pdist2([x2,y2], [xxx, yyy]);

[~ ,l1] = min(distances1(loopSearch));
[~, l2] = min(distances2(loopSearch)); 

if l2 < l1
    sortedSkeleton = seq2; 
else
    sortedSkeleton = seq1; 
end 
% need to get contours for possible combinations
%% GETGRAPH
function  [T,B, pxT, pxS, LABEL_IMAGE,maxLabelImg] = get_graph(rm_branch, branch_image)
% Returns a graph representation (G) of branch_image and the labeled image
% Possible bugs: NAN
image_size = size(rm_branch); 

% label the segments, connectivity of 8 is used. we create a label
% representation of rm_branch and branch_image separatly then reassign the
% numbers
labeld_obj = bwlabel(rm_branch, 8);
maxLabelImg = max(labeld_obj(:)); 
uniqu_fils = max(labeld_obj(:)); 
branch_connected = bwconncomp(rm_branch, 8); 
label_branch_points  = bwlabel(branch_image);
label_branch_points = label_branch_points + uniqu_fils; 
label_branch_points = label_branch_points.*branch_image;
LABEL_IMAGE = labeld_obj + label_branch_points; 

T = [];
B = [];

pxT = []; 
pxS=  [];

% Create the adjacency Matrix and return s and t 
for obj = 1:branch_connected.NumObjects
    % go over each object
    null_image = zeros(size(rm_branch)); 
    curr_indice = branch_connected.PixelIdxList{obj};
    null_image(curr_indice) = 1; 
    cur_ends = find(bwmorph(null_image, 'endpoints'));
    %if length(cur_ends) > 2
        %error('More than two ends for a unique branch')
    %end 
    %curr_label = labeld_obj(curr_indice(1)); 
    for e  = 1:length(cur_ends)
        sin_end = cur_ends(e);
        vincinty_coords = get_vincity(sin_end, image_size);
        connected_by = labeld_obj(sin_end); 
        connected_too = nonzeros((label_branch_points(vincinty_coords))); 
        px_tpp = vincinty_coords(find(label_branch_points(vincinty_coords))); 
        for c = 1:length(connected_too)
            next_image = zeros(size(rm_branch)); 
            next_image(label_branch_points == connected_too(c))= 1;
            next_ends = find(bwmorph(next_image, 'endpoints'));
            %if length(next_ends) >2
                %error('More than two ends for a adjacent branch')
            %end
            % for the adjacency matrix
            T = [T,  connected_too(c)];
            B = [B, connected_by]; 

            % for the sorting
            pxT = [pxT, px_tpp(c)];
            pxS = [pxS, sin_end];
            %adjacency_matrix = [adjacency_matrix; connected_too(c), connected_by]; 
        end 
    end 
end
%G = graph(s,t); 
end 

%% GETVINCITY 
function around_coords = get_vincity(single_end, image_size)
% Return indices around a single index in teh neighbourhood of 8
% Possbile Bugs: function will through error if any indice is present on
% the image boundary, use imclearborder on the skeleton image to avoid such
% error 
[yy, xx] = ind2sub(image_size, single_end); 
around_coords = [yy-1 xx-1;
    yy xx-1;
    yy+1 xx-1;
    yy-1 xx;
    yy+1 xx;
    yy-1 xx+1;
    yy xx+1
    yy+1 xx+1];

around_coords = sub2ind(image_size, around_coords(:,1), around_coords(:,2));
% Return Indices 
end 
%% SORTCOORDINATES
function contourIndice  = sortCoordinates(selectedPath,label_image, ...
    Trueloop, filamentTable)
% SORTCOORDINATES sorts coordinates t that pass through t+1 segment. Where
% t is the current label in a selected Path.  
% Core idea is to use the seed point of t contour and the end point of t+1
% contour to find the shortest path between t & t +1 through t. The indice
% that is visited first in t+1 forms the seed point for the next iteration. 
%tempPath = selectedPath;
%loopCrit = sum(selectedPath == Trueloop,2); 
imageSize = size(label_image); 

fromTooMat = filamentTable{:, ["From","To"]}; 
matSize = size(fromTooMat); 

contourIndice = zeros(0,0); 
for o = 1:length(selectedPath)-1
    spair = selectedPath(o:o+1); 
    loopCrit = any(ismember(spair,Trueloop)); 
    % Find the pair entry in the table 
    indxT = []; 
    for p = 1:matSize(1)
        if fromTooMat(p,:) == sort(spair)'
            indxT = [indxT, p]; 
        end 
    end 
    % Three cases, either one true, multiple true, non true 
    if ~loopCrit && length(indxT) == 1
        % initiate image 
        null_image = zeros(size(label_image)); 
        null_image(label_image == spair(1)) = 1;
        table_row = fromTooMat(indxT,:);
        val_f = find(table_row == spair(1)); 

        % If contour exits use the last point 
        if length(contourIndice) >= 1
            checkPreviousLoop = length(unique(nonzeros(contourIndice(end, :))))==2; 
            if checkPreviousLoop
                bothEnds = contourIndice(end,:); 
                possi1 = []; 
                possi2 = []; 

                for k = 1:2
                    seedPoint = bothEnds(k); 
                    if val_f == 1
                       fstart = filamentTable{indxT,"ToStart"}; 
                    else
                       fstart = filamentTable{indxT, "FromEnd"}; 
                    end
                    null_image(fstart) = 1; 
                    shortestPath = findShortestPath(seedPoint, fstart, imageSize, null_image);
                    if k == 1
                        possi1 = shortestPath;
                    else 
                        possi2 = shortestPath;
                    end
                end 

                while length(possi2) ~= length(possi1)
                    if length(possi1) > length(possi2)
                       possi2 = [possi2; possi2(end)];
                    else
                       possi1 = [possi1; possi1(end)];
                    end
                end 
                tempS = [possi1, possi2]; 
                contourIndice = [contourIndice; tempS];
            else
                seedPoint = contourIndice(end,1);
                % in case we are coming from loop region. 
                if val_f == 1
                    fstart = filamentTable{indxT, "ToStart"};
                else
                    fstart = filamentTable{indxT, "FromEnd"}; 
                end
                   % end point 
                null_image(fstart) = 1; 
                shortestPath = findShortestPath(seedPoint, fstart, imageSize, null_image);
                contourIndice = [contourIndice; shortestPath, zeros(length(shortestPath),1)]; 
            end 
        else
            E1 = find(bwmorph(null_image, 'endpoints'));
            % get row
            if val_f == 1
                fend = filamentTable{indxT, "FromEnd"};
                fstart = filamentTable{indxT, "ToStart"};
                seedPoint =setdiff(E1, fend);
            else
                fend = filamentTable{indxT, "ToStart"};
                fstart = filamentTable{indxT,"FromEnd"};
                seedPoint = setdiff(E1,fend);
            end
            if isempty(seedPoint)
                % empty when there is only one indice
                contourIndice = [contourIndice; fend, zeros(length(fend),1)];
            else
                null_image(fstart) = 1;
                shortestPath = findShortestPath(seedPoint, fstart, imageSize, null_image);
                contourIndice = [contourIndice; shortestPath, zeros(length(shortestPath),1)];

            end
            % end point
        end
    else
        % IF we are going towards or coming from a loop structure
        direc = ismember(spair(1),Trueloop);
        table_row = fromTooMat(indxT,:);


        if length(contourIndice) < 1

            null_image = zeros(size(label_image));
            null_image(label_image == spair(1)) = 1;
            %E1 = find(bwmorph(null_image, 'endpoints')); 
            table_row = fromTooMat(indxT,:);
            val_f = find(table_row(1,:) == spair(1));

            if val_f == 1
                % NEEED TO FIX 
                fend = filamentTable{indxT, "FromEnd"};
                fstart = filamentTable{indxT, "ToStart"};
                seedPoint =filamentTable{indxT, "FromStart"}; 
            else 
                fend = filamentTable{indxT, "ToStart"};
                fstart = filamentTable{indxT,"FromEnd"};
                seedPoint = filamentTable{indxT, "ToEnd"};
            end
            
            if isempty(seedPoint)
                seedPoint = E1; 
            end

            possi1 = [];
            possi2 = [];
            
            % We are going towards a loop
            for u = 1:length(indxT)
                null_image = zeros(size(label_image));
                null_image(label_image == spair(1)) = 1;
                val_f = find(table_row(u,:) == spair(1));
                if val_f == 1
                    fstart = filamentTable{indxT(u), "ToStart"};
                else
                    fstart = filamentTable{indxT(u), "FromEnd"};
                end
                null_image(fstart) =1;
                shortestPath = findShortestPath(seedPoint, fstart, imageSize, null_image);
                if u == 1
                    possi1 = shortestPath;
                else
                    possi2 = shortestPath;
                end
            end
            while length(possi2) ~= length(possi1)
                if length(possi1) > length(possi2)
                    possi2 = [possi2; possi2(end)];
                else
                    possi1 = [possi1; possi1(end)];
                end
            end
            tempS = [possi1, possi2];
            contourIndice = [contourIndice; tempS];

        else
            if ~direc
                % GOing towards the loop 
                possi1 = [];
                possi2 = [];
                for u = 1:length(indxT)
                    null_image = zeros(size(label_image));
                    null_image(label_image == spair(1)) = 1;
                    val_f = find(table_row(u,:) == spair(1));
                    seedPoint = contourIndice(end,1);
                    if val_f == 1
                        fstart = filamentTable{indxT(u), "ToStart"};
                    else
                        fstart = filamentTable{indxT(u), "FromEnd"};
                    end
                    null_image(fstart) =1;
                    shortestPath = findShortestPath(seedPoint, fstart, imageSize, null_image);
                    if u == 1
                        possi1 = shortestPath;
                    else
                        possi2 = shortestPath;
                    end
                end
                while length(possi2) ~= length(possi1)
                    if length(possi1) > length(possi2)
                        possi2 = [possi2; possi2(end)];
                    else
                        possi1 = [possi1; possi1(end)];
                    end
                end
                tempS = [possi1, possi2];
                contourIndice = [contourIndice; tempS];
            else
                % coming from a loop 
                possi1 = [];
                possi2 = [];
                
                val_f = find(table_row(1,:) == spair(1));

                if val_f == 1
                    bothPoints = filamentTable{indxT, "FromEnd"}; 
                    bothAppend = filamentTable{indxT, "ToStart"}; 
                else
                    bothPoints = filamentTable{indxT, "ToStart"};
                    bothAppend = filamentTable{indxT, "FromEnd"}; 
                end

               
                for u = 1:length(indxT)
                    null_image = zeros(size(label_image));
                    null_image(label_image == spair(1)) = 1;
                    val_f = find(table_row(u,:) == spair(1));
                    if u == 1
                        seedPoint = contourIndice(end,1); 
                        fstart = setdiff(bothPoints, seedPoint); 
                        burnI = bothAppend(2); 
                    else 
                        seedPoint = contourIndice(end,2); 
                        fstart = setdiff(bothPoints, seedPoint); 
                        burnI = bothAppend(1); 
                    end
                  
                    null_image(fstart) =1;
                    shortestPath = findShortestPath(seedPoint, fstart, imageSize, null_image);
                    shortestPath = [shortestPath; burnI];
                    if u == 1
                        possi1 = shortestPath;
                    else
                        possi2 = shortestPath;
                    end
                end
                while length(possi2) ~= length(possi1)
                    if length(possi1) > length(possi2)
                        possi2 = [possi2; possi2(end)];
                    else
                        possi1 = [possi1; possi1(end)];
                    end
                end
                tempS = [possi1, possi2];
                contourIndice = [contourIndice; tempS];
            end
        end

    end
end
% Final arrangement of the contour
null_image= zeros(size(label_image));
null_image(label_image==selectedPath(end))= 1;
E1 = find(bwmorph(null_image, 'endpoints'));
finalPoint = setdiff(E1, contourIndice(end,1));
null_image(contourIndice(end,1))=1; 
shortestPath = findShortestPath(contourIndice(end,1),finalPoint,imageSize,null_image); 
if ismember(selectedPath(end), Trueloop)
    contourIndice = [contourIndice; shortestPath, zeros(length(shortestPath),1)];
    contourIndice(end-length(shortestPath)+1:end,2) = flip(shortestPath); 
else
    contourIndice = [contourIndice; shortestPath, zeros(length(shortestPath),1)];

end

%dropDup = []; 
%for c = 1:length(contourIndice(:,1))-1
%if contourIndice(c) == contourIndice(c+1)
%    dropDup = [dropDup, c]; 
%end 
%end 
%contourIndice(dropDup) =[]; 
end 
 


%% FINDSHORTESETPATH
function shortestPath = findShortestPath(seedPoint,E1, imageSize, null_image)
% FINDSHORTESTPATH finds shortest path between seedPoint and E1 in null_image
% The shortest path starts at the seed point and ends at E1, this orients
% the contour label belonging to the seed point in ascending order. This
% function was implemented after reading path finding threads on MATHWORKS website,
% all credits to online resources. Funtion does not works with looped
% structure, therefore null_image needs to be preprocessed to avoid
% returing Inf values from the distance transform. 
% Ref: Mathworks maze shortest path tutorial
% https://blogs.mathworks.com/steve/2011/12/02/exploring-shortest-paths-part-3/


    [rr1 ,cc1] = ind2sub(imageSize, seedPoint); 
    [rr2, cc2] = ind2sub(imageSize, E1); 
    D1 = bwdistgeodesic(logical(null_image), cc1, rr1, 'quasi-euclidean'); 
    D2 = bwdistgeodesic(logical(null_image), cc2, rr2, 'quasi-euclidean'); 
    
    D = D1+ D2; % add both distance transforms
    D = round(D * 8) / 8;
    D(isnan(D))  = inf; % Assign nan to inf value
    paths = imregionalmin(D); % Find region minima 
    
    coordinates_short_path = find(paths); 
    distnce_path = D1(coordinates_short_path); % sort the coordinates
    [~, I] = sort(distnce_path); 
    
    coordinates_short_path = coordinates_short_path(I); 
    shortestPath = coordinates_short_path; % return sorted coordinates
end 
%% READPOINTS for manual input
function pts = readPoints(image, dim_desk)
% TAKEN FROM : https://in.mathworks.com/matlabcentral/answers/118724-how-to-get-onclick-coordinate-pixel-value-and-location-from-an-image
%readPoints   Read manually-defined points from image
%   POINTS = READPOINTS(IMAGE) displays the image in the current figure,
%   then records the position of each click of button 1 of the mouse in the
%   figure, and stops when another button is clicked. The track of points
%   is drawn as it goes along. The result is a 2 x NPOINTS matrix; each
%   column is [X; Y] for one point.
% 
%   POINTS = READPOINTS(IMAGE, N) reads up to N points only.
if nargin < 4
    n = Inf;
    pts = zeros(2, 0);
else
    pts = zeros(2, n);
end

pause(1.0)
% display image
xold = 0;
yold = 0;
k = 0;
fig = figure(1), imshow(imresize(image, [512, 512], "nearest"), 'InitialMagnification', 700, "Border", "tight"); hold on

title( 'Perform manual selection of the correct contour path' );

set(fig, 'Position',  [100, 100, dim_desk(3)/2, dim_desk(4)/2]) 
while 1

    [xi, yi, but] = ginput(1);      % get a point
    if ~isequal(but, 1)             % stop if not button 1
        break
    end
    k = k + 1;
    pts(1,k) = xi;
    pts(2,k) = yi;
      if xold
          plot([xold xi], [yold yi], 'go-');  % draw as we go
      else
          plot(xi, yi, 'go');         % first point on its own
      end
      if isequal(k, n)
          break
      end
      xold = xi;
      yold = yi;
  end
hold off;
if k < size(pts,2)
    pts = pts(:, 1:k);
end
end

end