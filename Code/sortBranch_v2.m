%% Written by Dhruv Khatri, IISER Pune. 
% Contributions from Prince Shah, Prof. Chaitanya Athale, IISER Pune,
% Cycels Lab. 

%% SORTBRANCH
function [sortedSkeleton, referenceEnd, allpossibilites,allScore] = sortBranch_v2(curContour, E1, previousSkeleton, ...
    referenceEnd,refSize,refOffset,Offset,manualEdit)
%SORTBRANCH find the most optimum arrangement of pixel indices based on a
% comparison with previous contour. 
% The functions are arranged in order getGraph->getPaths->selectBest.
% getGraph requires function get_vincity. getPaths requires inbuilt
% allpaths function that was introduced in 2022a ver. selectBest requires
% function sortCoordinates, which is the core function of this script. 
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
    [yr, xr] = ind2sub(refSize, referenceEnd);
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
[s,t,pxT, pxC, label_image] = get_graph(rmBranch, branchImage); 

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
%% Create all possible paths between start and endpoint 
possiblePaths = getPaths(G, startPoint, endPoint, group);
%disp(possiblePaths)
% Create coordinates indices of all possibilities and return the best match
% if we have a manual edit, get the correct path input manually; 
if manualEdit
    set(0,'units','pixels')
    dim_desk = get(0,'ScreenSize');
    pts = readPoints(label2rgb(label_image), dim_desk); 
    pts = round(pts'); 
    pts_2_ind = sub2ind(size(curContour), pts(:,2), pts(:,1)); 
    labels_correct = label_image(pts_2_ind); 
    %possiblePaths = {}; 
    possiblePaths = labels_correct';  % 
end 
%% Arrange coordinates
% Arrange the coordinates for all paths according to the filament table
% How to avoid redundant calculations? We will address this in Future :)
szP = size(possiblePaths); 
arrangement1 = cell(szP(1),1); 
for l = 1:szP(1)
    curPath = nonzeros(possiblePaths(l,:)); 
    contourIndice  = sortCoordinates(curPath,label_image, ...
    Trueloop,filamentTable); 
    arrangement1{l} = contourIndice; 
end 

%% Format coordinates and compute score 
% Currently we will only compute two combinates, i.e, the solution performs
% poorly when more than one loop is present, [Future Work]
scoreArray = zeros(length(arrangement1),1) ; 
orientation = zeros(length(arrangement1),1); 
for k = 1:length(arrangement1)
    skelCoord = arrangement1{k}; 
    loopr = skelCoord(:,2)>0;
    % might be useful for future 
    startIndx = [];
    endIndx = [];
    for w = 1 :length(loopr)-1
        if isequal(loopr(w:w+1), [0;1])
            startIndx = [startIndx; w]; 
        elseif isequal(loopr(w:w+1), [1;0])
            endIndx = [endIndx, w]; 
        end
    end  
    if ~isempty(startIndx)
        seq1 = skelCoord(:,1); 
        findRep = find(~skelCoord(:,2)); 
        seq2 = skelCoord(:,2);
        seq2(findRep) = seq1(findRep); 

        score1 = computeScore(previousSkeleton,refSize,refOffset, seq1,sizeImage,Offset); 
        score2 = computeScore(previousSkeleton, refSize, refOffset, seq2,sizeImage, Offset); 
        % Compute DTW Score
        if score2 > score1    
            % should we overwrite the arrangement here? 
           scoreArray(k) = score1; 
           orientation(k) = 1; 
           arrangement1{k} = seq1; 
           if manualEdit
               scoreArray(k) = score2; 
               orientation(k) = 2; 
               arrangement1{k} = seq2; 
           end
        else 
            scoreArray(k) = score2;
            orientation(k) = 2; 
            arrangement1{k} = seq2; 
            if manualEdit
               scoreArray(k) = score1; 
               orientation(k) = 1; 
               arrangement1{k} = seq1; 
            end 
        end 
    else 
        seq1 = skelCoord(:,1);
        arrangement1{k} = seq1; 
        % Compute DTW score 
        score1 = computeScore(previousSkeleton, refSize, refOffset, seq1,sizeImage, Offset); 
        scoreArray(k) =score1;
        orientation(k) = 1;
    end 
end 

%% Select best  [sortedSkeleton, referenceStart, allpossibilites,allScore]
[~, Indm] = min(scoreArray); 
sortedSkeleton = arrangement1{Indm}; 
allScore = scoreArray; 
allpossibilites = arrangement1; 

%figure(1), imshow(double(repmat(curContour,1,1,3))); hold on 
%[yy,xx]= ind2sub(size(label_image), sortedSkeleton);
%plot(xx,yy, 'r--'); hold off  
%pause(0.5)

%patch([xx' nan], [yy' nan], [(1:length(yy))*0.5 nan], ...
%'edgecolor', 'interp', 'LineWidth',4.0); hold off 
%axis on
%set(gca, 'FontSize', 16);
%pause(0.5)

%[sortedSkeleton, allpossibilites,allScore] = selectBest(label_image, possiblePaths, ...
%previousSkeleton,referenceStart,refSize,refOffset,Offset,Trueloop); 
end

%% GETGRAPH
function  [T,B, pxT, pxS, LABEL_IMAGE] = get_graph(rm_branch, branch_image)
% Returns a graph representation (G) of branch_image and the labeled image
% Possible bugs: NAN
image_size = size(rm_branch); 

% label the segments, connectivity of 8 is used. we create a label
% representation of rm_branch and branch_image separatly then reassign the
% numbers
labeld_obj = bwlabel(rm_branch, 8);
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

%% GETPATHS
function possiblePaths = getPaths(G,startPoint, endPoint,group)
% Function getPaths finds all possible paths between a startpoint and
% endpoint for a given graph component
% Possible bugs: NAN
% We need to improve this section. 
uniq_values = unique(group);
uniq_values = setdiff(uniq_values, [startPoint, endPoint]); 
possiblePaths = {}; 
i = 1 ; 
for p  = 1:length(uniq_values)
    path1 = allpaths(G,startPoint, uniq_values(p));
    if ~isempty(endPoint)
        path2 = allpaths(G,uniq_values(p), endPoint);

        for p1 = 1:length(path1)
            possiblePaths{i} = path1{p1}; % Feb/09/2023
            i = i +1 ; 
            for p2 = 1:length(path2)
                addPath = path2{p2}; 
                possiblePaths{i} = [path1{p1}, addPath(2:end)]; %#ok<AGROW> 
                %possiblePaths{i} = path1{p1};%, addPath(2:end)]; %#ok<AGROW> 

                i = i + 1 ; 
            end 
        end 
    else 
        for p1 = 1:length(path1)
            addPath = path1{p1}; 
            if addPath(end-1) ~= addPath(1) 
                %possiblePaths{i} = [path1{p1}, addPath(end-1)]; %#ok<AGROW>
                possiblePaths{i} = path1{p1};%, addPath(end-1)]; %#ok<AGROW>

                i = i + 1;
            end 
        end 
    end
end 

% Remove duplicates 
maxPath = max(cellfun(@length, possiblePaths)); 
resizePossible = cellfun(@(x) [x zeros(1, maxPath-numel(x))], possiblePaths,'UniformOutput',false);
resizePossible = vertcat(resizePossible{:});
[~,jj] = unique(resizePossible, 'rows'); 
possiblePaths = resizePossible(jj,:); 
end 

%% COMPUTESCORE 
function scoreM = computeScore(pContour,pSize,pOffset, cContour,cSize,cOffset)
% Clean the contour data remove repeating enteries
remI =zeros(length(cContour),1);  
for t  = 1:length(cContour)-1
    if cContour(t) == cContour(t+1)
        remI(t) = t;  
    end 
end 
cContour(nonzeros(remI)) =[]; 


[yref, xref] = ind2sub(pSize, pContour); 
yref = yref + pOffset(1);
xref = xref + pOffset(2); 

[y1,x1] = ind2sub(cSize, cContour);
y1 = y1 + cOffset(1); 
x1 = x1 + cOffset(2); 

scoreM = dtw([x1';y1'], [xref';yref']); %score1x+ score1y; 

end
%% SLECTBEST
function [bestMatch, allStore,allScore] = selectBest(label_image, possiblePaths, ...
    previousSkeleton,referenceStart,refSize,refOffset,Offset,Trueloop)
% SELECTBEST selects best matching template among all possbilepaths
% Function creates calculated dtw distance metric between the selected path
% contour list and the previous reference skeleton. We require a labelled
% image with differnet connected components (label_image). We also require prior
% knowledge of existence of any looped structure. 

[yref, xref] = ind2sub(refSize, previousSkeleton); 
yref = yref + refOffset(1);
xref = xref + refOffset(2); 
% Assign (x(1),y(1)) == (0,0)
%yref = yref;% -yref(1);
%xref = xref;% -xref(1); 
% Need to calculate an offset diff 
%adjustY = Offset(1) - refSize(1); 
%adjustX = Offset(2) - refSize(2);  

coordStore = {}; 
allStore= {}; 
counter = 1; 
score = {}; 
allScore = {}; 



imageSize = size(label_image); 
numPossi = size(possiblePaths); 

for p = 1:numPossi(1)
    % Iterate over possible paths 
    curPath =  possiblePaths(p,:); 
    curPath = curPath(curPath~=0); 

    arrangeCoordinates = sortCoordinates(curPath, label_image, ...
        referenceStart,Trueloop); % Most Imp. function, look for bugs here
    
    sizeArray = cellfun(@size, arrangeCoordinates,'UniformOutput', false); 
    
    % sortCoordinates returns arranged indices as they appear along a path,
    % for looped filament it return a Nx2 matrix instead of a single column
    % array 

    % We therefore create two possible paths considering the loop
    % orientation, CAUTION: if there are more than two loops within a
    % single skeleton, this code will return suboptimal result. 

    isLooped = false; 
    allPathContour = {}; 
    allPathContour{1} = []; 
    allPathContour{2} = [];
    for c = 1:length(sizeArray)
        curPathA =  arrangeCoordinates{c}; 
        sizePath = size(curPathA);
        if sizePath(2) > 1
           allPathContour{2} = [allPathContour{2}; nonzeros(curPathA(:,2))];
           allPathContour{1} = [allPathContour{1}; nonzeros(curPathA(:,1))]; 
           isLooped = true; 
        else
           allPathContour{1} =  [allPathContour{1}; nonzeros(curPathA)]; 
           allPathContour{2} = [allPathContour{2};nonzeros(curPathA)];
        end 
    end 

    % Find the distance score from a reference skeleton using dynamic time
    % wrapping function 

    debugg=false; 
    if debugg
        figure(10), subplot(2,2,1), imshow(label_image,[]); 
        figure(10), subplot(2,2,2), plot(xref,yref,'k--'); 
        figure(10), subplot(2,2,3), plot(xref,yref,'k--'); 
    end 

    if isLooped
        [y1,x1] = ind2sub(imageSize, allPathContour{1});
        y1 = y1 + Offset(1); 
        x1 = x1 + Offset(2); 
        %y1 = y1 - y1(1);
        %x1 = x1 - x1(1); 
        %score1 = dtw([x1, y1], [xref, yref]);
        %score1y = dtw(y1, yref); 
        %score1x = dtw(x1,xref);     
        score1 = dtw([x1';y1'], [xref';yref']); %score1x+ score1y; 

        [y2,x2] = ind2sub(imageSize, allPathContour{2}); 
        y2 = y2 + Offset(1); 
        x2 = x2 + Offset(2); 
        %y2 = y2 - y2(1);
        %x2 = x2 - x2(1); 
        %score2y = dtw(y2, yref);
        %score2x = dtw(x2, xref); 
        score2 = dtw([x2'; y2'], [xref';yref']);%score2y + score2x; 

        if debugg
            %linewidth1 = 5:numel(x1)+4; 
            %linewidth2 = 5:numel(x2)+4;
            figure(10), subplot(2,2,2); hold on, plot(x1, y1, 'r--'); 
            title(['DTW score: ' , num2str(score1)]); hold off
            figure(10), subplot(2,2,3); hold on, plot(x2,y2, 'b--'); 
            title(['DTW score : ', num2str(score2)]);  hold off
        end 

        if score1 < score2
            fScore = score1; 
            curPathA = allPathContour{1};
            
        else
            fScore  =  score2; 
            curPathA= allPathContour{2};
        end 
        allStore{counter} = allPathContour{1}; 
        allScore{counter} = score1; 
        counter = counter + 1;

        allStore{counter} = allPathContour{2}; 
        allScore{counter} = score2; 

        counter = counter + 1;

    else
        [y1,x1] = ind2sub(imageSize, allPathContour{1});
        %y1 = y1 -adjustY; 
        y1 = y1 + Offset(1); 
        x1 = x1 + Offset(2); 
        %y1 = y1 - y1(1);
        %x1 = x1 - x1(1); 
        %scorex = dtw(x1,xref); 
        %scorey = dtw(y1,yref); 
        fScore = dtw([x1';y1'], [xref';yref']);% scorey+scorex; 
        if debugg
            figure(10), subplot(2,2,2); hold on, plot(x1,y1,'b-');
            title(['DTW score: ' , num2str(fScore)]); hold off
        end 
        curPathA = allPathContour{1}; 
        allStore{counter} = curPathA; 
        allScore{counter} = fScore; 
        counter = counter + 1;
    end 
    %if debugg
        %pause(1.0)
    %end 
    score{p} = fScore; 
    coordStore{p} = curPathA; 


end 

scoreArray = vertcat(score{:});
[~,indx] = min(scoreArray); 
% Return the best match 
bestMatch = coordStore{indx}; 
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
fig = figure(1), imshow(image, 'InitialMagnification', 700); hold on
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
