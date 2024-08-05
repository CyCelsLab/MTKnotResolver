%% Written by Dhruv Khatri, IISER Pune. 
% Contributions from Prince Shah, Prof. Chaitanya Athale, IISER Pune,
% Cycels Lab. 


%% SORTBRANCH
function [sortedSkeleton, referenceEnd, allpossibilites,allScore] = sortBranch_v3(curContour, E1, previousSkeleton, ...
    referenceEnd,refSize,refOffset,Offset,manualEdit,frameInvertTrue)
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
%% Create all possible paths between start and endpoint 
roughArray = roughEstimate_v2(curContour, label_image,previousSkeleton, ...
    refSize, refOffset,sizeImage, Offset,startPoint, endPoint, G, ...
    Trueloop, filamentTable,frameInvertTrue); % get a rough estimate 


possiblePaths = getPaths(G, startPoint, endPoint, group, roughArray, maxLabelImg);
%disp(possiblePaths)
% Create coordinates indices of all possibilities and return the best match
% if we have a manual edit, get the correct path input manually; 
if manualEdit
    set(0,'units','pixels')
    dim_desk = get(0,'ScreenSize');
    pts = readPoints(label2rgb(label_image), dim_desk); 
    pts = round(pts'); 
    pts_2_ind = sub2ind([512, 512], pts(:,2), pts(:,1)); 
    resizeLabel = imresize(label_image, [512 512], "nearest"); 
    labels_correct = resizeLabel(pts_2_ind); 
    %possiblePaths = {}; 
    possiblePaths = labels_correct';
end 
%% Arrange coordinates
% Arrange the coordinates for all paths according to the filament table
% How to avoid redundant calculations? We will address this in Future :)
szP = size(possiblePaths); 
arrangement1 = cell(szP(1),1); 
parfor l = 1:szP(1)
    curPath = nonzeros(possiblePaths(l,:)); 
    %curPath = nonzeros(possiblePaths{l}); 
    contourIndice  = sortCoordinates(curPath,label_image, ...
    Trueloop,filamentTable); 
    arrangement1{l} = contourIndice; 
end 

%% Format coordinates and compute score 
% Currently we will only compute two combinates, i.e, the solution performs
% poorly when more than one loop is present, [Future Work]
[sortedSkeleton, allpossibilites,allScore] = selectBest(arrangement1, previousSkeleton, ...
    refSize,refOffset,Offset,sizeImage, frameInvertTrue,"all");
%{
figure(1), imshow(double(repmat(curContour,1,1,3)), 'InitialMagnification', 200); hold on 
[yy,xx]= ind2sub(size(label_image), arrangement1{9});
%plot(xx,yy, 'r--'); hold off  
%pause(0.5)

patch([xx' nan], [yy' nan], [(1:length(yy))*0.5 nan], ...
'edgecolor', 'interp', 'LineWidth',4.0); hold off 
axis on
set(gca, 'FontSize', 16);
%pause(0.5)

%[sortedSkeleton, allpossibilites,allScore] = selectBest(label_image, possiblePaths, ...
%previousSkeleton,referenceStart,refSize,refOffset,Offset,Trueloop); 
%}
end

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

%% GETPATHS

function possiblePaths = getPaths_v2(G,shortPath)
% Function getPaths finds all possible paths between a startpoint and
% endpoint for a given graph component
% Possible bugs: NAN
%% Use Rough Array to create path combinations. 
for s = 1:length(shortPath)
end 

end 


function possiblePaths = getPaths(G,startPoint, endPoint,group, roughArray, maxLabelImage)
% Function getPaths finds all possible paths between a startpoint and
% endpoint for a given graph component
% Possible bugs: NAN
%% Use Rough Array to create path combinations. 
possiblePaths = {};  
if ~isempty(roughArray)
possiblePaths{1} = roughArray; 
countP = 2;
else
    countP= 1; 
end 
% We need to improve this section. 
uniq_values = unique(group);
uniq_values = setdiff(uniq_values, [startPoint, endPoint]); 
i = countP ; 
for p  = 1:length(uniq_values)
    path1 = allpaths(G,startPoint, uniq_values(p),"MaxNumPaths",10);
    if ~isempty(endPoint)
        path2 = allpaths(G,uniq_values(p), endPoint,"MaxNumPaths",10);

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
        %continue 
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

function roughArray = roughEstimate_v2(curContour, label_image,previousSkeleton, ...
    refSize, refOffset,sizeImage, Offset, startPoint,endPoint, G, ...
    Trueloop, filamentTable,frameInvertTrue)
orderPx = find(curContour); 

[y1, x1] = ind2sub(sizeImage, orderPx); 
y1 = y1 + Offset(1); 
x1 = x1 + Offset(2); 

[y2,x2] = ind2sub(refSize, previousSkeleton);
y2 = y2 + refOffset(1); 
x2 = x2 + refOffset(2); 


[~ , IndxMat]= pdist2([x1,y1], [x2,y2], "euclidean", "Smallest", 3); 

% figure(1), plot(x1(IndxMat(1,:)), y1(IndxMat(1,:)), 'r--', "LineWidth",2); hold on
% figure(1), plot(x1, y1, 'r--', "LineWidth",2); hold on
% % 
%  figure(1), plot(x2,y2, 'b--'); hold off

roughArray = label_image(orderPx(IndxMat(1,:))); 

remI =zeros(length(roughArray),1);  
countArray = []; 
startCount = 1; 
for t  = 1:length(roughArray)-1
    if roughArray(t) == roughArray(t+1)
        remI(t) = t;  
        startCount  = startCount + 1; 
    else
        countArray = [countArray; startCount]; 
        startCount = 1;
    end 
end 

countArray = [countArray; startCount]; 

roughArray(nonzeros(remI)) =[]; 
roughArray = [roughArray, countArray]; 

threshmajor = max(roughArray(:,2)); 
if threshmajor > 8  
   threshmajor = 8; 
end
majorLabels = roughArray(roughArray(:,2) >= threshmajor);

if majorLabels(1) ~= startPoint
    majorLabels = [startPoint; majorLabels]; 
end 

if majorLabels(end) ~= endPoint
    majorLabels = [majorLabels; endPoint];
end

shortPath = {[]};
delonward = [];
for  s = 1:length(majorLabels)-1
    allPoss = allpaths(G,majorLabels(s), majorLabels(s+1));
    if ~isempty(allPoss)
        indxDel = [] ;
        if length(allPoss) >= 4
            % Check is any allPOss goes through another mjor label
            for t = 1:length(allPoss)
                currPoss = allPoss{t};
                remMajor = setdiff(majorLabels,[majorLabels(s+1),majorLabels(s)]);
                nParts = sum(ismember(setdiff(remMajor,[majorLabels(s), majorLabels(s+1)]), currPoss));
                if (nParts >= 2) || (length(currPoss) > 10)
                    indxDel = [indxDel, t];
                end
            end
        end
        allPoss(indxDel) =[] ;

        shortPath{s} = allPoss;% shortestpath(G,majorLabels(s), majorLabels(s+1));
        if isempty(allPoss) && isempty(delonward)
            delonward = s;
            disp('Template matching failed on current frame');
        end
    else
        if isempty(delonward)
            delonward = s;
            disp('Template matching failed on current frame');
        end
    end
end
shortPath(delonward:end)= [];

% Process shortPath and only select the best matching option
if length(shortPath) > 1
    firstC = [] ;
    for ss = 1:3:length(shortPath)
        first = shortPath{ss};
        oii = length(first);
        try
            second = shortPath{ss+1};
            ojj = length(second);
        catch
            second = {};
            ojj= 0;
        end
        try
            third = shortPath{ss+2};
            okk = length(third);

        catch
            third = {} ;
            okk= 0;

        end
        tempPossi = {};
        counter1 = 1;
        ii = oii;
        jj = ojj;
        kk = okk;
        while ii > 0
            part1 = first{ii};
            part1 = [firstC, part1];
            if jj == 0
                tempPossi{counter1} = part1;
                counter1 = counter1 + 1;
            else
                while jj > 0
                    part2 = [part1, second{jj}];
                    if kk == 0
                        tempPossi{counter1} = part2;
                        counter1 = counter1 + 1;
                    else
                        while kk > 0
                            part3 = [part2, third{kk}];
                            tempPossi{counter1} = part3;
                            counter1 = counter1 + 1;
                            kk = kk - 1;
                        end
                    end
                    jj = jj - 1;
                    kk = okk;
                end
            end
            ii =ii - 1;
            jj = ojj;
        end
        % Select the best combination from tempPossi and propagate it further
        scoreArray = zeros([length(tempPossi),1])';
        possiArray = {};
        parfor t  = 1:length(tempPossi)
            currPoss = tempPossi{t};
            currPoss = removeDuplicates(currPoss);
            currIndices = sortCoordinates(currPoss',label_image, ...
                Trueloop, filamentTable);
            possiArray{t} = currIndices;
        end
        % select and return top two
        [~, ~,allScore] = selectBest(possiArray, previousSkeleton, ...
            refSize,refOffset,Offset,sizeImage, frameInvertTrue,"rough");
        [~, II] = min(allScore);
        firstC = tempPossi{II(1)};
    end

else
    firstC = []; 
    %if isempty(firstC)
        %firstC =  roughArray(:,1); 
    %end
end
roughArray= removeDuplicates(firstC);

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
function [sortedSkeleton, allpossibilites,allScore] = selectBest(arrangement1, previousSkeleton, ...
    refSize,refOffset,Offset,sizeImage, frameInvertTrue,typeCase)
scoreArray = zeros(length(arrangement1),1) ;
orientation = zeros(length(arrangement1),1);
for k = 1:length(arrangement1)
    lengthPrev = length(previousSkeleton);
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
        switch typeCase
            case "rough"
                if length(seq1) < length(previousSkeleton)
                    previousSk = previousSkeleton(1:length(seq1)); 
                else
                    previousSk = previousSkeleton; 
                end 
            
            case "all"
                previousSk = previousSkeleton; 
        end 
        score1 = computeScore(previousSk,refSize,refOffset, seq1,sizeImage,Offset);
        score2 = computeScore(previousSk, refSize, refOffset, seq2,sizeImage, Offset);
        %if (isempty(endPoint)) && (length(seq1) > lengthPrev)
        %    seq2 = seq2(1:lengthPrev);
        %    seq1 = seq1(1:lengthPrev);
        %end

        % Compute DTW Score
        if score2 > score1
            % should we overwrite the arrangement here?
            scoreArray(k) = score1;
            orientation(k) = 1;
            arrangement1{k} = seq1;
            if frameInvertTrue
                scoreArray(k) = score2;
                orientation(k) = 2;
                arrangement1{k} = seq2;
            end
        else
            scoreArray(k) = score2;
            orientation(k) = 2;
            arrangement1{k} = seq2;
            if frameInvertTrue
                scoreArray(k) = score1;
                orientation(k) = 1;
                arrangement1{k} = seq1;
            end
        end
    else
        seq1 = skelCoord(:,1);
        switch typeCase
            case "rough"
                if length(seq1) < length(previousSkeleton)
                    previousSk = previousSkeleton(1:length(seq1)); 
                else
                    previousSk = previousSkeleton; 
                end 
            
            case "all"
                previousSk = previousSkeleton; 
        end 

        arrangement1{k} = seq1;
        % Compute DTW score

        score1 = computeScore(previousSk, refSize, refOffset, seq1,sizeImage, Offset);
        scoreArray(k) =score1;
        orientation(k) = 1;

        %if (isempty(endPoint)) && (length(seq1) > lengthPrev)
        %    seq1 = seq1(1:lengthPrev);
        %end

    end
end
%% Select best  [sortedSkeleton, referenceStart, allpossibilites,allScore]
[~, Indm] = min(scoreArray);
sortedSkeleton = arrangement1{Indm};
allScore = scoreArray;
allpossibilites = arrangement1;
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
function matcell = cell2mat(cellArray)
% Remove duplicates 
maxPath = max(cellfun(@length, cellArray)); 
resizePossible = cellfun(@(x) [x zeros(1, maxPath-numel(x))], cellArray,'UniformOutput',false);
resizePossible = vertcat(resizePossible{:});
[~,jj] = unique(resizePossible, 'rows'); 
matcell = resizePossible(jj,:); 
end
function remDup = removeDuplicates(roughLength)
repeatsArray = find(diff(roughLength,1) == 0); 
roughLength(repeatsArray) = []; 

delList = []; 
for s = 2:length(roughLength)-3
    curr = roughLength(s); 
    third = roughLength(s+2); 
    if curr == third
        prev = roughLength(s-1);
        fourth = roughLength(s+3); 
        if prev ~= fourth    
            delList = [delList, s,s+1]; 
        end 
    end 
end   
roughLength(delList) = [];
remDup = roughLength; 
end