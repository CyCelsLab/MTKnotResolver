%% Written by Dhruv Khatri, IISER Pune. 
% Ref: Mathworks website Exploring shortest path 
% https://blogs.mathworks.com/steve/2011/12/02/exploring-shortest-paths-part-3/
function [sortSkele, referenceEnd] = sortSkel(curContour, endpoint, referenceEnd,prevContour)
%% SORTSKEL , sorts a skeleton in indices as they appear along the contour path.
% The shortest path starts at the seed point and ends at E1, this orients
% the contour label belonging to the seed point in ascending order. This
% function was implemented after reading some threads on MATHWORKS website,
% all credits to online resources. Funtion does not works with looped
% structure, therefore null_image needs to be preprocessed to avoid
% returing Inf values from the distance transform. 

sizeImage = size(curContour); 
skeleton = find(curContour);

% set the reference point 
if isempty(referenceEnd)
    referenceEnd = endpoint(1); 
else 
    [ye, xe] = ind2sub(sizeImage, endpoint); 
    [yr, xr] = ind2sub(prevContour, referenceEnd);
    dist_matrix = pdist2([ye, xe], [yr, xr]); 
    [~, I] = min(dist_matrix, [],1);
    referenceEnd = endpoint(I); 
end 

% sort the indices 
[R,C] = ind2sub(sizeImage, referenceEnd) ; 
D1  =  bwdistgeodesic(curContour,C,R,'quasi-Euclidean'); 
distanceValue = D1(skeleton); 
[~,I] = sort(distanceValue); 
sortSkele = skeleton(I); 

%referenceEnd = sortSkele(1); 
end 



