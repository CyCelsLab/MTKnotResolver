function tangentPlot(resolveCoordinates)
sizemat = cellfun(@size, {resolveCoordinates.Skeleton}, 'UniformOutput',false);
maxSize = max(cellfun(@max, sizemat));
maxRow = max(cat(1,resolveCoordinates.FrameNumber));
matDisTangent = zeros(maxRow, maxSize);
matDisTangent(matDisTangent == 0) = nan;
endEndDistance= []; 
for l = 1:length(resolveCoordinates)
    %%
    singleCurve = resolveCoordinates(l).Skeleton;
    [yc, xc] = ind2sub(resolveCoordinates(l).smallSize, singleCurve);
    frameN = resolveCoordinates(l).FrameNumber;
    Offset = resolveCoordinates(l).Offset;

%     figure(1), plot(smooth(xc, 0.2), smooth(yc, 0.2), 'b-'); hold off
%     
%     xc  = flip(xc);
%     yc = flip(yc); 
% 
%     figure(2), plot(xc,yc); 
%     

    %xci = interp1(1:length(xc), xc, 1:0.01:length(xc), 'spline');
    %yci = interp1(1:length(yc), yc, 1:0.01:length(yc), 'spline');
    %figure(1), plot(xci, yci, 'r--');

    %figure(1), plot(smooth(xci,0.05), smooth(yci,0.05), 'k--'); hold off
    %legend(["Raw", "Smooth"]);
    
    xc = xc(:) + Offset(2);
    %xc = flip(xc);
    yc = yc(:) + Offset(1);
    %yc = flip(yc);

    xc = smooth(xc, 0.2, "sgolay");
    yc = smooth(yc, 0.2, "sgolay"); 

   
    yc1 = yc - yc(1); 
    xc1 = xc - xc(1);

    yc1 = yc1 * (106/1000); 
    xc1 = xc1 * (106/1000); 


    endEndDistance= [endEndDistance, pdist([xc1(1), yc1(1);xc1(end), yc1(end)])]; 

    g = figure(1),subplot(3,1,1); patch([xc1' nan], [yc1' nan], [1:length(xc) nan], ...
    'edgecolor', 'interp', 'LineWidth', 1.5); hold on 
    set(gca, 'FontSize', 14) ; hold off

    z = zeros(size(xc));
    % (1) Frenet-Serret orthonormal vector tangent (t) = dr/ds

    dx = gradient(xc);
    dy = gradient(yc);
    dz = gradient(z);


    ddx = gradient(dx);
    ddy = gradient(dy);
    ddz = gradient(dz);

    dr = [dx dy dz];
    ddr = [ddx ddy ddz];

    T = dr./mag(dr,3);
    matDisTangent(frameN, 1:length(xc1)-1) = tangentangle(xc1,yc1);
    % (2) Derivative of tangent
    dTx= gradient(T(:,1));
    dTy = gradient(T(:,2));
    dTz = gradient(T(:,3));

    dT = [dTx dTy dTz];

    % Normal
    N = dT./mag(dT,3);
    matDisNormal(frameN, 1:length(N(:,1))) = tanangle(N);
    % Binormal
    B = cross(T,N);

    % Curvature
    k = mag(cross(dr,ddr),1)./ ((mag(dr,1)).^3);
    matDisCurvature(frameN, 1:length(k)) = k;
end



colorbar; hold off
% print(gcf, '-dpdf', fullfile(figurePath, "Contour.pdf"),'-r600' )
% Show plots of matrices
%f = figure(1),subplot(3,1,1),  imagesc([1:maxSize]*(106/1000),[0:maxRow]*10,matDisCurvature);
%set(gca, 'FontSize', 16)
%cb1=colorbar;
%print(gcf, '-dpdf', fullfile("Figures", sprintf("Overlay%d.pdf",d)),'-r600')


%figure(2), imagesc(matDisNormal); colorbar
%print(gcf, '-dpdf', fullfile(inputPath, "Normal.pdf"),'-r600')

matDisTangent = matDisTangent/180; 
imAlpha = ones(size(matDisTangent)); 
imAlpha(isnan(matDisTangent)) = 0; 
matDisTangent(isnan(matDisTangent)) = 0 ;

tempMat = (matDisNormal> 0 ); 
tempMat = tempMat(~all(tempMat == 0, 2), :);
sizeTemp = size(tempMat); 
spreadZero = (prod(size(tempMat)) - length(nonzeros(tempMat)))/sizeTemp(1); 

figure(1),subplot(3,1,2), imagesc([1:maxSize]*(106/1000),[0:maxRow]*10,matDisTangent,'AlphaData',imAlpha);
set(gca,'color',0*[1 1 1]);
set(gca, 'FontSize', 14)
xlim([0 (maxSize-spreadZero)*(106/1000)])
cb=colorbar;
set(cb, 'Ticks', [-1, -0.5, 0 , 0.5 , 1], 'TickLabels', { '-\pi ','-\pi /2', '0', '\pi /2', '\pi '}, 'FontSize',16)
%     print(gcf, '-dpdf', fullfile(figurePath, "Tangent.pdf"),'-r600')


figure(1),subplot(3,1,3), plot(cat(1,resolveCoordinates.FrameNumber)*10, endEndDistance, 'k-', ...
    'LineWidth',2.0); 
set(gca, 'FontSize', 14)
g.Position = [1229 115 467 1200]; 
%     print(gcf, '-dpdf', fullfile(figurePath, "endEnd.pdf"),'-r600'


end 


function degreeArray  = tanangle(vectorarray)
[r,~] = size(vectorarray);
degreeArray = zeros(r, 1);
for t =  1:r
    val = vectorarray(t,1:2);
    i = val(1);
    j = val(2);
    if (i >=0 && j >=0)
        degreeArray(t) = abs(atand(j/i));
    else
        if (j < 0 && i < 0)
            degreeArray(t) = abs(atand(j/i)) + 180 ;

        else
            if (j < 0 && i >= 0)
                degreeArray(t) = 360-abs(atand(j/i));
            else
                degreeArray(t) = 180-abs(atand(j/i));
            end
        end

    end
end
end

function degreeArray  =  tangentangle(xc,yc)
degreeArray = atan2d(diff(yc), diff(xc)); % .*180/pi; 
end

function N = mag(T,n)
% MAGNITUDE OF A VECTOR (Nx3)
%  M = mag(U)
N = sum(abs(T).^2,2).^(1/2);
d = find(N==0);
N(d) = eps*ones(size(d));
N = N(:,ones(n,1));
end
