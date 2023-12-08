%% Created by Dhruv Khatri , IISER Pune. 
function dataArray = AnalysisFile(resolveCoordinates, filePath)
%AnalysisFile implements plotting of tangent angle along the contour and
% ouputs the csv file of framenumber , x and y coordinates. The output
% therefore are open figure of 3 subfigures (1) Overlay of contour with
% gradient coded along the length of the contour. (2) Tangent angle
% kymograph along. (3) End-to end distance between the two ends for the
% complete time series. 
% INPUT: 
% ResolveCoordinates: The struct output from the KnotResolver.m file. 
% filePath: File path of where the input is located. Creates a DemoOutput
% folder and saves the pdf of the fiugure and csv file. 
% Please be sure to correct the scaling factor in the script. 
% Scaling factor used on line 36,37; 


    % Get the size of each skeleton and find the maximum size
    sizemat = cellfun(@size, {resolveCoordinates.Skeleton}, 'UniformOutput', false);
    maxSize = max(cellfun(@max, sizemat));

    % Find the maximum frame number
    maxRow = max(cat(1, resolveCoordinates.FrameNumber));

    % Initialize a matrix to store tangent angles
    matDisTangent = zeros(maxRow, maxSize);
    matDisTangent(matDisTangent == 0) = nan;

    % Initialize an empty array to store end-to-end distances
    endEndDistance = [];

    % Loop through each entry in resolveCoordinates
    for l = 1:length(resolveCoordinates)
        % Extract coordinates of a single curve
        singleCurve = resolveCoordinates(l).Skeleton;
        [yc, xc] = ind2sub(resolveCoordinates(l).smallSize, singleCurve);
        frameN = resolveCoordinates(l).FrameNumber;
        Offset = resolveCoordinates(l).Offset;
        xc = xc(:) + Offset(2);
        yc = yc(:) + Offset(1);
        
        % Smooth the x and y coordinates using a Savitzky-Golay filter
        try 
        xc = smooth(xc, 0.2, "sgolay");
        yc = smooth(yc, 0.2, "sgolay");
        catch
            continue 
        end 
        % Normalize coordinates to a specific scale
        yc1 = yc - yc(1);
        xc1 = xc - xc(1);
        yc1 = yc1 * (106 / 1000);
        xc1 = xc1 * (106 / 1000);

        % Calculate and store end-to-end distances
        endEndDistance = [endEndDistance, pdist([xc1(1), yc1(1); xc1(end), yc1(end)])];

        % Plot the curve and set axis properties
        g = figure(1), subplot(3, 1, 1);
        patch([xc1' nan], [yc1' nan], [1:length(xc) nan], 'edgecolor', 'interp', 'LineWidth', 1.5);
        hold on;
        set(gca, 'FontSize', 14);
        hold off;

        % Calculate the Frenet-Serret orthonormal vector tangent (t) = dr/ds
        z = zeros(size(xc));
        dx = gradient(xc);
        dy = gradient(yc);
        dz = gradient(z);

        ddx = gradient(dx);
        ddy = gradient(dy);
        ddz = gradient(dz);
        dr = [dx dy dz];
        ddr = [ddx ddy ddz];
        T = dr ./ mag(dr, 3);

        % Store tangent angles in the matDisTangent matrix
        matDisTangent(frameN, 1:length(xc1) - 1) = tangentangle(xc1, yc1);

        % Calculate the derivative of the tangent (dT/ds)
        dTx = gradient(T(:, 1));
        dTy = gradient(T(:, 2));
        dTz = gradient(T(:, 3));
        dT = [dTx dTy dTz];

        % Calculate and store normal angles in the matDisNormal matrix
        % (Not used currently, may be incorrect) 
        N = dT ./ mag(dT, 3);
        matDisNormal(frameN, 1:length(N(:, 1))) = tanangle(N);

        % Calculate the binormal (B) and curvature (k)
        B = cross(T, N);
        k = mag(cross(dr, ddr), 1) ./ ((mag(dr, 1)).^3);
        matDisCurvature(frameN, 1:length(k)) = k;
    end

    % Create a colorbar and plot the tangent angles
    colorbar;
    xlabel('X (\mu m)');
    ylabel('Y (\mu m)'); 
    hold off
    matDisTangent = matDisTangent * (pi/ 180); % convert degree to rad/ might be 
    % some errors in plotting (adjust the colorbar accordingly)
    imAlpha = ones(size(matDisTangent));
    imAlpha(isnan(matDisTangent)) = 0;
    matDisTangent(isnan(matDisTangent)) = 0;

    % Calculate and plot spreadZero
    tempMat = (matDisNormal > 0);
    tempMat = tempMat(~all(tempMat == 0, 2), :);
    sizeTemp = size(tempMat);
    spreadZero = (prod(size(tempMat)) - length(nonzeros(tempMat))) / sizeTemp(1);

    figure(1), subplot(3, 1, 2);
    imagesc([1:maxSize] * (106 / 1000), [0:maxRow] * 10, matDisTangent, 'AlphaData', imAlpha);
    xlabel('MT (\mu m)')
    ylabel('Time (s)')
    set(gca, 'color', 0 * [1 1 1]);
    set(gca, 'FontSize', 14);
    xlim([0.1 (maxSize - spreadZero-1) * (106 / 1000)]);
    ylim([0 maxRow*10]);

    cb = colorbar;
    %set(cb, 'Ticks', [-1, -0.5, 0, 0.5, 1], 'TickLabels', {'-\pi ', '-\pi /2', '0', '\pi /2', '\pi '}, 'FontSize', 16);
    
    % Plot end-to-end distances
    figure(1), subplot(3, 1, 3);
    plot(cat(1, resolveCoordinates.FrameNumber) * 10, endEndDistance, 'k-', 'LineWidth', 2.0);
    ylabel('d_{e} (\mu m)') 
    xlabel('Time (s)')
    set(gca, 'FontSize', 14);
    g.Position = [1229 115 467 1200];

    
    %% Save the CSV file with tangent angles
    sizemat = cellfun(@size, {resolveCoordinates.Skeleton}, 'UniformOutput', false);
    maxSize = max(cellfun(@max, sizemat));
    dataL = sum(cellfun(@max, sizemat));
    maxRow = max(cat(1, resolveCoordinates.FrameNumber));
    dataArray = zeros(dataL, 3);
    start = 1;

    % Prepare data and save to a CSV file
    for l = 1:length(resolveCoordinates)
        singleCurve = resolveCoordinates(l).Skeleton;
        [yc, xc] = ind2sub(resolveCoordinates(l).smallSize, singleCurve);
        frameN = resolveCoordinates(l).FrameNumber;
        Offset = resolveCoordinates(l).Offset;
        yc = yc + Offset(1);
        xc = xc + Offset(2);
        dataArray(start:start + length(xc) - 1, 1:3) = [repmat(frameN, length(xc), 1), xc, yc];
        start = start + length(xc);
    end

    [figurePath, ~, ~] = fileparts(filePath);
    writematrix(dataArray, fullfile(figurePath, 'MTdata.csv'));
end

% Helper function to calculate tangent angles
function degreeArray = tanangle(vectorarray)
%DEGREEARRAY plots the tangent angle between -2\pi and 2\pi
    [r, ~] = size(vectorarray);
    degreeArray = zeros(r, 1);
    for t = 1:r
        val = vectorarray(t, 1:2);
        i = val(1);
        j = val(2);
        if (i >= 0 && j >= 0)
            degreeArray(t) = abs(atand(j / i));
        else
            if (j < 0 && i < 0)
                degreeArray(t) = abs(atand(j / i)) + 180;
            else
                if (j < 0 && i >= 0)
                    degreeArray(t) = 360 - abs(atand(j / i));
                else
                    degreeArray(t) = 180 - abs(atand(j / i));
                end
            end
        end
    end
end

% Helper function to calculate tangent angles
function degreeArray = tangentangle(xc, yc)
%TANGENTANGLE standaed tan angle calculation using in built atan2d
    degreeArray = atan2d(diff(yc), diff(xc));
end

% Helper function to calculate magnitude of a vector
function N = mag(T, n)
%MAG calculates the magnitude of a vector. 
    N = sum(abs(T).^2, 2).^(1 / 2);
    d = find(N == 0);
    N(d) = eps * ones(size(d));
    N = N(:, ones(n, 1));
end

