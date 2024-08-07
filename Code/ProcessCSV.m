function outMat = ProcessCSV(mtcsv, endEndDistance,scalingFactor,endDistancemicron,timeStep, figurePath)
    % Check if the scalingFactor is provided, otherwise set a default value
    if nargin < 2
        scalingFactor = 106/1000; % Default scaling factor value (you can change it as needed)
    end
    
    % Set the endDistance (This distance is from the end)
    endDistancemicron = ceil(endDistancemicron/scalingFactor);  % convert micron to px
    % Read the CSV file
    data = mtcsv;  % readmatrix(mtcsv);

    % Extract frame number, x-coordinates, and y-coordinates from the data
    frameNumber = data(:, 1);
    xCoordinates = data(:, 2:end-1);
    yCoordinates = data(:, end);

    % Calculate tip angles for each frame
    tipAngles = CalculateTipAngles(xCoordinates, yCoordinates, frameNumber,endDistancemicron);

    tipAnglesRaw = tipAngles; % Replace with your distance end  

    % Assign the tipAngles array to the output variable (outMat)
    outMat = tipAngles;
    try 
    tipAngles = smooth(tipAngles, 0.2, 'sgolay');
    catch
        tipAngles = tipAngles
    end 
    % Perform FFT analysis on the tip angles
    fftResult = fft(tipAngles);
    amplitude = abs(fftResult);
    Ts = timeStep;
    fs = 1/Ts; 

    frequencies = (0:length(tipAngles) - 1)*fs / length(tipAngles); % Normalized frequencies
     % Plot the data of tip angles
    g = figure(10);subplot(2,1,1)

    plot(unique(frameNumber).*timeStep, tipAngles, 'b-', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Tip Angle (degrees)');
    title('Tip Angle Variation');
    set(gca, 'FontSIze', 16)
    grid on;
    
    % Consider only the first half of the FFT result
    n = length(tipAngles);
    halfN = floor(n / 2) + 1;
    amplitudeDom = amplitude(2:halfN);
    halfFreq = frequencies(2:halfN); 
    % Find the most dominant frequency
    [~, maxIndex] = max(amplitudeDom); % Get the index of the maximum amplitude
    dominantFrequency = halfFreq(maxIndex); % Get the corresponding frequency



    % Plot the FFT analysis results
    figure(10);subplot(2,1,2);
    stem(frequencies, amplitude, 'r-', 'LineWidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    title('FFT Analysis of Tip Angle Variation');
    grid on;
    %% GENERATE SUMMARY
    inputData = data;        % Replace with your input file name
    outputFileName = 'statsSummary.txt'; 
    outputFilePath = fullfile(figurePath,outputFileName) ;% Replace with your output file name
    generateStatsSummary(inputData,tipAnglesRaw,tipAngles, dominantFrequency, outputFilePath,scalingFactor, timeStep, endDistancemicron, endEndDistance);
end

function tipAngles = CalculateTipAngles(xCoordinates, yCoordinates, frameNumber,d)

    tipAngles = []; % Initialize an empty array to store tip angles
    fnum = unique(frameNumber); 
    for frame = 1:length(fnum)
        indxFrame = find(frameNumber==fnum(frame));  
        x = xCoordinates(indxFrame, :);
        y = yCoordinates(indxFrame, :);

        % Calculate the vector connecting the tip and a point d distance away
        % from the tip using the provided scalingFactor
        try 
        xx = smooth(x, 0.2, 'sgolay');
        yy = smooth(y, 0.2, 'sgolay'); 
        catch
            xx = x;
            yy= y; 
        end 
        xxa =  xx(1:6); 
        yya = yy(1:6); 

        P = polyfit(xxa,yya,1); 
        slope1 = P(1);


        if d ==0 
            xe = xx(end);
            ye = yy(end); 
        
        else
            xe = xx(d);
            ye = yy(d); 
        end 
        slope2 = (ye-yya(1)) / (xe - xxa(1)); 

        tangle = atand((slope1-slope2)/(1+ slope1*slope2)); 

        % Calculate the angle between the two vectors using the dot product formula
        tipAngles = [tipAngles; tangle];
    end
end

% Function to generate statistics summary
function generateStatsSummary(inputData, tipAnglesRaw, tipAngles, dominantFrequency, outputFileName, scalingFactor, timeStep, distanceEnd, endEndDistance)
    
    % Calculate frame-wise length statistics
    allF = length(unique(inputData(:,1)));
    outputData  = zeros(allF, 5); % Frame , Length , Tip Angle, Tip angle smooth, endToEndDistance
    for l = 1:length(unique(inputData(:,1)))
        curFrameData = inputData(inputData(:,1)==l, :); 
        xc = curFrameData(:,2) * scalingFactor;
        yc = curFrameData(:,3) * scalingFactor; 
        lengthFilament = lengthFil(xc,yc);  
        tipAngleList = tipAnglesRaw(l);
        tipAngleListSmooth = tipAngles(l); 
        endToEndList = endEndDistance(l); 

        outputData(l,:) = [unique(curFrameData(:,1)), lengthFilament, ...
            tipAngleList, tipAngleListSmooth, endToEndList];
    end 

    % Calculate statistics
    frameLengths = outputData(:,2); 
    meanLength = mean(frameLengths);
    medianLength = median(frameLengths);
    stdLength = std(frameLengths);
    minLength = min(frameLengths);
    maxLength = max(frameLengths);
    
    % Write header information and metadata to the file
    fileID = fopen(outputFileName, 'w');
    
    fprintf(fileID, 'Statistics Summary\n');
    fprintf(fileID, '------------------\n');
%     fprintf(fileID, 'Input Filename: %s\n', inputFileName);
    fprintf(fileID, 'Scaling Factor: %.5f\n', scalingFactor);
    fprintf(fileID, 'Time Step: %.2f\n', timeStep);
    fprintf(fileID, 'Distance ignored from end (microns): %.2f\n', distanceEnd);
    fprintf(fileID, 'Dominant Frequency (Hz): %.5f\n', dominantFrequency);
%     fprintf(fileID, 'Segmentation Parameters: %s\n', segmentationParams);
    fprintf(fileID, '\n');
    
    % Write the table header
    fprintf(fileID, 'Mean Length, Median Length, Std Length, Min Length, Max Length\n');
    
    % Write the statistics
    fprintf(fileID, '%.2f, %.2f, %.2f, %.2f, %.2f\n', meanLength, medianLength, stdLength, minLength, maxLength);
    

    fprintf(fileID, 'Frame, Length (micron), Tip Angle, Tip Angle (smooth), EndToEndDistance (micron)\n');
    for i = 1:length(outputData(:,1))
            fprintf(fileID, '%d, %.3f, %.3f, %.3f, %.3f\n', ...
            outputData(i,1), outputData(i,2), outputData(i,3), outputData(i,4), outputData(i,5));
    end 
    fclose(fileID);
    
    % Display completion message
    fprintf('Statistics summary has been written to %s\n', outputFileName);
end
function getLength = lengthFil(xc,yc)
getLength = sum(sqrt(diff(xc).^2 + diff(yc).^2)); 
end 
