function outMat = ProcessCSV(mtcsv, scalingFactor,endDistancemicron,timeStep)
    close all;
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
    g = figure(1); subplot(1,2,1); 

    plot(unique(frameNumber).*timeStep, tipAngles, 'b-', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Tip Angle (degrees)');
    title('Tip Angle Variation');
    set(gca, 'FontSIze', 16)
    grid on;


    % Plot the FFT analysis results
    figure(1); subplot(1,2,2); 
    stem(frequencies, amplitude, 'r-', 'LineWidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    title('FFT Analysis of Tip Angle Variation');
    grid on;

    g.Position = [995 909 1016 413];
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
