function fun_logProgress(i, nPoints)

    increment = 2.5; % Percentage increment
    keyPoints = increment:increment:100;
    
    currentStatus = round(i / nPoints * 100, 10);
    previousStatus = round((i - 1) / nPoints * 100, 10);
    
    if i == 1
        fprintf('Simulation started.\n');
    elseif i == nPoints
        fprintf('Simulation ended.\n');
    else
        crossed = (previousStatus < keyPoints) & (currentStatus >= keyPoints);
        if any(crossed)
            fprintf('Simulation progress: %.1f %%\n', keyPoints(find(crossed, 1)));
        end
    end

end
