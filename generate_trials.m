function generate_trials(standardParams, deviantParams1, deviantProbability1, interstimulusInterval, numTrials, paramsDir)
    % Calculate number of trials for deviants 
    numDeviants = ceil(numTrials * deviantProbability1);

    % Ensure deviants don't exceed 10% of total trials
    if numDeviants > 0.1 * numTrials
        numDeviants = floor(0.1 * numTrials);
    end

    % Initialize cell array to store trial types
    trialTypes = cell(1, numTrials);


    minimumSpacing = 4:1:10;

    % Initialize the first trial index
    trialIdx = 1;
    deviantIndices = [];
    for i = 1:numDeviants
        % Add a random spacing from the minimumSpacing range
        spacing = randi([min(minimumSpacing), max(minimumSpacing)], 1, 1); 
        
        % Compute the next trial index by adding the spacing
        trialIdx = trialIdx + spacing;

        % Assign the deviant to the calculated trial index
        trialTypes{trialIdx} = 'D';
        deviantIndices(i) = trialIdx;
    end
    for i = 1:numTrials
        if isempty(trialTypes{i})
            trialTypes{i} = 'S';
        end
    end

    % Open the text files for each parameter
    toneAmpFile = fopen(fullfile(paramsDir, 'ToneAmp.txt'), 'w');
    toneFreqFile = fopen(fullfile(paramsDir, 'ToneFreq.txt'), 'w');
    toneDurFile = fopen(fullfile(paramsDir, 'ToneDur.txt'), 'w');
    modAmpFile = fopen(fullfile(paramsDir, 'ModDepth.txt'), 'w');
    modFreqFile = fopen(fullfile(paramsDir, 'ModFreq.txt'), 'w');
    sweepTimeFile = fopen(fullfile(paramsDir, 'FMSweepTime.txt'), 'w');
    f1File = fopen(fullfile(paramsDir, 'FM1.txt'), 'w');
    f2File = fopen(fullfile(paramsDir, 'FM2.txt'), 'w');
    stimTypeFile = fopen(fullfile(paramsDir, 'StimType.txt'), 'w');
    isiFile = fopen(fullfile(paramsDir, 'ISI.txt'), 'w');
    deviantFile = fopen(fullfile(paramsDir, 'Deviant.txt'), 'w'); % New file for deviant info

    try
        for trial = 1:numTrials
            % Determine the current trial type
            currentTrialType = trialTypes{trial};

            % Set the parameters for the current trial
            switch currentTrialType
                case 'S'
                    currentParams = standardParams;
                    deviantFlag = 0; % Standard trial
                case 'D'
                    currentParams = deviantParams1;
                    deviantFlag = 1; % Deviant trial
            end

            % Write the parameters to their respective text files
            fprintf(toneAmpFile, '%f\n', currentParams.ToneAmp);
            fprintf(toneFreqFile, '%f\n', currentParams.ToneFreq);
            fprintf(toneDurFile, '%f\n', currentParams.ToneDur);
            fprintf(modAmpFile, '%f\n', currentParams.ModAmp);
            fprintf(modFreqFile, '%f\n', currentParams.ModFreq);
            fprintf(sweepTimeFile, '%f\n', currentParams.ID_SweepTime);
            fprintf(f1File, '%f\n', currentParams.ID_F1);
            fprintf(f2File, '%f\n', currentParams.ID_F2);
            fprintf(stimTypeFile, '%f\n', currentParams.StimType);
            fprintf(isiFile, '%f\n', interstimulusInterval);

            % Write the deviant flag to the deviant file
            fprintf(deviantFile, '%d\n', deviantFlag);
        end

        % Print indices for deviant at the end
        fprintf('Indices of deviants: %s\n', num2str(find(strcmp(trialTypes, 'D'))));

    catch exception
        % Display the error message
        disp('An error occurred: something stopped the loop before it could finish. Closing files and saving data.');
        disp(exception.message);
    end

    % Close the text files
    fclose(toneAmpFile);
    fclose(toneFreqFile);
    fclose(toneDurFile);
    fclose(modAmpFile);
    fclose(modFreqFile);
    fclose(sweepTimeFile);
    fclose(f1File);
    fclose(f2File);
    fclose(stimTypeFile);
    fclose(isiFile);
    fclose(deviantFile); % Close the new deviant file
end