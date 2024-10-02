function generate_trials(standardParams, deviantParams1, deviantParams2, deviantProbability1, deviantProbability2,interstimulusInterval, numTrials, paramsDir)
    % Check if deviantParams2 exists to determine if we need one or two deviants
    hasSecondDeviant = exist('deviantParams2', 'var') && ~isempty(deviantParams2);

    % Calculate the number of trials for each deviant type based on specified probabilities
    numDeviants1 = floor(deviantProbability1 * numTrials);

    if hasSecondDeviant
        numDeviants2 = floor(deviantProbability2 * numTrials);
    else
        numDeviants2 = 0;  % No second deviant
    end

    % Initialize cell array to store trial types
    trialTypes = cell(1, numTrials);

    % Randomly assign deviant positions ensuring at least 4 standards between deviants
    deviant1Indices = sort(randperm(numTrials, numDeviants1));
    if hasSecondDeviant
        deviant2Indices = sort(randperm(numTrials, numDeviants2));
    end

    % Ensure at least 4 standards between consecutive deviants
    lastD1Index = 0;
    lastD2Index = 0;

    % Initialize lists to hold deviant indices
    dev1I = [];
    dev2I = [];

% Available positions for deviants
    availablePositions = 1:numTrials;
    
    % Ensure at least 4 standards between consecutive deviants
    minSpacing = 4;

    % Assign Type 1 Deviants (D1) ensuring spacing
    for i = 1:numDeviants1
        % Check if there are still enough positions left
        potentialPositions = availablePositions(availablePositions > minSpacing & availablePositions <= numTrials - minSpacing);
        if isempty(potentialPositions)
            warning('Not enough valid positions left for deviant 1 placement.');
            break;
        end
        
        % Randomly pick a position from the valid ones
        idx = randperm(length(potentialPositions), 1);
        trialIdx = potentialPositions(idx);
        
        % Remove any positions that would violate the spacing constraints
        availablePositions = setdiff(availablePositions, trialIdx - minSpacing : trialIdx + minSpacing);
        
        trialTypes{trialIdx} = 'D1';
        dev1I(end+1,:) = trialIdx;
    end

    % Assign Type 2 Deviants (D2) ensuring spacing
    if hasSecondDeviant
        for i = 1:numDeviants2
            % Check if there are still enough positions left
            potentialPositions = availablePositions(availablePositions > minSpacing & availablePositions <= numTrials - minSpacing);
            if isempty(potentialPositions)
                warning('Not enough valid positions left for deviant 2 placement.');
                break;
            end

            % Randomly pick a position from the valid ones
            idx = randperm(length(potentialPositions), 1);
            trialIdx = potentialPositions(idx);
            
            % Remove any positions that would violate the spacing constraints
            availablePositions = setdiff(availablePositions, trialIdx - minSpacing : trialIdx + minSpacing);
            
            trialTypes{trialIdx} = 'D2';
            dev2I(end+1,:) = trialIdx;
        end
    end
    % Fill the remaining trials with standard stimuli
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
                case 'D1'
                    currentParams = deviantParams1;
                    deviantFlag = 1; % Deviant type 1
                case 'D2'
                    currentParams = deviantParams2;
                    deviantFlag = 2; % Deviant type 2
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

        % Print indices for deviant types
        disp('Indices of deviant type 1:');
        dev1I
        if hasSecondDeviant
            disp('Indices of deviant type 2:');
            dev2I
        end

    catch exception
        % Handle errors and save data safely
        disp('An error occurred: stopping the loop early. Saving data.');
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
    fclose(deviantFile);
end
