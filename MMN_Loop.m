% Define the parameters for standard stimulus
standardParams = struct('ToneAmp', 0.2, 'ToneFreq', 1000, 'ToneDur', 250, 'ModAmp', 1, 'ModFreq', 20, 'ID_SweepTime', 150, 'ID_F1', 1000, 'ID_F2', 4000, 'StimType', 2);

% Define the parameters for deviant stimulus
deviantParams = struct('ToneAmp', 0.2, 'ToneFreq', 2000, 'ToneDur', 250, 'ModAmp', 1, 'ModFreq', 25, 'ID_SweepTime', 250, 'ID_F1', 1000, 'ID_F2', 4000, 'StimType', 2);

% Define the probability of a deviant stimulus
deviantProbability = 0.1;

% Initialize the previous stimulus type
previousStimulusType = 'standard';

% Define interstimulus interval (in milliseconds)
interstimulusInterval = 600; % 

% Specify the number of trials
numTrials = 100;

% Initialize the previous stimulus type
previousStimulusType = 'standard';

% Generate a list of trial types with the desired proportion of deviants
deviantIndices = randperm(numTrials, round(deviantProbability * numTrials));

% Ensure there are no consecutive deviants by shuffling again if necessary
sorted_dev = sort(deviantIndices);
while any(diff(sorted_dev) == 1)
    deviantIndices = randperm(numTrials, round(deviantProbability * numTrials));
    sorted_dev = sort(deviantIndices);
end

% Create a table to store parameters
paramTable = table('Size', [numTrials, 12], 'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'Trial', 'ToneAmp', 'ToneFreq', 'ToneDur', 'ModAmp', 'ModFreq', 'ID_SweepTime', 'ID_F1', 'ID_F2', 'StimType', 'DeviantProbability', 'InterstimulusInterval'});


% Add constant values to the table
paramTable.Trial = (1:numTrials)';
paramTable.DeviantProbability(:) = deviantProbability;

try
    for trial = 1:numTrials
        % Determine whether the current trial should be a standard or deviant
        if any(trial == deviantIndices)
            currentParams = deviantParams;
            previousStimulusType = 'deviant';
        else
            currentParams = standardParams;
            previousStimulusType = 'standard';
        end
        
        % Set the parameters for the current trial
        RP.SetTagVal('ToneAmp', currentParams.ToneAmp);
        RP.SetTagVal('ToneFreq', currentParams.ToneFreq);
        RP.SetTagVal('ToneDur', currentParams.ToneDur);
        RP.SetTagVal('ModAmp', currentParams.ModAmp);
        RP.SetTagVal('ModFreq', currentParams.ModFreq);
        RP.SetTagVal('ID_SweepTime', currentParams.ID_SweepTime);
        RP.SetTagVal('ID_F1', currentParams.ID_F1);
        RP.SetTagVal('ID_F2', currentParams.ID_F2);
        RP.SetTagVal('StimType', currentParams.StimType);
        
        % Trigger the stimulus presentation
        RP.SoftTrg(1);
        
        % Pause for the interstimulus interval
        pause(interstimulusInterval / 1000); % Convert to seconds
        
        % Update the table with current trial parameters
        paramTable(trial, 2:10) = struct2table(currentParams);
    end
    
    % Generate a date and time string
    dateTimeString = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    
    % Save the entire table to a CSV file with the date and time in the name
    writetable(paramTable, ['trial_parameters_' dateTimeString '.csv']);
    
catch exception
    % Display the error message
    disp('An error occurred: something stopped the loop before it could finish. SAVING DATA!');
    disp(exception.message);
    
    % Save the current data to a CSV file before exiting
    dateTimeString = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    writetable(paramTable, ['error_data_' dateTimeString '.csv']);
    
    % Re-throw the exception to exit the script
    rethrow(exception);
end

