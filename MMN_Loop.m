% Define the parameters for standard stimulus
standardParams = struct('ToneAmp', 0.15, 'ToneFreq', 1000, 'ToneDur', 300, 'ModAmp', 1, 'ModFreq', 20, 'ID_SweepTime', 200, 'ID_F1', 1000, 'ID_F2', 4000, 'StimType', 2);

% Define the parameters for deviant stimulus
deviantParams = struct('ToneAmp', 0.15, 'ToneFreq', 1200, 'ToneDur', 300, 'ModAmp', 1, 'ModFreq', 40, 'ID_SweepTime', 300, 'ID_F1', 1000, 'ID_F2', 4000, 'StimType', 2);

% Define the probability of a deviant stimulus
deviantProbability = 0.1;

% Initialize the previous stimulus type
previousStimulusType = 'standard';

% Define interstimulus interval (in milliseconds)
interstimulusInterval = 600; % 1 second

% Specify the number of trials
numTrials = 100;

for trial = 1:numTrials
    
    % Determine whether the current trial should be a standard or deviant
    if rand() < deviantProbability
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

end
