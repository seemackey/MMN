tic
filename='xxx';
% set up actx server/control 
handles.RP = actxcontrol('RPco.x');
RP	= handles.RP;

% connect to the device and halt any ongoing processes
RP.ConnectRX6('USB',1);
RP.Halt;
RP.ClearCOF;
    
% load our rcx file and run it
%RP.LoadCOF('C:\TDT\RPvdsEx\Examples\TDT_test_soft_trig.rcx');
RP.LoadCOF('C:\MMN-main\MMN_withTRig.rcx');
RP.Run; %

% Define the parameters for standard stimulus. %Stimtype #0 is umodulated tone,#1 is AM and #2 is FM
standardParams = struct('ToneAmp', .070, 'ToneFreq', 2000, 'ToneDur', 100, 'ModAmp', 1, 'ModFreq', 20, 'ID_SweepTime', 100, 'ID_F1', 2000, 'ID_F2', 12000, 'StimType', 0); 

% Define the parameters for deviant stimulus. %Stimtype #0 is umodulated tone,#1 is AM and #2 is FM
deviantParams1 = struct('ToneAmp', 0, 'ToneFreq', 2000, 'ToneDur', 100, 'ModAmp', 1, 'ModFreq', 20, 'ID_SweepTime', 100, 'ID_F1', 2000, 'ID_F2', 12000, 'StimType', 0);

% Define the probability of a deviant stimulus
deviantProbability1 = 0.1; %should be 0.1

% Define interstimulus interval (in milliseconds)
interstimulusInterval = 624; % 

% Specify the number of trials
numTrials = 1000;
generate_stimuli(RP, standardParams, deviantParams1, deviantProbability1, interstimulusInterval, numTrials,filename);
toc

function generate_stimuli(RP, standardParams, deviantParams1, deviantProbability1, interstimulusInterval, numTrials,filename)
    % Calculate number of trials for deviants (silent tones)
    numDeviants = ceil(numTrials * deviantProbability1);

    % Ensure deviants don't exceed 10% of total trials
    if numDeviants > 0.1 * numTrials
        numDeviants = floor(0.1 * numTrials);
    end

    % Initialize cell array to store trial types
    trialTypes = cell(1, numTrials);

    % Randomly assign deviant positions ensuring at least 2 standards between deviants
    deviantIndices = sort(randperm(numTrials, numDeviants));

    % Ensure at least 2 standards between consecutive deviants
    lastDeviantIndex = 0;

    % Assign deviant types while ensuring order and spacing
    for i = 1:numDeviants
        trialIdx = deviantIndices(i);
        while trialIdx <= numTrials && (trialIdx - lastDeviantIndex < 4 || strcmp(trialTypes{trialIdx}, 'D'))
            trialIdx = trialIdx + 1;
        end
        if trialIdx <= numTrials
            trialTypes{trialIdx} = 'D';
            lastDeviantIndex = trialIdx;
            dev1I(i,:)= trialIdx;
        end
    end

    % Fill in the rest with standard trials
    for i = 1:numTrials
        if isempty(trialTypes{i})
            trialTypes{i} = 'S';
        end
    end

    % Create a table to store parameters
    paramTable = table('Size', [numTrials, 13], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double','double', 'double'}, ...
        'VariableNames', {'Trial', 'ToneAmp', 'ToneFreq', 'ToneDur', 'ModAmp', 'ModFreq', 'ID_SweepTime', 'ID_F1', 'ID_F2', 'StimType', 'DeviantProbability', 'InterstimulusInterval', 'SorD'});

    % Add constant values to the table
    paramTable.Trial = (1:numTrials)';

    % Assign interstimulus interval to each row
%     paramTable.InterstimulusInterval = repmat(interstimulusInterval, numTrials, 1);
    paramTable.DeviantProbability(:) = deviantProbability1;
    paramTable.SorD(dev1I) = 2; %this will tell if stimulus delivered was Dev(2) or Standard (0)

    try
        for trial = 1:numTrials
            % Determine the current trial type
            currentTrialType = trialTypes{trial};

            % Set the parameters for the current trial
            switch currentTrialType
                case 'S'
                    currentParams = standardParams;
                case 'D'
                    currentParams = deviantParams1;
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
%             pause(interstimulusInterval / 1000); % Convert to seconds
               % Pause for the interstimulus interval
        if any(trial == deviantIndices)
            pause((interstimulusInterval+deviantParams1.ToneDur)/1000); % Convert to seconds
            paramTable.InterstimulusInterval(trial, :) = interstimulusInterval+deviantParams1.ToneDur;
        else
            pause((interstimulusInterval+standardParams.ToneDur)/1000); % Convert to seconds
            paramTable.InterstimulusInterval(trial, :) = interstimulusInterval+standardParams.ToneDur;
        end

            % Update the table with current trial parameters
            paramTable(trial, 2:10) = struct2table(currentParams);
        end

        % Print indices for deviant at the end
        fprintf('Indices of deviants: %s\n', num2str(find(strcmp(trialTypes, 'D'))));

        % Assign trial types to the table
%         paramTable.DeviantType = trialTypes';


%         % Generate a date and time string
%         dateTimeString = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

        % Save the entire table to a CSV file with the date and time in the name
%         writetable(paramTable, [ filename '.csv'], 'WriteVariableNames', true);
          writetable(paramTable, [ filename '.csv']);
        

        % Open the generated Excel file
%         winopen(['trial_parameters_' dateTimeString '.xlsx']);

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
end
