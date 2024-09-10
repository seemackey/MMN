tic

clear
close all
clc


% Set up actx server/control
handles.RP = actxcontrol('RPco.x');
RP = handles.RP;

% Connect to the device and halt any ongoing processes
RP.ConnectRX6('USB', 1);
RP.Halt;
RP.ClearCOF;

% Load your rcx file and run it
RP.LoadCOF('C:\MMN-main\MMN_NewTiming.rcx');


% Directory for the text files
paramsDir = 'C:\MMN-main\';  %  your directory 
gimmefiggies = 1; % plots of the stimulus parameters as a check

% Define the parameters for the standard stimulus. 
% Stimtype #0 is unmodulated tone, #1 is AM, and #2 is FM
standardParams = struct('ToneAmp', .1, 'ToneFreq', 2000, 'ToneDur', 100, 'ModAmp', 1, 'ModFreq', 20, 'ID_SweepTime', 100, 'ID_F1', 2000, 'ID_F2', 12000, 'StimType', 1);

% Define the parameters for the deviant stimulus.
deviantParams1 = struct('ToneAmp', 0.1, 'ToneFreq', 2000, 'ToneDur', 100, 'ModAmp', 1, 'ModFreq', 80, 'ID_SweepTime', 100, 'ID_F1', 2000, 'ID_F2', 12000, 'StimType', 1);

% Define the probability of a deviant stimulus
deviantProbability1 = 0.1;

% Define interstimulus interval (in milliseconds)
interstimulusInterval = 624;

% Specify the number of trials
numTrials = 50;

%% CHECK THE PARAMETER SETTINGS WITH THIS PLOT
if gimmefiggies == 1
    % List of parameter filenames
    paramFiles = {'ToneAmp.txt', 'ToneFreq.txt', 'ToneDur.txt', 'ModAmp.txt', 'ModFreq.txt', ...
                  'FMSweepTime.txt', 'FM1.txt', 'FM2.txt', 'StimType.txt','ISI.txt'};
              
    % Corresponding parameter names for plot titles
    paramNames = {'Tone Amplitude', 'Tone Frequency', 'Tone Duration', 'Modulation Amplitude', ...
                  'Modulation Frequency', 'FM Sweep Time', 'FM1 Frequency', 'FM2 Frequency', ...
                  'Stimulus Type','Interstimulus Interval'};
    
    % Number of bins for histograms
    numBins = 20;
    
    % Create a figure for the histograms
    figure('Name', 'Parameter Histograms', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);
    
    % Loop through each parameter file
    for i = 1:length(paramFiles)
        % Read the parameter values from the text file
        paramValues = load(fullfile(paramsDir, paramFiles{i}));
    
        % Create a subplot for each parameter
        subplot(3, 4, i);  % Adjust the layout as needed 
        
        % Plot the histogram
        histogram(paramValues, numBins);
        
        % Set the title and labels
        title(paramNames{i});
        xlabel(paramNames{i});
        ylabel('Frequency');
    end
end


%  write to text files
generate_trials(standardParams, deviantParams1, deviantProbability1, interstimulusInterval, numTrials, paramsDir);

% run the circuit and record outputs
RP.Run;
RP.SoftTrg(1);

%  "TrialParameters" directory 
futureDir = fullfile(paramsDir, 'TrialParameters');
if ~exist(futureDir, 'dir')
    mkdir(futureDir);
end

% Define the list of parameter files and their corresponding names
paramFiles = {'ToneAmp.txt', 'ToneFreq.txt', 'ToneDur.txt', 'ModAmp.txt', 'ModFreq.txt', ...
              'FMSweepTime.txt', 'FM1.txt', 'FM2.txt', 'StimType.txt', 'ISI.txt'};
paramNames = {'Tone Amplitude', 'Tone Frequency', 'Tone Duration', 'Modulation Amplitude', ...
              'Modulation Frequency', 'FM Sweep Time', 'FM1 Frequency', 'FM2 Frequency', ...
              'Stimulus Type', 'Interstimulus Interval'};

% Initialize an empty matrix to hold all parameters
allParams = [];

% Loop through each parameter file and concatenate the data into allParams
for i = 1:length(paramFiles)
    % Read the parameter values from the text file
    paramValues = load(fullfile(paramsDir, paramFiles{i}));
    
    % Concatenate the parameter values as a new column in allParams
    allParams = [allParams, paramValues];
end

% Generate the output file names based on the current date and time
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
ev2FileName = fullfile(futureDir, [timestamp '.ev2']);
csvFileName = fullfile(futureDir, [timestamp '.csv']);
xlsxFileName = fullfile(futureDir, [timestamp '.xlsx']);

% Add column headers for CSV and Excel files
headers = strjoin(paramNames, ',');

% Save the data to a single .ev2 file (as plain text with space-separated values)
dlmwrite(ev2FileName, allParams, 'delimiter', ' ');

% Save the data to a single .csv file with headers
fid = fopen(csvFileName, 'w');
fprintf(fid, '%s\n', headers); % Write the headers
fclose(fid);
dlmwrite(csvFileName, allParams, '-append');

% Save the data to a single .xlsx file with headers
dataCell = [paramNames; num2cell(allParams)];
writecell(dataCell, xlsxFileName);

% pause(20)
% RP.Halt;

toc



function generate_trials(standardParams, deviantParams1, deviantProbability1, interstimulusInterval, numTrials, paramsDir)
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
        end
    end

    % Fill in the rest with standard trials
    for i = 1:numTrials
        if isempty(trialTypes{i})
            trialTypes{i} = 'S';
        end
    end

    % Open the text files for each parameter
    toneAmpFile = fopen(fullfile(paramsDir, 'ToneAmp.txt'), 'w');
    toneFreqFile = fopen(fullfile(paramsDir, 'ToneFreq.txt'), 'w');
    toneDurFile = fopen(fullfile(paramsDir, 'ToneDur.txt'), 'w');
    modAmpFile = fopen(fullfile(paramsDir, 'ModAmp.txt'), 'w');
    modFreqFile = fopen(fullfile(paramsDir, 'ModFreq.txt'), 'w');
    sweepTimeFile = fopen(fullfile(paramsDir, 'FMSweepTime.txt'), 'w');
    f1File = fopen(fullfile(paramsDir, 'FM1.txt'), 'w');
    f2File = fopen(fullfile(paramsDir, 'FM2.txt'), 'w');
    stimTypeFile = fopen(fullfile(paramsDir, 'StimType.txt'), 'w');
    isiFile = fopen(fullfile(paramsDir, 'ISI.txt'), 'w');

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
             % 
            fprintf(isiFile, '%f\n', interstimulusInterval);

  
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
end