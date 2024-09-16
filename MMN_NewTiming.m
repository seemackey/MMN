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




% Directory for the text files
paramsDir = 'C:\MMN-main\';  %  your directory 
gimmefiggies = 1; % plots of the stimulus parameters as a check

% Define the parameters for the standard stimulus. 
% Stimtype #0 is unmodulated tone, #1 is AM, and #2 is FM
standardParams = struct('ToneAmp', .08, 'ToneFreq', 1000, 'ToneDur', 100, 'ModAmp', 1, 'ModFreq', 20, 'ID_SweepTime', 100, 'ID_F1', 2000, 'ID_F2', 12000, 'StimType',0);

% Define the parameters for the deviant stimulus.
deviantParams1 = struct('ToneAmp', 0.04, 'ToneFreq', 8000, 'ToneDur', 100, 'ModAmp', 1, 'ModFreq', 80, 'ID_SweepTime', 100, 'ID_F1', 2000, 'ID_F2', 4000, 'StimType',0);

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



%  "TrialParameters" directory 
futureDir = fullfile(paramsDir, 'TrialParameters');
if ~exist(futureDir, 'dir')
    mkdir(futureDir);
end

% Define the list of parameter files and their corresponding names
paramFiles = {'ToneAmp.txt', 'ToneFreq.txt', 'ToneDur.txt', 'ModAmp.txt', 'ModFreq.txt', ...
              'FMSweepTime.txt', 'FM1.txt', 'FM2.txt', 'StimType.txt', 'ISI.txt', 'Deviant.txt'};
paramNames = {'Tone Amplitude', 'Tone Frequency', 'Tone Duration', 'Modulation Amplitude', ...
              'Modulation Frequency', 'FM Sweep Time', 'FM1 Frequency', 'FM2 Frequency', ...
              'Stimulus Type', 'Interstimulus Interval', 'Deviant'}; % Add 'Deviant' as a new column name

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

% Pause to allow the user to review the parameters
pause(2);

% Ask the user if the parameters are okay and if they are ready to proceed
while true
    userResponse = input('Are the parameters correct? Type "yes" to proceed or "no" to review again: ', 's');
    if strcmpi(userResponse, 'yes')
        % Load the rcx file and run it
        RP.LoadCOF('C:\MMN-main\MMN_NewTiming.rcx');
        RP.Run;
        disp('Circuit is now running...');
        break; % Exit the loop and proceed
    elseif strcmpi(userResponse, 'no')
        disp('Please review the parameters and run the script again.');
        return; % Exit the script
    else
        disp('Invalid response. Please type "yes" or "no".');
    end
end

% Notify the user about the expected run time
totalDuration = (numTrials * (interstimulusInterval + standardParams.ToneDur)) / 1000; % Convert from ms to seconds
fprintf('The circuit is expected to run for approximately %.2f seconds; type "RP.Halt" to halt the circuit.\n', totalDuration);

% Initialize the timer
startTime = tic;


% If the loop completes without manual intervention, halt the circuit
if toc(startTime) >= totalDuration
    RP.Halt;
    disp('The expected circuit duration has elapsed. The circuit has been halted.');
end






