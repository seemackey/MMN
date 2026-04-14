function generate_trials2Devs(standardParams, deviantParams1, deviantParams2, deviantProbability1, deviantProbability2, interstimulusInterval, numTrials, paramsDir)

    % Calculate number of trials for deviants 
    numDeviants1 = ceil(numTrials * deviantProbability1);
    numDeviants2 = ceil(numTrials * deviantProbability2);

    % Ensure deviants don't exceed 10% of total trials
    
    if numDeviants1 > deviantProbability1 * numTrials
        numDeviants1 = floor(deviantProbability1 * numTrials);
    end

    if numDeviants2 > deviantProbability2 * numTrials
        numDeviants2 = floor(deviantProbability2 * numTrials);
    end
totalDuration = 0;
%     % Initialize cell array to store trial types
%     trialTypes = cell(1, numTrials);
% 
% 
%     minimumSpacing = 3:1:10;
% 
%     % Initialize the first trial index
%     trialIdx = 1;
%     deviantIndices = [];
%     for i = 1:numDeviants
%         % Add a random spacing from the minimumSpacing range
%         spacing = randi([min(minimumSpacing), max(minimumSpacing)], 1, 1); 
%         
%         % Compute the next trial index by adding the spacing
%         trialIdx = trialIdx + spacing;
% 
%         % Assign the deviant to the calculated trial index
%         trialTypes{trialIdx} = 'D';
%         deviantIndices(i) = trialIdx;
%     end
%     for i = 1:numTrials
%         if isempty(trialTypes{i})
%             trialTypes{i} = 'S';
%         end
%     end



 % Initialize cell array to store trial types
    trialTypes = cell(1, numTrials);
      trialTypes2 = cell(1, numTrials);
    % Ensure at least 2 standards between consecutive deviants
    lastD1Index = 0;
    lastD2Index = 0;



%%%% If want random number of standards between Devs - Comment/Uncomment - then save
                                                %%Randomly assign deviant positions ensuring at least 2 standards between deviants
                                                deviant1Indices = sort(randperm(numTrials+35, numDeviants1));
                                                deviant2Indices = sort(randperm(numTrials+35, numDeviants2));
                                                 RR=[4 5 6 7 9];
                                                % Assign deviant types while ensuring order and spacing
                                                for i = 1:max(numDeviants1, numDeviants2)
                                                    if i <= numDeviants1
                                                        trialIdx = deviant1Indices(i);
                                                         x=randsample(RR,1); %randomly selects number of standards between the previous dev and this one
                                                        while trialIdx <= numTrials && (trialIdx - lastD2Index < x || strcmp(trialTypes{trialIdx}, 'D2'))
                                                            trialIdx = trialIdx + 1;
                                                        end
                                                        if trialIdx <= numTrials
                                                            trialTypes{trialIdx} = 'D1';
                                                            lastD1Index = trialIdx;
                                                             deviant1Indices_real(i,:)=trialIdx;
                                                        end
                                                    end

                                                    if i <= numDeviants2
                                                        trialIdx = deviant2Indices(i);
                                                         x1=randsample(RR,1);%randomly selects number of standards between the previous dev and this one
                                                        while trialIdx <= numTrials && (trialIdx - lastD1Index < x1 || strcmp(trialTypes{trialIdx}, 'D1'))
                                                            trialIdx = trialIdx + 1;
                                                        end
                                                        if trialIdx <= numTrials
                                                            trialTypes{trialIdx} = 'D2';
                                                            lastD2Index = trialIdx;
                                                             deviant2Indices_real(i,:)=trialIdx;
                                                        end
                                                    end
                                                end

                                                E= find(deviant1Indices_real ==0 | deviant1Indices_real>=numTrials);
                                                deviant1Indices_real(E)=[];
                                                deviant1Indices_real=unique(deviant1Indices_real);
                                                [trialTypes2{deviant1Indices_real}] = deal('D1');
                                                clear E;E= find(deviant2Indices_real ==0 | deviant2Indices_real>=numTrials);
                                                deviant2Indices_real(E)=[];
                                                deviant2Indices_real=unique(deviant2Indices_real);
                                                [trialTypes2{deviant2Indices_real}] = deal('D2');


    %%%% If want 4 standards between Devs - temporal regularity but get more Devs - Comment/Uncomment - then save
                            % deviant1Indices = sort(randperm(numTrials, numDeviants1));
                            % deviant2Indices = sort(randperm(numTrials, numDeviants2));
                            % % Assign deviant types while ensuring order and spacing
                            % for i = 1:max(numDeviants1, numDeviants2)
                            %     if i <= numDeviants1
                            %         trialIdx = deviant1Indices(i);
                            %         while trialIdx <= numTrials && (trialIdx - lastD2Index < 4 || strcmp(trialTypes{trialIdx}, 'D2'))
                            %             trialIdx = trialIdx + 1;
                            %         end
                            %         if trialIdx <= numTrials
                            %             trialTypes{trialIdx} = 'D1';
                            %             lastD1Index = trialIdx;
                            %             deviant1Indices_real(i,:)=trialIdx;
                            %         end
                            %     end
                            % 
                            %     if i <= numDeviants2
                            %         trialIdx = deviant2Indices(i);
                            %         while trialIdx <= numTrials && (trialIdx - lastD1Index < 4 || strcmp(trialTypes{trialIdx}, 'D1'))
                            %             trialIdx = trialIdx + 1;
                            %         end
                            %         if trialIdx <= numTrials
                            %             trialTypes{trialIdx} = 'D2';
                            %             lastD2Index = trialIdx;
                            %             deviant2Indices_real(i,:)=trialIdx;
                            %         end
                            %     end
                            % end
                            % 
                            % E= find(deviant1Indices_real>=numTrials);
                            % deviant1Indices_real(E)=[];
                            % deviant1Indices_real=unique(deviant1Indices_real);
                            % [trialTypes2{deviant1Indices_real}] = deal('D1');
                            % clear E;E= find( deviant2Indices_real>=numTrials);
                            % deviant2Indices_real(E)=[];
                            % deviant2Indices_real=unique(deviant2Indices_real);
                            % [trialTypes2{deviant2Indices_real}] = deal('D2');



      % Fill in the rest with standard trials
    for i = 1:numTrials
        if isempty(trialTypes2{i})
            trialTypes2{i} = 'S';
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
    totalDur= fopen(fullfile(paramsDir, 'totalDur.txt'), 'w');
    Dev1Prob= fopen(fullfile(paramsDir, 'Dev1Prob.txt'), 'w');
    Dev2Prob= fopen(fullfile(paramsDir, 'Dev2Prob.txt'), 'w');
    %% Calculate the expected run time considering different durations for standard and deviant trials
    try
        for trial = 1:numTrials
            % Determine the current trial type
            currentTrialType = trialTypes2{trial};

            % Set the parameters for the current trial
            switch currentTrialType
                case 'S'
                    currentParams = standardParams;
                    deviantFlag = 0; % Standard trial
                      totalDuration = totalDuration + (interstimulusInterval + standardParams.ToneDur);
                case 'D1'
                    currentParams = deviantParams1;
                    deviantFlag = 1; % Deviant trial
                     totalDuration = totalDuration + (interstimulusInterval + deviantParams1.ToneDur);
                case 'D2'
                    currentParams = deviantParams2;
                    deviantFlag = 2; % Deviant trial
                     totalDuration = totalDuration + (interstimulusInterval + deviantParams2.ToneDur);
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

        % Print indices for deviant1 and deviant2 together at the end
        fprintf('Indices of deviants1: %s\n', num2str(find(strcmp(trialTypes2, 'D1'))));
         fprintf('Indices of deviants2: %s\n', num2str(find(strcmp(trialTypes2, 'D2'))));
      

    catch exception
        % Display the error message
        disp('An error occurred: something stopped the loop before it could finish. Closing files and saving data.');
        disp(exception.message);
    end

% Convert to seconds
    totalDuration = totalDuration / 1000;
 fprintf(totalDur, '%d\n', totalDuration);
  fprintf( Dev1Prob, '%d\n', length(deviant1Indices_real)/numTrials);
  fprintf( Dev2Prob, '%d\n', length(deviant2Indices_real)/numTrials);
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
     fclose(totalDur); % Close the new totalDuration file
     fclose(Dev1Prob); % Close the new Probability file
     fclose(Dev2Prob); % Close the new Probability file
end