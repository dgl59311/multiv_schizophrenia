% Scripts to calculate the Lempel-Ziv Complexity (LZC)

% Directory management
clc; 
clear; 
close all;

% get current directory
% use .../amplitude_features
CurrDir = pwd;                                    

% Paths to the data
SubjectsDir = uigetdir([],'Path to the preprocessed data');

% Subjects pre-processed data files
Subject_data = dir(fullfile(SubjectsDir,'*.mat'));

% The pool of subjects
Subject_pool = {Subject_data(:).name}';

% Path to save .mat results
ResultsDir = uigetdir([], 'Save .mat files with results');

%% 

% Loop across subjects
for iSubject = 1:size(Subject_pool,1)
    
    % Load participant data
    EEGFile = fullfile(SubjectsDir,Subject_pool{iSubject});
    EEG = pop_loadset(EEGFile);
     
    % For the first run of subjects, claim some memory to store individual
    % subjects values for the Lempel-Ziv Complexity (LZC)
    % subject x channels
    if iSubject == 1
       LZC_exhaustive = zeros(size(Subject_pool,1),size(EEG.data,1)); 
       LZC_primitive = zeros(size(Subject_pool,1),size(EEG.data,1)); 
    end
    
    % Claim some memory for the LZC calculation
    % channels x epochs
    tmpLZC_exhaustive = zeros(size(EEG.data,1),size(EEG.data,3));
    tmpLZC_primitive = zeros(size(EEG.data,1),size(EEG.data,3));
    
    % Loop for each channel
    for ichan = 1:size(EEG.data,1)
        % Loop for each segment to have a better estimate
        for iepoch = 1:size(EEG.data,3)
            
            % Calculate the LZC for each channel and epoch
            m_data = median(EEG.data(ichan,:,iepoch));
            logical_data = double(EEG.data(ichan,:,iepoch) > m_data);
            
            [C, ~, ~] = calc_lz_complexity(logical_data, 'exhaustive', true);
            tmpLZC_exhaustive(ichan,iepoch) = C;
            
            [C_2, ~, ~] = calc_lz_complexity(logical_data, 'primitive', true);
            tmpLZC_primitive(ichan,iepoch) = C_2;
            
        end
    end
    
    % Assign the biweight estimate of the mean to avoid influence of
    % outliers
    [LZC_exhaustive(iSubject,:),~] = myBiweight(tmpLZC_exhaustive);
    [LZC_primitive(iSubject,:),~] = myBiweight(tmpLZC_primitive);
    
end


%%
% select where to save data
cd(ResultsDir)
save('LZC.mat', 'LZC_exhaustive','LZC_primitive');
%%




