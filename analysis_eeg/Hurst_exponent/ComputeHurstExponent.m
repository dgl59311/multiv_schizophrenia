% Scripts to calculate the Hurst Exponent (HE)

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
    
    % Get and load the subject
    EEGFile = fullfile(SubjectsDir,Subject_pool{iSubject});
    load(EEGFile);
    
    % For the first run of subjects, claim some memory to store individual
    % subjects values for the Hurst Exponent HE
    % subject x channels
    if iSubject == 1
       HE = zeros(size(Subject_pool,1),size(EEG.data,1));
       HE_std = zeros(size(Subject_pool,1),size(EEG.data,1));
    end
    
    % Claim some memory for the HE calculation
    % channels x epochs
    tmpHE = zeros(size(EEG.data,1),size(EEG.data,3));
    
    % Loop for each channel
    for ichan = 1:size(EEG.data,1)
        % Loop for each segment to have a better estimate
        for iepoch = 1:size(EEG.data,3)
            
            % Calculate the HE for each channel and epoch
            tmpHE(ichan,iepoch) = estimate_hurst_exponent(EEG.data(ichan,:,iepoch));
            
        end
    end
    
    % Assign the biweight estimate of the mean to avoid influence of
    % outliers
    [HE(iSubject,:),HE_std(iSubject,:)] = myBiweight(tmpHE);

end

%%

save('HE.mat','hurst_exponent')

%%