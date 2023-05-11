% Scripts to calculate relative spectral amplitudes

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

epoch = 4; % epoch length in seconds

% Loop across participants
for iSubject = 1:size(Subject_pool,1)
    
    % Load participant data
    EEGFile = fullfile(SubjectsDir,Subject_pool{iSubject});
    EEG = pop_loadset(EEGFile);
    
    % For the first run of subjects, claim some memory to store individual
    % subjects values for different bands
    % subject x channels
    if iSubject == 1
       % for mean
       Delta = zeros(size(Subject_pool,1),size(EEG.data,1)); 
       Theta = zeros(size(Subject_pool,1),size(EEG.data,1));
       Alpha = zeros(size(Subject_pool,1),size(EEG.data,1));
       Beta = zeros(size(Subject_pool,1),size(EEG.data,1));
       Gamma = zeros(size(Subject_pool,1),size(EEG.data,1));
       
       % for std
       Delta_std = zeros(size(Subject_pool,1),size(EEG.data,1)); 
       Theta_std = zeros(size(Subject_pool,1),size(EEG.data,1));
       Alpha_std = zeros(size(Subject_pool,1),size(EEG.data,1));
       Beta_std = zeros(size(Subject_pool,1),size(EEG.data,1));
       Gamma_std = zeros(size(Subject_pool,1),size(EEG.data,1));
    end
    
    % Claim some memory for the calculations for different bands
    % channels x epochs
    tmpDelta = zeros(size(EEG.data,1),size(EEG.data,3));
    tmpTheta = zeros(size(EEG.data,1),size(EEG.data,3));
    tmpAlpha = zeros(size(EEG.data,1),size(EEG.data,3));
    tmpBeta = zeros(size(EEG.data,1),size(EEG.data,3));
    tmpGamma = zeros(size(EEG.data,1),size(EEG.data,3));
    
    % Loop for each channel
    for ichan = 1:size(EEG.data,1)
        % Loop for each segment to have a better estimate
        for iepoch = 1:size(EEG.data,3)
            % for delta - from 1 to 4 Hz
            tmpDelta(ichan,iepoch) = relAmplitude(EEG.data(ichan,:,iepoch),EEG.srate,1,4,1,70);
            
            % for theta - from 4 to 8 Hz
            tmpTheta(ichan,iepoch) = relAmplitude(EEG.data(ichan,:,iepoch),EEG.srate,4,8,1,70);
            
            % for alpha - from 8 to 13 Hz
            tmpAlpha(ichan,iepoch) = relAmplitude(EEG.data(ichan,:,iepoch),EEG.srate,8,13,1,70);
            
            % for beta - from 13 to 30 Hz
            tmpBeta(ichan,iepoch) = relAmplitude(EEG.data(ichan,:,iepoch),EEG.srate,13,30,1,70);
            
            % for gamma - from 30 to 45 Hz
            tmpGamma(ichan,iepoch) = relAmplitude(EEG.data(ichan,:,iepoch),EEG.srate,30,70,1,70);
        end
    end
    
    % Assign the biweight estimate of the mean to avoid influence of
    % outliers
    [Delta(iSubject,:),Delta_std(iSubject,:)] = myBiweight(tmpDelta);
    [Theta(iSubject,:),Theta_std(iSubject,:)] = myBiweight(tmpTheta);
    [Alpha(iSubject,:),Alpha_std(iSubject,:)] = myBiweight(tmpAlpha);
    [Beta(iSubject,:),Beta_std(iSubject,:)] = myBiweight(tmpBeta);
    [Gamma(iSubject,:),Gamma_std(iSubject,:)] = myBiweight(tmpGamma);
    
end

%%
% select where to save data
cd(ResultsDir)

save('relative_amplitudes.mat',...
        'Delta','Theta',...
        'Alpha','Beta','Gamma')
