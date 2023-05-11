% Scripts to calculate Hjorth parameters

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
% Analysis settings

epoch = 4; % epoch length in seconds

% number of channels
nchan = 64;

%% 

% To store data

hjorth_activity = zeros(length(Subject_pool),nchan);
hjorth_mobility = zeros(length(Subject_pool),nchan);
hjorth_complexity = zeros(length(Subject_pool),nchan);

for iSubject = 1:size(Subject_pool,1)
    
    % Load participant data
    EEGFile = fullfile(SubjectsDir,Subject_pool{iSubject});
    EEG = pop_loadset(EEGFile);
    
    % Number of epochs
    n_eps = size(EEG.data,3);

    % Features to obtain
    features_s = zeros(n_eps,nchan,3);
    
    for n_epochs = 1:n_eps
        
        for chan = 1:nchan
            
            % Obtain time-series
            x_c=squeeze(EEG.data(chan,:,n_epochs));
            
            % Hjorth parameters
            
            % Activity
            activity = var(x_c);
   
            % Mobility
            mobility = std(diff(x_c))./std(x_c);

            % Complexity
            complexity = std(diff(diff(x_c)))./std(diff(x_c))./mobility;
 
            features_s(n_epochs,chan,1) = activity;
            features_s(n_epochs,chan,2) = mobility;
            features_s(n_epochs,chan,3) = complexity;
                
        end        
    end
    
    hjorth_activity(iSubject,:)=myBiweight(squeeze(features_s(:,:,1))');
    hjorth_mobility(iSubject,:)=myBiweight(squeeze(features_s(:,:,2))');
    hjorth_complexity(iSubject,:)=myBiweight(squeeze(features_s(:,:,3))');
    
end


%%
% select where to save data
cd(ResultsDir)
save('hjorth_parameters.mat','hjorth_activity','hjorth_mobility','hjorth_complexity');

%%
