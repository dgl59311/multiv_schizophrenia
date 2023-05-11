% Scripts to calculate Higuchi and Katz FD     

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

for iSubject = 1:size(Subject_pool,1)    
    
    % Load participant data
    EEGFile = fullfile(SubjectsDir,Subject_pool{iSubject});
    EEG = pop_loadset(EEGFile);
    
    % Filter the data for different frequency bands
    Delta = pop_eegfiltnew(EEG, 1, 4, [], 0, [], 0);
    Theta = pop_eegfiltnew(EEG, 4, 8, [], 0, [], 0);
    Alpha = pop_eegfiltnew(EEG, 8, 13, [], 0, [], 0);
    Beta = pop_eegfiltnew(EEG, 13, 30, [], 0, [], 0);
    Gamma = pop_eegfiltnew(EEG, 30, 70, [], 0, [], 0);

    
    if iSubject == 1
        
       % HFD
       HFD_Delta = zeros(size(Subject_pool,1),size(EEG.data,1)); 
       HFD_Theta = zeros(size(Subject_pool,1),size(EEG.data,1));
       HFD_Alpha = zeros(size(Subject_pool,1),size(EEG.data,1));
       HFD_Beta = zeros(size(Subject_pool,1),size(EEG.data,1));
       HFD_Gamma = zeros(size(Subject_pool,1),size(EEG.data,1));
       
       % Katz
       KFD_Delta = zeros(size(Subject_pool,1),size(EEG.data,1)); 
       KFD_Theta = zeros(size(Subject_pool,1),size(EEG.data,1));
       KFD_Alpha = zeros(size(Subject_pool,1),size(EEG.data,1));
       KFD_Beta = zeros(size(Subject_pool,1),size(EEG.data,1));
       KFD_Gamma = zeros(size(Subject_pool,1),size(EEG.data,1));
       
    end
    
    % Claim some memory for the calculations for different bands
    % channels x epochs
    % for HFD
    
    tmp_HFD_Delta = zeros(size(EEG.data,1),size(EEG.data,3));
    tmp_HFD_Theta = zeros(size(EEG.data,1),size(EEG.data,3));
    tmp_HFD_Alpha = zeros(size(EEG.data,1),size(EEG.data,3));
    tmp_HFD_Beta = zeros(size(EEG.data,1),size(EEG.data,3));
    tmp_HFD_Gamma = zeros(size(EEG.data,1),size(EEG.data,3));
    
    % for KFD
    tmp_KFD_Delta = zeros(size(EEG.data,1),size(EEG.data,3));
    tmp_KFD_Theta = zeros(size(EEG.data,1),size(EEG.data,3));
    tmp_KFD_Alpha = zeros(size(EEG.data,1),size(EEG.data,3));
    tmp_KFD_Beta = zeros(size(EEG.data,1),size(EEG.data,3));
    tmp_KFD_Gamma = zeros(size(EEG.data,1),size(EEG.data,3));
    
     
    % Loop for each channel
    for ichan = 1:size(EEG.data,1)
        % Loop for each segment to have a better estimate
        for iepoch = 1:size(EEG.data,3)
            % Delta
            tmp_HFD_Delta(ichan,iepoch) = Higuchi_FD(Delta.data(ichan,:,iepoch), 25);
            tmp_KFD_Delta(ichan,iepoch) = Katz_FD(Delta.data(ichan,:,iepoch));
            
            % Theta
            tmp_HFD_Theta(ichan,iepoch) = Higuchi_FD(Theta.data(ichan,:,iepoch), 25);
            tmp_KFD_Theta(ichan,iepoch) = Katz_FD(Theta.data(ichan,:,iepoch));
            
            % Alpha
            tmp_HFD_Alpha(ichan,iepoch) = Higuchi_FD(Alpha.data(ichan,:,iepoch), 25);
            tmp_KFD_Alpha(ichan,iepoch) = Katz_FD(Alpha.data(ichan,:,iepoch));
            
            % Beta
            tmp_HFD_Beta(ichan,iepoch) = Higuchi_FD(Beta.data(ichan,:,iepoch), 25);
            tmp_KFD_Beta(ichan,iepoch) = Katz_FD(Beta.data(ichan,:,iepoch));
            
            % Gamma
            tmp_HFD_Gamma(ichan,iepoch) = Higuchi_FD(Gamma.data(ichan,:,iepoch), 25);
            tmp_KFD_Gamma(ichan,iepoch) = Katz_FD(Gamma.data(ichan,:,iepoch));

        end
    end
    
    % Assign the biweight estimate of the mean to avoid influence of
    % outliers
    % for amplitude total power
    [HFD_Delta(iSubject,:), ~] = myBiweight(tmp_HFD_Delta);
    [HFD_Theta(iSubject,:), ~] = myBiweight(tmp_HFD_Theta);
    [HFD_Alpha(iSubject,:), ~] = myBiweight(tmp_HFD_Alpha);
    [HFD_Beta(iSubject,:), ~] = myBiweight(tmp_HFD_Beta);
    [HFD_Gamma(iSubject,:), ~] = myBiweight(tmp_HFD_Gamma);
    
    % for mean of the envelope
    [KFD_Delta(iSubject,:), ~] = myBiweight(tmp_KFD_Delta);
    [KFD_Theta(iSubject,:), ~] = myBiweight(tmp_KFD_Theta);
    [KFD_Alpha(iSubject,:), ~] = myBiweight(tmp_KFD_Alpha);
    [KFD_Beta(iSubject,:), ~] = myBiweight(tmp_KFD_Beta);
    [KFD_Gamma(iSubject,:), ~] = myBiweight(tmp_KFD_Gamma);
    
    n_eps = size(EEG.data,3);

end

%%
% select where to save data
cd(ResultsDir)

save('KFD.mat', 'KFD_Delta', 'KFD_Theta', 'KFD_Alpha',...
     'KFD_Beta', 'KFD_Gamma')
 
save('HFD.mat', 'HFD_Delta', 'HFD_Theta', 'HFD_Alpha',...
     'HFD_Beta', 'HFD_Gamma')
 
%%

