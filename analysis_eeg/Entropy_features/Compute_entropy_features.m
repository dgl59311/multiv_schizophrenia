% Scripts to calculate entropies

% Directory management
clc; 
clear; 
close all;
%%
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

% Number of EEG channels
nchan = 64;

%%

% To store data

sample_entropy = zeros(length(Subject_pool),nchan);
approximate_entropy = zeros(length(Subject_pool),nchan);
correlation_dimension = zeros(length(Subject_pool),nchan);
lyapunov_exponent = zeros(length(Subject_pool),nchan);

for iSubject = 1:size(Subject_pool,1)
   
    % Load participant data
    EEGFile = fullfile(SubjectsDir,Subject_pool{iSubject});
    EEG = pop_loadset(EEGFile);
    
    % Number of available epochs
    n_epochs = size(EEG.data,3);
    
    % Store entropy features
    features_s = zeros(n_epochs,nchan,4);
    
    for n_ep = 1:n_epochs
        
        signals = squeeze(EEG.data(:,:,n_ep));

        parfor (chan = 1:nchan,5)
        %for chan = 1:nchan
            
            % Get time series
            x_c = signals(chan,:);
            
            % Embedding dimension = 3
            [~, lag, dim] = phaseSpaceReconstruction(x_c);
            get_samp_en = sampen(x_c,dim,0.2);
            get_approx_entr = approximateEntropy(x_c,lag,dim);
            get_corr_dim = correlationDimension(x_c,lag,dim);
            try
                % Calculate the LyapExp for each channel and epoch
                % https://www.sciencedirect.com/science/article/pii/S2352914819303090
                [get_lyap, ~] = lyaprosen(x_c,0,0); close all;
                
            catch
                get_lyap = nan;
            end
            % Save entropy features
            ENS = [get_samp_en;get_approx_entr;get_corr_dim;get_lyap];
            features_s(n_ep,chan,:) = ENS;
                
        end    
    end
    
    % Store data
    
    sample_entropy(iSubject,:) = myBiweight(squeeze(features_s(:,:,1))');
    approximate_entropy(iSubject,:) = myBiweight(squeeze(features_s(:,:,2))');
    correlation_dimension(iSubject,:) = myBiweight(squeeze(features_s(:,:,3))');
    lyapunov_exponent(iSubject,:) = nanmean(squeeze(features_s(:,:,4)));
    
end


%%
% select where to save data
cd(ResultsDir)

save('entropies_fullband.mat',...
     'sample_entropy','approximate_entropy', 'correlation_dimension', 'lyapunov_exponent');

