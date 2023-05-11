% Life and waiting times 

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

% Number of channels
nchan = 64;

%% 

% To store data
LifeTimes = zeros(5,length(Subject_pool),nchan);
WaitingTimes = zeros(5,length(Subject_pool),nchan);

for iSubject = 1:size(Subject_pool,1)
    
    % Load participant data
    EEGFile = fullfile(SubjectsDir,Subject_pool{iSubject});
    EEG = pop_loadset(EEGFile);

    % Time series
	signals = EEG.data;
    % Total data length
    largeset = length(signals);   

    for band = 1:5
        
        % Set bandpass filter
        [locutoff,hicutoff] = FiltLims(band);
        filtorder = 2*fix(EEG.srate/locutoff);
        b = fir1(filtorder, [locutoff, hicutoff]./(EEG.srate/2),'bandpass');
        
        % Store life and waiting times
        LT = [];
        WT = [];
    
        for es = 1:nchan
           
            x = signals(es,:);
            xf = filter(b,1,double(x));
            A = abs(hilbert(xf));
            
            % Use the median of the signal as threshols
            Th = median(A);
            
            % Find bursts above threshold
            EnvL = find(A>Th);
            dL = find(diff(EnvL)>1);
            LT = [dL length(EnvL)]-[0 dL];
            
            % Find bursts below threshold
            EnvW = find(A<Th);
            dW = find(diff(EnvW)>1);
            WT = [dW length(EnvW)]-[0 dW];
            
            % Calculate cumulative probability distributions
            [P95L] = cumpdf(1000*(LT/EEG.srate));  %ms
            [P95W] = cumpdf(1000*(WT/EEG.srate));  %ms
            
            % Store data
            LifeTimes(band,iSubject,es) = P95L;
            WaitingTimes(band,iSubject,es) = P95W;
     
            LT = [];
            WT = [];   
        end
    end
end

%%
% Data to save

lifetime_delta = squeeze(LifeTimes(1,:,:));
lifetime_theta = squeeze(LifeTimes(2,:,:));
lifetime_alpha = squeeze(LifeTimes(3,:,:));
lifetime_beta = squeeze(LifeTimes(4,:,:));
lifetime_gamma = squeeze(LifeTimes(5,:,:));

waitingtime_delta = squeeze(WaitingTimes(1,:,:));
waitingtime_theta = squeeze(WaitingTimes(2,:,:));
waitingtime_alpha = squeeze(WaitingTimes(3,:,:));
waitingtime_beta = squeeze(WaitingTimes(4,:,:));
waitingtime_gamma = squeeze(WaitingTimes(5,:,:));

%%
% select where to save data
cd(ResultsDir)
    
    save('lifetimes.mat',...
        'lifetime_delta','lifetime_theta',...
        'lifetime_alpha','lifetime_beta','lifetime_gamma')
    
    save('waitingtimes.mat',...
        'waitingtime_delta','waitingtime_theta',...
        'waitingtime_alpha','waitingtime_beta','waitingtime_gamma')   

%%


