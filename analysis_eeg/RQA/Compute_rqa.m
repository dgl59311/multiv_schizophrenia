% RQA using CPR toolbox and Pred. Maintenance Toolbox

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
epoch=4; % epoch length in seconds
nchan=64;

%%
for iSubject = 1:size(Subject_pool,1)

    % Load participant data
    EEGFile = fullfile(SubjectsDir,Subject_pool{iSubject});
    EEG = pop_loadset(EEGFile);
       
    % Dataset info
    n_eps = size(EEG.data,3);
    n_ch = size(EEG.data,1);
    sampl = EEG.srate;
	RQA_C = zeros(n_ch,13,n_eps);
    
    for ne = 1:n_eps
        
        signals = squeeze(EEG.data(:,:,ne));
        
        parfor (c = 1:n_ch, 5)
        
            x = double(signals(c,:));
            [~, lag, dim] = phaseSpaceReconstruction(x);
            recrqa = crqa(x',dim,lag,0.1,'nonormalize','nogui','rr','silent');
            RQA_C(c,:,ne) = recrqa;   
            
        end
    end
    
    RQA_data{iSubject} = RQA_C;
    fprintf(' %d', iSubject); 

end

%%

for i=1:size(Subject_pool,1)

    Subj_data_1=RQA_data{i};
    
    for nchan = 1:nchan
        
        % RQA metrics
        RQA_determinism(i,nchan)=myBiweight(squeeze(Subj_data_1(nchan,2,:))');
        RQA_mean_diag(i,nchan)=myBiweight(squeeze(Subj_data_1(nchan,3,:))');
        RQA_max_diag(i,nchan)=myBiweight(squeeze(Subj_data_1(nchan,4,:))');
        RQA_entropy(i,nchan)=myBiweight(squeeze(Subj_data_1(nchan,5,:))');
        RQA_laminarity(i,nchan)=myBiweight(squeeze(Subj_data_1(nchan,6,:))');
        RQA_trappingtime(i,nchan)=myBiweight(squeeze(Subj_data_1(nchan,7,:))');
        RQA_max_vert(i,nchan)=myBiweight(squeeze(Subj_data_1(nchan,8,:))');     
        RQA_rte(i,nchan)=myBiweight(squeeze(Subj_data_1(nchan,11,:))');
        
    end
end

%%
% select where to save data
cd(ResultsDir)
save('rqa_pred.mat','RQA_determinism',...
     'RQA_mean_diag','RQA_max_diag','RQA_entropy',...
     'RQA_laminarity','RQA_trappingtime','RQA_max_vert','RQA_rte')
    
%%
