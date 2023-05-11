% this script is used to apply the CSD to the resting SZ data
%   Kayser J, Tenke CE. Principal components analysis of Laplacian 
%   waveforms as a generic method for identifying ERP generator patterns: 
%   I. Evaluation with auditory oddball tasks. Clin Neurophysiol 2006; 117: 348–68.
% Janir Ramos on 31/01/2020

clc;
clear;
CurrDir = pwd;                                  % gets current directory
mkdir(pwd,'CSD Results');                       % creates a directory to save data
SaveDir = strcat(pwd,'\CSD Results');           % directory where to save data
% DataDir = ; % where data is

% load the transformation matrices
load('64-BiosemiCSD_TransfMatrixes.mat');       % H and G
Group = {'SZPatients';'Controls'};

for iGroup = 1:size(Group,1)
    
    % Group directory with the data
    GroupDir = strcat(DataDir,'\',char(Group{iGroup}));
    cd(GroupDir);
    % Creates the directory to save data for each group
    mkdir(SaveDir,Group{iGroup});
    
    % Subjects mat files
    Subject_mat = dir('*.mat*');
    Subject_pool = {Subject_mat(:).name}';
    
    for iSubject = 1:size(Subject_pool,1)
        % load the EEG data of the subject
        load(char(Subject_pool(iSubject)));
        % claim data memory to speed computations
        data = single(zeros(size(EEG.data)));
        
        for itrial = 1:EEG.trials
            EEGtmp = single(squeeze(EEG.data(:,:,itrial)));         % reduce data precision to reduce memory
            EEGCSD = CSD(EEGtmp,G,H);                               % Compute the CSD
            EEGCSD_z = zscore(EEGCSD,1,2);                          % 0 mean and unit variance for connectivity
            data(:,:,itrial) = EEGCSD_z;                            % Assign data output
        end
        
        EEG.data = double(data);                                    % final CSD data
        
        cd(strcat(SaveDir,'\',char(Group{iGroup})));                % change to save directory
        save(char(Subject_pool{iSubject}),'EEG');              % save data
        cd(GroupDir);                                               % change to group directory
    end
    
end

