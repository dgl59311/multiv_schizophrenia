classdef eeg_conn_matrix
    
    properties
        data=[];
    end
    
    % The method for acquiring the data path
    % Adding one persistent variable sub_path to remember the data path.
    methods(Static)
        function path = data_path
            persistent sub_path
            if isempty(sub_path)
                % Paths to the data
                path.SubjectsDir = uigetdir([],'Path to the data of Subjects');
                sub_path = path.SubjectsDir;
            else
                path.SubjectsDir = uigetdir(sub_path,'Path to the data of Subjects');
                sub_path = path.SubjectsDir;
            end
            % Subjects pre-processed data files
            Subject_data = dir(fullfile(path.SubjectsDir,'*preprocess.mat'));
            % The pool of subjects
            path.Subject_pool = {Subject_data(:).name}';
        end        
    end
    
    methods
        % The function to build the object for the following analysis
        function obj = eeg_conn_matrix
            obj.data=[];
        end
        
        % The function to perform the frequency analysis
        function  freq=freq_ana(obj,EEG,cfg_freq)
            % The function transfrom the data format from EEGlab to fieldtrip
            data_trans=eeglab2fieldtrip(EEG,'preprocessing','none');
            freq      = ft_freqanalysis(cfg_freq, data_trans);
        end
        
        % The function to perform the connectivity analysis
        function freq_conn = connectivity_ana(obj,cfg_conn,freq)
            
            freq_conn = ft_connectivityanalysis(cfg_conn, freq);
        end
        
        % The function to perform the network analysis
        function network_output = network_ana(obj,cfg_network,freq_conn)            
            network_full = ft_networkanalysis(cfg_network,freq_conn);
            network_output=eval(['network_full.',cfg_network.method]);
        end
    end
    
end

