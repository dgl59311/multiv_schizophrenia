clc
clearvars -except sub_path
% add fieldtrip to path
%%
%%% frequency parameters
% more can be found from the following link:
% http://www.fieldtriptoolbox.org/reference/ft_freqanalysis/
cfg_freq        = [];
cfg_freq.method = 'mtmfft';
cfg_freq.output = 'fourier';
cfg_freq.taper  = 'dpss';
cfg_freq.foi = 1:1:100;
cfg_freq.tapsmofrq = 1;

%%
%%% connectivity parameters
% Calculating different types of connectivity
% more can be found from the following link:
% http://www.fieldtriptoolbox.org/reference/ft_connectivityanalysis/

% Method 1 coherence
%cfg_conn           = [];
%cfg_conn.method    = 'coh';
%cfg_conn.complex = 'absimag';

%
% Method 2 phase locking value
cfg_conn           = [];
cfg_conn.method    = 'plv';

% Method 3 dtf
%cfg_conn           = [];
%cfg_conn.method    = 'dtf';

%%
cfg_network           = [];
cfg_network.method    = 'clustering_coef';
cfg_network.parameter = [cfg_conn.method 'spctrm']  ;
%%
% Getting data path and data object to perform analysis
data_info = eeg_conn_matrix.data_path;
data_obj = eeg_conn_matrix;

%% Run connectivity analysis
for i=1:size(data_info.Subject_pool,1)
    
    EEGFile = fullfile(data_info.SubjectsDir,data_info.Subject_pool{i});
    load(EEGFile);
    
    % Performing frequency analysis
    freq = data_obj.freq_ana(EEG,cfg_freq);
    % Performing connectivity analysis
    eval(['data_obj.data(i,1).',cfg_conn.method,' = data_obj.connectivity_ana(cfg_conn,freq);'])
    
end

%% Run network analysis
% Allocate data
method_nx = {'clustering_coef', 'betweenness'};
el_betw_ = zeros(5, length(data_info.Subject_pool), length(freq.label));
el_cc_ = zeros(5, length(data_info.Subject_pool), length(freq.label));
el_nodestr_ = zeros(5, length(data_info.Subject_pool), length(freq.label));

%%
for ii =1:size(data_info.Subject_pool,1)
    %Performing network analysis
    for jj = 1:length(method_nx)
        cfg_network.method = method_nx{jj};
        data_obj.data(ii,1).net = eval(['data_obj.network_ana(cfg_network,data_obj.data(ii,1).',cfg_conn.method,')']);
        for kk = 1:6
            [l_, h_] = FiltLims(kk);    
            if jj == 1
                el_cc_(kk, ii, :) = mean(data_obj.data(ii).net(:, l_:h_), 2);
            else
                el_betw_(kk, ii, :) = mean(data_obj.data(ii).net(:, l_:h_), 2);
            end
            if cfg_conn.method == 'coh'
                el_nodestr_(kk, ii, :) = sum(mean(data_obj.data(ii).coh.cohspctrm(:, :, l_:h_), 3), 2);
            elseif cfg_conn.method == 'plv'
                el_nodestr_(kk, ii, :) = sum(mean(data_obj.data(ii).plv.plvspctrm(:, :, l_:h_), 3), 2);
            else 
                 el_nodestr_(kk, ii, :) = sum(mean(data_obj.data(ii).dtf.dtfspctrm(:, :, l_:h_), 3), 2);               
            end
        end
    end
end
%% Write csv tables with results
% Specify the correct algorithm and group 
bands = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
alg = cfg_conn.method;
f_names = fieldnames(data_obj.data(1));
group = 'patients';

for i = 1:5
    
    btw = array2table(squeeze(el_betw_(i, :, :)), 'VariableNames', data_obj.data(1).(f_names{1}).label);
    writetable(btw, [bands{i} '_' alg '_betw_' group '.csv'])
    nodestr = array2table(squeeze(el_nodestr_(i, :, :)), 'VariableNames', data_obj.data(1).(f_names{1}).label);
    writetable(nodestr, [bands{i} '_' alg '_cc_' group '.csv'])
    cc = array2table(squeeze(el_cc_(i, :, :)), 'VariableNames', data_obj.data(1).(f_names{1}).label);
    writetable(cc, [bands{i} '_' alg '_str_' group '.csv'])
    
end
