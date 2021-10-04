%% PATHS AND SETTING UP FIELDTRIP AND PATHS 
clear classes;
master_dir = 'W:\PhD\msterdir';
main_path = 'W:\PhD\PatternGlareData\participants\participant_';
results_dir = 'W:\PhD\PatternGlareData\Results';
%rmpath 'W:\PhD\MatlabPlugins\spm8';
%addpath 'W:\PhD\MatlabPlugins\spm12';
addpath 'W:\PhD\MatlabPlugins\fieldtrip-20210906';
ft_defaults;
cd(master_dir);

%% WHAT TYPE OF EXPERIMENT(s) ARE WE RUNNING?
experiment_types = {'onsets-2-8-explicit'};   
desired_design_mtxs = {'no-factor'}; %
start_latency = -0.20;
end_latency = 3.9;
n_participants = 40

%% SHALL WE APPLY A ROI, IF SO HOW?
region_of_interest = 0;
roi_applied = 'two-tailed';
weight_roi = 0;
roi_to_apply = 0;
grand_avgs = {};
cfg_parameter = {'avg','thin','thick','med'};

%% GENERATE ERPS AND COMPUTE CONFIDENCE INTERVALS
generate_erps = 1;
weight_erps = 1; % weights based on quartiles
weighting_factor = 0.00; % weights based on quartiles

%% CHOOSE THE TYPE OF ANALYSIS EITHER 'frequency_domain' or 'time_domain'
type_of_analysis = 'time_domain';

if strcmp(type_of_analysis, 'frequency_domain')
    disp('RUNNING A FREQUENCY-DOMAIN ANALYSIS');
    compute_frequency_data = 0; % compute the freq data per particpant else load
    frequency_type = 'fourier'; % compute inter trial coherence
    run_mua = 0; % run a MUA in the frequnecy domain?
    analysis_on_aggr_data = 1; % analysis on the aggregate power data?
elseif strcmp(type_of_analysis, 'time_domain')
    disp('RUNNING A TIME-DOMAIN ANALYSIS');
end

data_file = 'time_domain_mean_intercept_onsets_2_3_4_5_6_7_8_grand-average.mat';
partition.is_partition = 0;
partition.partition_number = 0;

[data, participant_order_1] = load_postprocessed_data(main_path, n_participants, ...
               data_file, partition);


for index = 1:length(cfg_parameter)
    cfg = [];
    cfg.channel   = 'all';
    cfg.latency   = 'all';
    cfg.parameter = cfg_parameter{index};

    grand_avgs{index} = ft_timelockgrandaverage(cfg, data{:});
end

cfg = [];
cfg.linewidth = 2;
cfg.showlegend    = {'thin','thick','med'};
cfg.baseline = [2.8 3.0]
cfg.channel = 'A22';

ft_singleplotER(cfg, grand_avgs{1}, grand_avgs{2}, grand_avgs{3});

%% load post-processed fildtrip data
function [ft_regression_data, participant_order] = ...
    load_postprocessed_data(main_path, n_participants, filename, partition)

    ft_regression_data = {};  
    participant_order = {};

    idx_used_for_saving_data = 1;
    for i=1:n_participants
        disp(strcat('LOADING PARTICIPANT...', int2str(i)));
        participant_main_path = strcat(main_path, int2str(i));

        if exist(participant_main_path, 'dir')
            cd(participant_main_path);
            
            if isfile(filename)
                load(filename);
            else
                continue;
            end
            
            ft.label = data.label;
            ft.time = data.time{1};
            ft.trialinfo = [1];
            ft.elec = data.elec;
            ft.dimord = 'chan_time';
            
            % find the condition labels used to match up the data
            
            if partition.is_partition
               if partition.partition_number == 1
                   if isfield(data, 'p1_pgi')
                        pgi = data.p1_pgi;
                   end
                   thin = data.p1_thin;
                   med = data.p1_med;
                   thick = data.p1_thick;   
               elseif partition.partition_number == 2
                   if isfield(data, 'p2_pgi')
                        pgi = data.p2_pgi;
                   end
                   thin = data.p2_thin;
                   med = data.p2_med;
                   thick = data.p2_thick;
               elseif partition.partition_number == 3
                   if isfield(data, 'p3_pgi')
                        pgi = data.p3_pgi;
                   end
                   thin = data.p3_thin;
                   med = data.p3_med;
                   thick = data.p3_thick;
               end
            elseif ~partition.is_partition
                pgi = data.med - (data.thin + data.thick)/2;
                thin = data.thin;
                med = data.med;
                thick = data.thick;
                ft.avg = pgi;
            end
            
            if isfield(data, 'p1_pgi') || isfield(data, 'p2_pgi') || isfield(data, 'p3_pgi') 
                ft.avg = pgi;
            end
            
            ft.thin = thin;
            ft.med = med;
            ft.thick = thick;
            
            ft_regression_data{idx_used_for_saving_data} = ft;
            participant_order{idx_used_for_saving_data} = i;
            idx_used_for_saving_data = idx_used_for_saving_data + 1;
        end
    end
end

%% related to bootstrapping the erps
function ci = bootstrap_erps(data, e_idx)
    [~, n_participants] = size(data);
    [all_med, all_thick, all_thin] = deal({}, {}, {});
    
    % get all of the participant trials into one matrix of each type
    for participant=1:n_participants
        p = data{participant};
        fields = fieldnames(p);
        for k=1:numel(fields)
           time_series_name = fields{k}; 
           participant_level.series = p.(fields{k});
           participant_level.weighting = p.weighting;
           
           if contains(time_series_name, 'med')
                participant_level.series = participant_level.series(e_idx,:);
                all_med{end+1} = participant_level;
           elseif contains(time_series_name, 'thick')
                participant_level.series = participant_level.series(e_idx,:);
                all_thick{end+1} = participant_level;
           elseif contains(time_series_name, 'thin')
                participant_level.series = participant_level.series(e_idx,:);
                all_thin{end+1} = participant_level;
           end
           
        end
    end
    
    [~, n_participants] = size(data);
    
    
    % start the bootstrapping process to create the plots with CIs
    n_iterations = 3000;
    [dist_med, dist_thin, dist_thick, dist_pgi] = deal([ ], [], [], []);
    for n =1:n_iterations
        % sample trials with replacement
        sampled_med = datasample(all_med,n_participants);
        sampled_thick = datasample(all_thick,n_participants);
        sampled_thin = datasample(all_thin,n_participants);
        
        % weight the ERPs using the arithmetic mean amd create a
        % bootstrapped ERP
        avg_med = calculate_aritmetic_mean(sampled_med);
        avg_thick = calculate_aritmetic_mean(sampled_thick);
        avg_thin = calculate_aritmetic_mean(sampled_thin);
        avg_pgi = avg_med - (avg_thin + avg_thick)/2;
        
        % add to our distribution of ERPs
        dist_med(:,:,n) = avg_med;
        dist_thin(:,:,n) = avg_thin;
        dist_thick(:,:,n) = avg_thick;
        dist_pgi(:,:,n) = avg_pgi;
    end
    
    [dist_med_low, dist_med_high, dist_med_avg] = deal([], [], []);
    [dist_thick_low, dist_thick_high, dist_thick_avg] = deal([], [], []);
    [dist_thin_low, dist_thin_high, dist_thin_avg] = deal([], [], []);
    [dist_pgi_low, dist_pgi_high, dist_pgi_avg] = deal([], [], []);
    
    
    n_samples = size(dist_med, 2);
    for i=1:n_samples
        
        % medium 2.5% and 95% CI
        med_at_time_t = dist_med(1,i,:);
        lower = prctile(med_at_time_t, 2.5);
        upper = prctile(med_at_time_t, 97.5);
        dist_med_low(i) = lower;
        dist_med_high(i) = upper;
        dist_med_avg(i) = mean(med_at_time_t);
        
        % thick 2.5% and 95% CI
        thick_at_time_t = dist_thick(1,i,:);
        lower = prctile(thick_at_time_t, 2.5);
        upper = prctile(thick_at_time_t, 95);
        dist_thick_low(i) = lower;
        dist_thick_high(i) = upper;
        dist_thick_avg(i) = mean(thick_at_time_t);
        
        % thin 2.5% and 95% CI
        thin_at_time_t = dist_thin(1,i,:);
        lower = prctile(thin_at_time_t, 2.5);
        upper = prctile(thin_at_time_t, 95);
        dist_thin_low(i) = lower;
        dist_thin_high(i) = upper;
        dist_thin_avg(i) = mean(thin_at_time_t);
        
        % pgi 2.5% and 95% CI
        pgi_at_time_t = dist_pgi(1,i,:);
        lower = prctile(pgi_at_time_t, 2.5);
        upper = prctile(pgi_at_time_t, 95);
        dist_pgi_low(i) = lower;
        dist_pgi_high(i) = upper;
        dist_pgi_avg(i) = mean(pgi_at_time_t);
    end
    
    ci.dist_pgi_low = dist_pgi_low;
    ci.dist_pgi_high = dist_pgi_high;
    ci.dist_pgi_avg = dist_pgi_avg;
    
    ci.dist_thin_low = dist_thin_low;
    ci.dist_thin_high = dist_thin_high;
    ci.dist_thin_avg = dist_thin_avg;
    
    ci.dist_med_high = dist_med_high;
    ci.dist_med_low = dist_med_low;
    ci.dist_med_avg = dist_med_avg;
    
    ci.dist_thick_low = dist_thick_low;
    ci.dist_thick_high = dist_thick_high;
    ci.dist_thick_avg = dist_thick_avg;
    
end
