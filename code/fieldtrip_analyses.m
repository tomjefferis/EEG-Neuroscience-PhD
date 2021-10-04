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

%% SHALL WE APPLY A ROI, IF SO HOW?
region_of_interest = 0;
roi_applied = 'two-tailed';
weight_roi = 0;
roi_to_apply = 0;

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
    
%% OFF TO THE RACES WE GO
for i = 1:numel(experiment_types)
    for j = 1:numel(desired_design_mtxs)
        experiment_type = experiment_types{i};
        desired_design_mtx = desired_design_mtxs{j};
        %% create the results save path depending on the experiment
        if strcmp(experiment_type, 'partitions-2-8')
            save_path = strcat(results_dir, '\', 'partitions', '\', desired_design_mtx);
        elseif contains(experiment_type, 'erps-23-45-67')
            save_path = strcat(results_dir, '\', 'onsets', '\', desired_design_mtx);
        elseif contains(experiment_type,'onsets-2-8-explicit')
            save_path = strcat(results_dir, '\', 'mean_intercept', '\', desired_design_mtx);
        elseif strcmp(experiment_type, 'partitions_vs_onsets')
            save_path = strcat(results_dir, '\', 'partitions_vs_onsets', '\', desired_design_mtx);
        end
        
        %% Are we looking at onsets 2-8 or partitions
        % set up the experiment as needed
        if strcmp(experiment_type, 'onsets-2-8-explicit')
            data_file = 'time_domain_mean_intercept_onsets_2_3_4_5_6_7_8_grand-average_baseline_3.mat';
            regressor = 'ft_statfun_depsamplesT';
            n_participants = 40;
            start_latency = -0.2;
            end_latency = 3.9;

            partition.is_partition = 0;
            partition.partition_number = 0;

            [data, participant_order_1] = load_postprocessed_data(main_path, n_participants, ...
                data_file, partition);
            n_part = numel(data);
            design_matrix =  [1:n_part 1:n_part; ones(1,n_part) 2*ones(1,n_part)]; 


        elseif strcmp(experiment_type, 'partitions-2-8')
            partition1.is_partition = 1; % partition 1
            partition1.partition_number = 1;
            partition2.is_partition = 1; % partition 2
            partition2.partition_number = 2;
            partition3.is_partition = 1; % partition 3
            partition3.partition_number = 3;
            type_of_effect = 'habituation';
            regressor = 'ft_statfun_indepsamplesregrT';
            regression_type = desired_design_mtx;
            n_participants = 39;
            
            if strcmp(type_of_analysis,'time_domain')
                data_file = 'partitions_partitioned_onsets_2_3_4_5_6_7_8_grand-average.mat';
                
                [data1, participant_order_1] = load_postprocessed_data(main_path, n_participants, ...
                    data_file, partition1);
                [data2, participant_order_2] = load_postprocessed_data(main_path, n_participants, ...
                    data_file, partition2);
                [data3, participant_order_3] = load_postprocessed_data(main_path, n_participants, ...
                    data_file, partition3);
                
                partition = 1;
                [design1, new_participants1] = create_design_matrix_partitions(participant_order_1, data1, ...
                        regression_type, partition, type_of_effect);
                partition = 2;
                [design2, new_participants2] = create_design_matrix_partitions(participant_order_2, data2, ...
                        regression_type, partition, type_of_effect);
                partition = 3;
                [design3, new_participants3] = create_design_matrix_partitions(participant_order_3, data3, ...
                        regression_type, partition, type_of_effect);
                
                data = [new_participants1, new_participants2, new_participants3];
                design_matrix = [design1, design2, design3];
                n_part = numel(data);
                n_part_per_desgin = numel(design1);

                if strcmp(desired_design_mtx, 'no-factor')
                    design1(1:numel(design1)) = 2.72;
                    design2(1:numel(design1)) = 1.65;
                    design3(1:numel(design1)) = 1.00;
                    design_matrix = [design1, design2, design3];
                end

                design_matrix = design_matrix - mean(design_matrix);
                save_desgin_matrix(design_matrix, n_part_per_desgin, save_path, 'habituation')

                if region_of_interest == 1
                    if strcmp(experiment_type, 'partitions-2-8') 
                        if strcmp(roi_applied, 'one-tailed')
                            load('D:\PhD\fieldtrip\roi\one_tailed_roi_28.mat');
                        elseif strcmp(roi_applied, 'two-tailed')
                            load('D:\PhD\fieldtrip\roi\two_tailed_roi_28.mat');
                        end
                    end
                    data = create_hacked_roi(data, roi, weight_roi);
                end
                
            elseif strcmp(type_of_analysis, 'frequency_domain')
                if compute_frequency_data == 1
                    analysis = 'preprocess';
                    data_file = 'frequency_domain_partitions_partitioned_onsets_2_3_4_5_6_7_8_trial-level.mat';
                else
                    analysis = 'load';
                    data_file = 'frequency_domain_partitions_partitioned_onsets_2_3_4_5_6_7_8_grand-average.mat';
                end
                
                electrode = return_mua_electrode(desired_design_mtx);
                
                if strcmp(frequency_type, 'fourier')
                    fnames = {'3_fourier_med.mat','3_fourier_thin.mat','3_fourier_thick.mat'};
                elseif strcmp(frequency_type, 'pow')
                    fnames = {'pow_med.mat','pow_thin.mat','pow_thick.mat'};
                end
                
                [data1, participant_order1] = load_postprocessed_data(main_path, n_participants, ...
                    data_file, partition1);
                [data2, participant_order2] = load_postprocessed_data(main_path, n_participants, ...
                    data_file, partition2);
                [data3, participant_order3] = load_postprocessed_data(main_path, n_participants, ...
                    data_file, partition3);          
                
                p1_freq = to_frequency_data(data1, main_path, 1, ...
                    participant_order1, analysis,fnames,frequency_type, electrode);                
                p2_freq = to_frequency_data(data2, main_path, 2, ...
                    participant_order2, analysis,fnames,frequency_type, electrode);
                p3_freq = to_frequency_data(data3, main_path, 3, ...
                    participant_order3, analysis,fnames,frequency_type, electrode);
                                
                if strcmp(desired_design_mtx, 'no-factor') && compute_frequency_data == 0
                   [~] = compute_spectrogram(p1_freq,1,save_path, electrode, desired_design_mtx, 'all', frequency_type);
                   [~] = compute_spectrogram(p2_freq,2,save_path, electrode, desired_design_mtx, 'all', frequency_type);
                   [~] = compute_spectrogram(p3_freq,3,save_path, electrode, desired_design_mtx, 'all', frequency_type);
                elseif ~strcmp(desired_design_mtx, 'no-factor') && compute_frequency_data == 0
                    
                    % get the partitioned data
                    [p1_freq_h, p1_freq_l, ~, ~] = get_partitions_medium_split(p1_freq, participant_order1,...
                        desired_design_mtx, 1, type_of_effect, 0, 0);
                    f_data.p1_freq_h = compute_spectrogram(p1_freq_h,1,save_path, electrode, desired_design_mtx, 'high', frequency_type);
                    f_data.p1_freq_l = compute_spectrogram(p1_freq_l,1,save_path, electrode, desired_design_mtx, 'low', frequency_type);
                    plot_spectrogram(f_data.p1_freq_h,save_path,1,frequency_type,electrode, 'High')
                    plot_spectrogram(f_data.p1_freq_l,save_path,1,frequency_type,electrode, 'Low')
                    
                    [p2_freq_h, p2_freq_l, ~, ~] = get_partitions_medium_split(p2_freq, participant_order2,...
                        desired_design_mtx, 1, type_of_effect, 0, 0);  
                    f_data.p2_freq_h = compute_spectrogram(p2_freq_h,2,save_path, electrode, desired_design_mtx, 'high', frequency_type);
                    f_data.p2_freq_l = compute_spectrogram(p2_freq_l,2,save_path, electrode, desired_design_mtx, 'low', frequency_type);
                    plot_spectrogram(f_data.p2_freq_h,save_path,2,frequency_type,electrode, 'High')
                    plot_spectrogram(f_data.p2_freq_l,save_path,2,frequency_type,electrode, 'Low')
                    
                    [p3_freq_h, p3_freq_l, ~, ~] = get_partitions_medium_split(p3_freq, participant_order3,...
                        desired_design_mtx, 1, type_of_effect, 0, 0);
                    f_data.p3_freq_h = compute_spectrogram(p3_freq_h,3,save_path, electrode, desired_design_mtx, 'high', frequency_type);
                    f_data.p3_freq_l = compute_spectrogram(p3_freq_l,3,save_path, electrode, desired_design_mtx, 'low', frequency_type);
                    plot_spectrogram(f_data.p3_freq_h,save_path,3,frequency_type,electrode, 'High')
                    plot_spectrogram(f_data.p3_freq_l,save_path,3,frequency_type,electrode, 'Low')
                    
                    if analysis_on_aggr_data == 1          
                        % compute the aggregate freq-pow
                        aggr_data = aggregate_freq_data(f_data, frequency_type);

                        % create the grand average plot for each stimulus type
                        plot_spectrogram(aggr_data,save_path,123,frequency_type, electrode, 'All')

                        % based on the grand average, apply an ROI to each
                        % participant and extract a value, this is manually
                        % written in the design matrix fn.
                        time = [0.090, 0.175];
                        freq = [5, 8];
                        average_power_values(p1_freq, freq, time, electrode);
                        average_power_values(p2_freq, freq, time, electrode);
                        average_power_values(p3_freq, freq, time, electrode);
                    end
                end
                
                % let arnold schwarzenegger tell me when the analysis is
                % complete
                continue;
            end
                
        elseif strcmp(experiment_type, 'erps-23-45-67') 
            data_file23 = 'mean_intercept_onsets_2_3_grand-average.mat';
            data_file45 = 'mean_intercept_onsets_4_5_grand-average.mat';
            data_file67 = 'mean_intercept_onsets_6_7_grand-average.mat';
            type_of_effect = 'sensitization';
            regressor = 'ft_statfun_indepsamplesregrT';

            n_participants = 40;
            partitions.is_partition = 0;
            partitions.partition_number = 0;
            regression_type = desired_design_mtx;

            [data1, participant_order_1] = load_postprocessed_data(main_path, n_participants, ...
                data_file23, partitions);
            [data2, participant_order_2] = load_postprocessed_data(main_path, n_participants, ...
                data_file45, partitions);
            [data3, participant_order_3] = load_postprocessed_data(main_path, n_participants, ...
                data_file67, partitions);


            partition = 1;
            [design1, new_participants1] = create_design_matrix_partitions(participant_order_1, data1, ...
                    regression_type, partition, type_of_effect);
            partition = 2;
            [design2, new_participants2] = create_design_matrix_partitions(participant_order_2, data2, ...
                    regression_type, partition, type_of_effect);
            partition = 3;
            [design3, new_participants3] = create_design_matrix_partitions(participant_order_3, data3, ...
                    regression_type, partition, type_of_effect);

            data = [new_participants1, new_participants2, new_participants3];
            design_matrix = [design1, design2, design3];
            n_part = numel(data);
            n_part_per_desgin = numel(design1);


            if strcmp(desired_design_mtx, 'no-factor')
                design_matrix_a(1:numel(design1)) = 1.00;
                design_matrix_b(1:numel(design1)) = 1.65;
                design_matrix_c(1:numel(design1)) = 2.72;
                design_matrix = [design_matrix_a, design_matrix_b, design_matrix_c];
                min_x = min(design_matrix);
                design_matrix = design_matrix - min_x;
            end

            if strcmp(desired_design_mtx, 'headache-scores-entire-median-split') ... 
                    || strcmp(desired_design_mtx, 'discomfort-scores-entire-median-split')
                newdesign1 = design1 * 1.00;
                newdesign2 = design2 * 1.65;
                newdesign3 = design3 * 2.72;
                design_matrix = [newdesign1, newdesign2, newdesign3];
                min_x = min(design_matrix);
                design_matrix = design_matrix - min_x;  
            end

            design_matrix = design_matrix - mean(design_matrix);
            save_desgin_matrix(design_matrix, n_part_per_desgin, save_path, 'sensitization')

            if region_of_interest == 1
                if strcmp(experiment_type, 'erps-23-45-67') || strcmp(experiment_type,  'erps-23-45-67-no-factor')
                    if strcmp(roi_applied, 'one-tailed')
                        load('D:\PhD\fieldtrip\roi\one_tailed_roi_28.mat');
                    elseif strcmp(roi_applied, 'two-tailed')
                        load('D:\PhD\fieldtrip\roi\two_tailed_roi_28.mat');
                    end
                end
                data = create_hacked_roi(data, roi, weight_roi);
            end
        elseif strcmp(experiment_type, 'partitions_vs_onsets')
            onsets_2_3 = 'time_domain_partitions_partitioned_onsets_2_3_grand-average.mat';
            onsets_4_5 = 'time_domain_partitions_partitioned_onsets_4_5_grand-average.mat';
            onsets_6_7 = 'time_domain_partitions_partitioned_onsets_6_7_grand-average.mat';
            type_of_effect = 'habituation';
            regressor = 'ft_statfun_indepsamplesregrT';
            
            regression_type = desired_design_mtx;
            n_participants = 39;
            partition1.is_partition = 1; 
            partition1.partition_number = 1;
            partition2.is_partition = 1; 
            partition2.partition_number = 2;
            partition3.is_partition = 1; 
            partition3.partition_number = 3;
            
            [p1_23, po_p1_23] = load_postprocessed_data(main_path, n_participants, ...
                onsets_2_3, partition1);
            [p2_23, po_p2_23] = load_postprocessed_data(main_path, n_participants, ...
                onsets_2_3, partition2);
            [p3_23, po_p3_23] = load_postprocessed_data(main_path, n_participants, ...
                onsets_2_3, partition3);
                                    
            partition = 1;
            [design_p1_23, p1_23] = create_design_matrix_partitions(po_p1_23, p1_23, ...
                    regression_type, partition, type_of_effect);
            partition = 2;
            [design_p2_23, p2_23] = create_design_matrix_partitions(po_p2_23, p2_23, ...
                    regression_type, partition, type_of_effect);
            partition = 3;
            [design_p3_23, p3_23] = create_design_matrix_partitions(po_p3_23, p3_23, ...
                    regression_type, partition, type_of_effect);
                
            design_p1_23 = design_p1_23 * 1;
            design_p2_23 = design_p2_23 * 1;
            design_p3_23 = design_p3_23 * 1;
            
            [p1_45, po_p1_45] = load_postprocessed_data(main_path, n_participants, ...
                onsets_4_5, partition1);
            [p2_45, po_p2_45] = load_postprocessed_data(main_path, n_participants, ...
                onsets_4_5, partition2);
            [p3_45, po_p3_45] = load_postprocessed_data(main_path, n_participants, ...
                onsets_4_5, partition3);
            
            partition = 1;
            [design_p1_45, p1_45] = create_design_matrix_partitions(po_p1_45, p1_45, ...
                    regression_type, partition, type_of_effect);
            partition = 2;
            [design_p2_45, p2_45] = create_design_matrix_partitions(po_p2_45, p2_45, ...
                    regression_type, partition, type_of_effect);
            partition = 3;
            [design_p3_45, p3_45] = create_design_matrix_partitions(po_p3_45, p3_45, ...
                    regression_type, partition, type_of_effect);
                
            design_p1_45 = design_p1_45 * 1.65;
            design_p2_45 = design_p2_45 * 1.65;
            design_p3_45 = design_p3_45 * 1.65;
            
            [p1_67, po_p1_67] = load_postprocessed_data(main_path, n_participants, ...
                onsets_2_3, partition1);
            [p2_67, po_p2_67] = load_postprocessed_data(main_path, n_participants, ...
                onsets_2_3, partition2);
            [p3_67, po_p3_67] = load_postprocessed_data(main_path, n_participants, ...
                onsets_2_3, partition3);
            
            partition = 1;
            [design_p1_67, p1_67] = create_design_matrix_partitions(po_p1_67, p1_67, ...
                    regression_type, partition, type_of_effect);
            partition = 2;
            [design_p2_67, p2_67] = create_design_matrix_partitions(po_p2_67, p2_67, ...
                    regression_type, partition, type_of_effect);
            partition = 3;
            [design_p3_67, p3_67] = create_design_matrix_partitions(po_p3_67, p3_67, ...
                    regression_type, partition, type_of_effect);
            
            design_p1_67 = design_p1_67 * 2.72;
            design_p2_67 = design_p2_67 * 2.72;
            design_p3_67 = design_p3_67 * 2.72;
            
            design_matrix = [
                design_p1_23, design_p1_45, design_p1_67, ...
                design_p2_23, design_p2_45, design_p2_67, ...
                design_p3_23, design_p3_45, design_p3_67, ...
            ];
        
            design_matrix = design_matrix - mean(design_matrix);
            
            
            first_partition_size = size([design_p1_23, design_p1_45, design_p1_67],2);
            plot(design_matrix(1:first_partition_size));
            hold;
            second_partition_size = size([design_p2_23, design_p2_45, design_p2_67],2);
            plot(design_matrix(first_partition_size+1:(second_partition_size*2)));
            third_partition_size = size([design_p3_23, design_p3_45, design_p3_67],2);
            plot(design_matrix((second_partition_size*2)+1:third_partition_size*3));
            xlabel('Participants');
            ylabel('3-way Interaction');
            legend({'P1 (Onsets 2:3, 4:5, 6:7)', 'P2 (Onsets 2:3, 4:5, 6:7)', 'P3 (Onsets 2:3, 4:5, 6:7)'},'Location','northwest')
            set(gcf,'Position',[100 100 1000 1000])
            save_dir = strcat(save_path, '\', 'design_matrix.png');
            exportgraphics(gcf,save_dir,'Resolution',500);
            close; 
               
            data = [
                p1_23, p1_45, p1_67, ...
                p2_23, p2_45, p2_67, ...
                p3_23, p3_45, p3_67
            ];
        
            compute_peak
               
        end
        %% setup FT analysis
        % we have to switch to SPM8 to use some of the functions in FT
        addpath 'W:\PhD\MatlabPlugins\spm8';
        rmpath 'W:\PhD\MatlabPlugins\spm12';

        % we need to tell fieldtrip how our electrodes are structured
        cfg = [];
        cfg.feedback = 'no';
        cfg.method = 'distance';
        cfg.elec = data{1}.elec;
        neighbours = ft_prepare_neighbours(cfg);

        % all experiment configuraitons:
        cfg = [];
        cfg.latency = [start_latency, end_latency];
        cfg.channel = 'eeg';
        cfg.statistic = regressor;
        cfg.method = 'montecarlo';
        cfg.correctm = 'cluster';
        cfg.neighbours = neighbours;
        cfg.clusteralpha = 0.025;
        cfg.numrandomization = 1;
        cfg.tail = roi_to_apply; 
        cfg.design = design_matrix;
        cfg.computeprob = 'yes';
        cfg.alpha = 0.05; 
        cfg.correcttail = 'alpha'; 
        
        
        %% run the fieldtrip analyses
        if contains(experiment_type, 'onsets-2-8-explicit') || contains(experiment_type, 'onsets-1-explicit')
            cfg.uvar = 1;
            cfg.ivar = 2;
            null_data = set_values_to_zero(data); % create null data to hack a t-test
            stat = ft_timelockstatistics(cfg, data{:},null_data{:});
            desired_cluster =1;
            %get_region_of_interest_electrodes(stat, desired_cluster, experiment_type, roi_applied);
        elseif contains(experiment_type, 'partitions') || contains(experiment_type, 'onsets-2-8-factor') ...
                || contains(experiment_type, 'onsets-1-factor') || contains(experiment_type, 'erps-23-45-67') ...
                || contains(experiment_type, 'coarse-vs-fine-granularity') || contains(experiment_type, 'Partitions')
            cfg.ivar = 1;
            stat = ft_timelockstatistics(cfg, data{:});
        end

        %% get peak level stats
        [pos_peak_level_stats, pos_all_stats] = get_peak_level_stats(stat, 1, 'positive');
        [neg_peak_level_stats, neg_all_stats] = get_peak_level_stats(stat, 1, 'negative');

        %% function that plots the t values through time and decides whcih electrode to plot
        if numel(pos_all_stats) > 0
            pos_peak_level_stats = compute_best_electrode_from_t_values(stat,pos_all_stats,save_path, 'positive', pos_peak_level_stats);
        end
        if numel(neg_all_stats) > 0
            neg_peak_level_stats = compute_best_electrode_from_t_values(stat,neg_all_stats,save_path, 'negative', neg_peak_level_stats);
        end
        
        %% generate ERPs using the stat information
        if generate_erps == 1
            generate_peak_erps(master_dir, main_path, experiment_type, ...
                stat, pos_peak_level_stats, 'positive', desired_design_mtx, 1, ...
                save_path, weight_erps, weighting_factor);
            generate_peak_erps(master_dir, main_path, experiment_type, ...
                stat, neg_peak_level_stats, 'negative', desired_design_mtx, 1, ...
                save_path, weight_erps, weighting_factor);
        end
        
        %% get cluster level percentage through time
        % 1 for the most positive going cluster
        xlim = 256;
        title = 'Most positive going cluster through time as a % of entire volume';
        calculate_cluster_size(stat, 1, title, xlim, 'positive', ...
            save_path);
        title = 'Most negative going cluster through time as a % of entire volume';
        calculate_cluster_size(stat, 1, title, xlim, 'negative', ...
            save_path);
        
        %% make pretty plots
        create_viz_topographic_maps(data, stat, start_latency, end_latency, ...
            0.05, 'positive', save_path)
        create_viz_topographic_maps(data, stat, start_latency, end_latency, ...
            0.05, 'negative', save_path)
    end
end

%% plot the t values through time and select the best electrode
function peak_stat_info = compute_best_electrode_from_t_values(stat, electrode_stats, save_dir, tail, peak_stat_info)
    significant_masks = stat.mask;
    electrodes = electrode_stats.electrodes;
    
    if strcmp(tail, 'positive')
        cluster_labels = stat.posclusterslabelmat;
        fname = '\positive_t_values_through_time.png';
    else
        cluster_labels = stat.negclusterslabelmat;
        fname = '\negative_t_values_through_time.png';
    end
    
    time = stat.time;
    most_significant_electrode = electrode_stats.electrodes{1};
    previous_best = 0;
    
    for i=1:numel(electrodes)
        electrode_idx = find(strcmp(stat.label,electrodes(i)));
        significant_time_series = significant_masks(electrode_idx,:);
        cluster_time_series = cluster_labels(electrode_idx,:);
        raw_time_series = zeros(1,numel(cluster_time_series));
        significance_count = 0;
        for j=1:numel(significant_time_series)
            if cluster_time_series(j)==1
                is_significant = significant_time_series(j);
                significance_count = significance_count + is_significant;
                raw_time_series(j) = significance_count;
            else
                raw_time_series(j) = significance_count;
            end
        end
        if previous_best < significance_count && strcmp(tail, 'positive')
            e = electrodes(i);
            most_significant_electrode = e;
            previous_best = significance_count;
        elseif previous_best > significance_count && strcmp(tail, 'negative')
            e = electrodes(i);
            most_significant_electrode = e;
            previous_best = significance_count;
        end
        if sum(raw_time_series) > 0
            plot(time, raw_time_series);
            hold on;
        end
    end
    
   [col, ~] = size(electrode_stats);
    
   for i = 1:col
      if strcmp(electrode_stats.electrodes{i}, most_significant_electrode)
          peak_stat_info.electrode = most_significant_electrode;
          peak_stat_info.time = electrode_stats.time(i);
          peak_stat_info.t_value = electrode_stats.t_value(i);
      end
   end
       
    
    xlabel('Time (ms)');
    ylabel('Cumulative Frequency (If T-value is significant');
    legend(electrodes,'Location','northwest')
    title(strcat('Sustained Significant Electrode:',most_significant_electrode)); 
    save_dir = strcat(save_dir, fname);
    set(gcf,'Position',[100 100 1000 1000])
    exportgraphics(gcf,save_dir,'Resolution',500);
    close;
end

%% plot and save the design matrix
function save_desgin_matrix(design_matrix, n_participants, save_path, experiment_type)
    plot(design_matrix(1:n_participants), 'color', 'r', 'LineWidth', 1);
    hold on;
    plot(design_matrix(n_participants+1:n_participants*2), 'color', 'g', 'LineWidth', 1);
    plot(design_matrix((n_participants*2)+1:n_participants*3), 'color', 'b', 'LineWidth', 1);
    xlabel('Participants');
    ylabel('Interaction');
    if strcmp(experiment_type, 'habituation')
        legend({'P1', 'P2', 'P3'},'Location','northwest')
    else
        legend({'Onsets 2:3', 'Onsets 4:5', 'Onsets 6:7'},'Location','northwest')
    end
    set(gcf,'Position',[100 100 1000 1000])
    save_dir = strcat(save_path, '\', 'design_matrix.png');
    exportgraphics(gcf,save_dir,'Resolution',500);
    close;   
end

%% generate ERPs
function generate_peak_erps(master_dir, main_path, experiment_type, ...
    stat, peak_information, effect_type, regression_type, desired_cluster, ...
    save_dir, weight_erps, weighting_factor)
    
    df = stat.df;
    time = stat.time;
    peak_electrode = peak_information.electrode;
    peak_time = peak_information.time;
    peak_t_value = peak_information.t_value;
    
    if strcmp(effect_type, 'positive')
        if numel(stat.posclusters) < 1 || strcmp(peak_electrode, '')
           return ; 
        end
        
        labels = stat.posclusterslabelmat;
        labels(labels>desired_cluster) =0;
        pvalue = round(stat.posclusters(desired_cluster).prob,4);
        cluster_size = stat.posclusters(desired_cluster).clusterstat;
        plot_desc = 'positive_peak_erp.png';
    else
        if numel(stat.negclusters) < 1 || strcmp(peak_electrode, '')
            return ; 
        end
        
        
        labels = stat.negclusterslabelmat;
        peak_t_value = abs(peak_t_value);
        labels(labels>desired_cluster) =0;
        pvalue = round(stat.negclusters(desired_cluster).prob,4);
        cluster_size = stat.negclusters(desired_cluster).clusterstat;
        plot_desc = 'negative_peak_erp.png';
    end

    through_time = sum(labels);
    start_idx = find(through_time(:), desired_cluster, 'first');
    end_idx = find(through_time(:), desired_cluster, 'last');
    start_of_effect = time(start_idx);
    end_of_effect = time(end_idx);
    
    save_dir = strcat(save_dir, '\', plot_desc);
    
    generate_plots(master_dir, main_path, experiment_type, start_of_effect,...
        end_of_effect, peak_electrode, peak_time, peak_t_value, df, ...
        regression_type, pvalue, cluster_size, save_dir, effect_type, weight_erps, weighting_factor)
    
    close;
end

%% Create a hacked ROI in the data
function new_data = create_hacked_roi(data, roi, weight_roi)
    roi_clusters = roi.clusters;
    roi_time = roi.time;
    
    if weight_roi == 1
        [electrodes, ~] = size(roi_clusters);
        total_sum = sum(roi_clusters);
        weighted_template = total_sum/electrodes;
    end

    start_latency = NaN;
    end_latency = NaN;
    new_data = {};
    for idx_i = 1:numel(data)
        each_participant = data{idx_i};
        participant_data = each_participant.avg;
        participant_time = each_participant.time;
        [electrodes, time] = size(participant_data);
        new_participant_data = NaN(electrodes, time);
        time_x = 1:1:time;
        time_x = time_x * 2;
     
        for roi_idx = 1:numel(roi_time)
            t = roi_time(roi_idx)*1000;
            [~,idx]=min(abs(time_x-t));
            idx = idx + 100; % for baselining period
            clusters_at_t = roi_clusters(:,roi_idx);
            
            if isnan(start_latency) && sum(clusters_at_t)>=1
                start_latency =t;
            elseif roi_idx == numel(roi_time)
                end_latency = t; 
            end
            
            for electrode_idx=1:numel(clusters_at_t)
                if clusters_at_t(electrode_idx) == 1
                    new_participant_data(electrode_idx, idx) = participant_data(electrode_idx, idx);
                    if weight_roi == 1
                        weighting_at_t =  weighted_template(roi_idx);
                        new_participant_data(electrode_idx, idx) = new_participant_data(electrode_idx, idx) * weighting_at_t;
                    end      
                end
            end
            
        end
        

        
        
        %[~,start_idx]=min(abs(time_x-start_latency));
        %[~,end_idx]=min(abs(time_x-end_latency));
        %new_participant_data = new_participant_data(:,start_idx:end_idx);
        each_participant.avg = new_participant_data;  
        
        %time = each_participant.time;
        %[~,start_idx]=min(abs(time-start_latency/1000));
        %[~,end_idx]=min(abs(time-end_latency/1000));
        
        new_data{idx_i} = each_participant;
    end
end

%% baed on onsets 2-8 mean intercept effect, get the start and end latency
%% based on our ROI
function get_start_and_end_latency_of_mean_intercept(stat)
    load(roi_fn);
    roi_labels = roi.label;
    stat_labels = stat.label;
    time = stat.time;
    cluster_labels = stat.posclusterslabelmat;
    cluster_labels(cluster_labels>1) = 0; 
    cluster_sum_by_time = sum(cluster_labels);
    
    start_found = 0;
    end_found=0;
    start_latecy = 0;
    end_latency = 0;
    
    for i = 1:numel(cluster_sum_by_time)
       if cluster_sum_by_time(i) > 0
          electrodes =  cluster_labels(:,i);
          electrodes_with_effect = find(electrodes==1);
          t = time(i);
            
          for electrode_idx = 1:numel(stat_labels)
              if any(electrodes_with_effect(:)==electrode_idx)
                electrode_name = stat_labels(electrode_idx);
                if numel(intersect(electrode_name, roi_labels)) ==1
                   if start_found == 0
                      start_latency =  t;
                      start_found =1;
                   end          
                end
              end
          end
       end
    end
end

%% based on onsets 2-8 mean intercept effect, save the electrodes 
% which are at the peak level - used for ROI analysis
% saves the neighbourhood of electrodes for use later on
function get_region_of_interest_electrodes(stat, desired_cluster, experiment_type, roi_applied)
    roi_clusters = stat.posclusterslabelmat;
    roi_clusters(roi_clusters>desired_cluster) = 0;
    time = stat.time;
    roi.clusters = roi_clusters;
    roi.time = time;
    
    if contains(experiment_type, 'onsets-2-8')  && contains(roi_applied, 'one-tailed')
        save 'C:\ProgramFiles\PhD\fieldtrip\roi\one_tailed_roi_28.mat' roi; 
    elseif contains(experiment_type, 'onsets-1') && contains(roi_applied, 'one-tailed')
        save 'C:\ProgramFiles\PhD\fieldtrip\roi\one_tailed_roi_1.mat' roi; 
    elseif contains(experiment_type, 'onsets-2-8')  && contains(roi_applied, 'two-tailed')
        save 'C:\ProgramFiles\PhD\fieldtrip\roi\two_tailed_roi_28.mat' roi; 
    elseif contains(experiment_type, 'onsets-1') && contains(roi_applied, 'two-tailed')
        save 'C:\ProgramFiles\PhD\fieldtrip\roi\two_tailed_roi_1.mat' roi; 
    end
   
end

    %% plots the cluster effect through time
function create_viz_cluster_effect(stat, alpha)
    cfg = [];
    cfg.alpha = alpha;
    cfg.parameter = 'stat';
    cfg.zlim = [-4 4];
    cfg.xlim = [0.56, 0.256];
    cfg.highlightcolorpos = [1,0,0]; 
    cfg.highlightsymbolseries =  ['*', '*','*','*','*'];
    %cfg.toi = [0.56, 0.256];
    ft_clusterplot(cfg, stat);
end

%% used to calculate the cluster size through time
function calculate_cluster_size(stat, desired_cluster, ptitle, xlim_t, type, save_dir)
    if contains(type, 'positive')
        cluster_labelling_mtx = stat.posclusterslabelmat;
        save_dir = strcat(save_dir, '\', 'positive_cluster.png');
    elseif contains(type, 'negative')
        cluster_labelling_mtx = stat.negclusterslabelmat;
        save_dir = strcat(save_dir, '\', 'negative_cluster.png');
    end
    
    time_mtx = stat.time;
    
    [electrodes, time] = size(cluster_labelling_mtx);
    
    cluster_size_through_time = zeros(2,time);
    
    for i = 1:time
        t = time_mtx(i);
        electrodes_in_time = cluster_labelling_mtx(:,i);
        clusters_in_time = find(electrodes_in_time==desired_cluster);
        cluster_size_through_time(1,i) = t;
        cluster_size_through_time(2,i) = numel(clusters_in_time)/electrodes;
    end
    
    area(cluster_size_through_time(1,:)*1000, cluster_size_through_time(2,:)*100);
    ylim([0 100]);
    grid on;
    xlabel('Time (ms)');
    ylabel('Percentage of cluster');
    %xlim([0,260])
    xlim([0,xlim_t])
    title(ptitle, 'FontSize', 14); 
   
    set(gcf,'Position',[100 100 1000 1000])
    exportgraphics(gcf,save_dir,'Resolution',500);
    close;
end

%% used to generate 'peak' level stats
function [peak_level_electrode, all_electrode_stats] = get_peak_level_stats(stat, cluster, type)
    peak_level_stats = containers.Map;

    time = stat.time;
    electrodes = stat.label;
    t_values = stat.stat;
    if contains(type, 'positive')
        cluster_matrix_locations = stat.posclusterslabelmat;
    elseif contains(type, 'negative')
        cluster_matrix_locations = stat.negclusterslabelmat;
    end
    
    [rows, columns] = size(cluster_matrix_locations);
    
    for col = 1:columns
        for row = 1:rows
            
            t_value = t_values(row, col);
            t = time(col);
            electrode = electrodes{row};

            which_cluster = cluster_matrix_locations(row, col);
            if which_cluster == cluster
                keys = peak_level_stats.keys;
                if any(strcmp(keys,electrode))
                    previous_values = peak_level_stats(electrode);
                    previous_t = previous_values{2};
                    if t_value > previous_t && contains(type, 'positive')
                        peak_level_stats(electrode) = {t, t_value};
                    elseif t_value < previous_t && contains(type, 'negative')
                        peak_level_stats(electrode) = {t, t_value};
                    end
                else
                    peak_level_stats(electrode) = {0, 0};
                end
            end
        end 
    end
    peak_level_electrode = get_most_significant_electrode(peak_level_stats, type);
    all_electrode_stats = get_all_electrode_stats(peak_level_stats, type);
end

function mtx = get_all_electrode_stats(peak_level_stats, type)
    keys = peak_level_stats.keys;
    electrodes = cell(numel(keys), 1);
    time = zeros(numel(keys), 1);
    t_value = zeros(numel(keys), 1);
    
    for i = 1:numel(keys)
        electrode = keys{i};
        stats = peak_level_stats(electrode);
        tv = stats{2};
        t = stats{1};
        
        electrodes{i} = electrode;
        time(i) = t;
        t_value(i) = tv;
    end
    
    mtx = table(electrodes, time, t_value);
    if strcmp(type, 'positive')
        mtx = sortrows(mtx, {'t_value'}, {'descend'});
    elseif strcmp(type, 'negative')
        mtx = sortrows(mtx, {'t_value'}, {'ascend'});
    end
end

function most_significant = get_most_significant_electrode(peak_level_stats, type)
    keys = peak_level_stats.keys;
    
    most_significant.electrode = '';
    most_significant.time = 0;
    most_significant.t_value = 0;
    for i = 1:numel(keys)
        electrode = keys{i};
        stats = peak_level_stats(electrode);
        t_value = stats{2};
        time = stats{1};
        if t_value > most_significant.t_value && contains(type, 'positive')
            most_significant.electrode = electrode;
            most_significant.time = time;
            most_significant.t_value = t_value;
        elseif t_value < most_significant.t_value && contains(type, 'negative')
            most_significant.electrode = electrode;
            most_significant.time = time;
            most_significant.t_value = t_value; 
        end
    end
end
%% used to create the topographic maps
function create_viz_topographic_maps(data1, stat, start_latency, ... 
    end_latency, alpha, type, save_dir)

    if numel(stat.negclusters) < 1 && strcmp(type, 'negative')
        return;
    elseif numel(stat.posclusters) < 1 && strcmp(type,'positive')
        return;
    end

    % fixes the interpolation issue of NaNs
    stat.stat(isnan(stat.stat))=0;
        
    
    cfg = [];
    cfg.channel   = 'all';
    cfg.latency   = 'all';
    cfg.parameter = 'avg';
    
    grand_avg1 = ft_timelockgrandaverage(cfg, data1{:});
    grand_avg1.elec = data1{1}.elec;
    [i1,i2] = match_str(grand_avg1.label, stat.label);
    
    figure;
    timestep = 0.01; % in seconds usually 0.01
    sampling_rate = 512;
    sample_count  = length(stat.time);
    
    j = [start_latency:timestep:end_latency]; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
    m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in M/EEG samples
    
    if contains(type, 'positive')
        pos_cluster_pvals = [stat.posclusters(:).prob];
        pos_clust = find(pos_cluster_pvals < alpha);
        clust = ismember(stat.posclusterslabelmat, pos_clust);
        save_dir = strcat(save_dir, '\', 'positive_topographic.png');
    elseif contains(type, 'negative')
        neg_cluster_pvals = [stat.negclusters(:).prob];
        neg_clust = find(neg_cluster_pvals < alpha);
        clust = ismember(stat.negclusterslabelmat, neg_clust);    
        save_dir = strcat(save_dir, '\', 'negative_topographic.png');
    end
    max_iter = numel(m)-1;
    
    for k = 1:max_iter
         subplot(4,5,k);
         cfg = [];
         cfg.xlim = [j(k) j(k+1)];
         cfg.zlim = [-5e-14 5e-14];
         pos_int = zeros(numel(grand_avg1.label),1);
         pos_int(i1) = all(clust(i2, m(k):m(k+1)), 2);
         cfg.highlight = 'on';
         cfg.highlightchannel = find(pos_int);
         cfg.highlightcolor = {'r',[0 0 1]};
         cfg.comment = 'xlim';
         cfg.commentpos = 'title';
         cfg.parameter = 'stat';
         cfg.zlim=[-4,4];
         %cfg.colorbar = 'SouthOutside';
         %cfg.layout = 'biosemi128.lay';
         ft_topoplotER(cfg, stat);
         %ft_clusterplot(cfg, stat)
    end
    
    %%cfg = [];
    %%cfg.alpha  = 0.025;
    %%cfg.parameter = 'stat';
    %%cfg.zlim   = [-4 4];
    %%ft_clusterplot(cfg, stat);
    set(gcf,'Position',[100 100 1000 1000])
    exportgraphics(gcf,save_dir,'Resolution',500);
    close;
end

%% set all my values to 0 to hack a T test
function data = set_values_to_zero(data)
    for idx = 1:numel(data)
        participant_data = data{idx};
        series = participant_data.avg;
        series(:) = 0;
        participant_data.avg = series;
        data{idx} = participant_data;
    end
end
%% function to create the design matrix for all experiments
function [design, new_participants] = create_design_matrix_partitions(participant_order, participants, ...
    regression_type, partition, type_of_effect)

    scores = return_scores(regression_type, type_of_effect);
    
    if strcmp(regression_type, 'headache-scores-entire-median-split') ...
            || strcmp(regression_type, 'discomfort-scores-entire-median-split') ...
            || strcmp(regression_type, 'vs-scores-entire-median-split')
        score = scores.one;
        calc(:,1) = score(:,2);
        calc(:,2) = score(:,1);
        calc = sortrows(calc);
        [rows, ~] = size(calc);
        calc(1:rows/2,1) = -1;
        calc(((rows+1)/2):rows,1) =1;
        ratings(:,1) = calc(:,2);
        ratings(:,2) = calc(:,1);      
    else
        if partition == 1
            ratings = scores.one;
        elseif partition == 2
            ratings = scores.two;
        elseif partition == 3
            ratings = scores.three;
        end
    end
    
    new_participants = {};
    cnt = 1;
    for j = 1:numel(participant_order)
        participant = participant_order{j};
        score = ratings(find(ratings(:,1)==participant),2);

        if numel(score) > 0
            design(1,cnt)=score;
            new_participants{cnt} = participants{j};
            cnt = cnt + 1;
        end
    end
    
end

%% return scores
function scores = return_scores(regression_type, type_of_effect)
    if strcmp(regression_type, 'no-factor')
        dataset = [
        1,1;2,1;3,1;4,1;5,1;6,1;7,1;8,1;9,1;10,1;11,1;12,1;13,1;14,1;15,1;
        16,1;17,1;18,1;19,1;20,1;21,1;22,1;23,1;24,1;25,1;26,1;27,1;28,1;
        29,1;30,1;31,1;32,1;33,1;34,1;35,1;36,1;37,1;38,1;39,1;40,1;
        ];

        scores.one = dataset;
        scores.two = dataset;
        scores.three = dataset;
    
    elseif strcmp(regression_type, 'headache_pow')
        scores.one = [
            1,-0.2885;2,-1.8595;3,10.9627;4,6.9245;5,8.6422;6,-0.1374;7,-0.3267;
            8,5.1259;9,-0.5626;10,0.8134;11,2.8146;12,3.2084;14,11.449;16,-1.9692;
            17,14.0609;20,2.289;21,5.2967;22,0.8743;23,0.2754;24,-3.1482;26,-8.6788;
            28,-5.1555;29,2.0135;30,3.1227;31,-1.608;32,-0.3722;33,-0.9248;36,-1.3025;
            37,-8.4987;38,-3.7692;39,1.4244;
        ];
    
        scores.two = [
            1,-1.23;2,-6.5894;3,-3.7494;4,0.0093;5,1.1629;6,0.2218;7,0.0675;
            8,6.4038;9,1.3281;10,-2.5226;11,0.0225;12,0.1993;14,12.0323;16,-3.813;
            17,6.2715;20,2.082;21,5.7632;22,0.538;23,-0.3641;24,20.5365;26,-3.53;
            28,-0.4075;29,-0.9858;30,4.5612;31,-9.9556;32,-0.252;33,1.3657;36,-0.2272;
            37,8.9155;38,1.533;39,-5.7536;
        ];
    
        scores.three = [
            1,-2.5159;2,6.135;3,2.3024;4,-6.4444;5,-4.9224;6,-1.026;7,-0.2189;
            8,-2.0872;9,-0.792;10,-3.7518;11,-4.4468;12,-0.3598;14,9.0351;16,0.7133;
            17,-12.8217;20,-1.192;21,-1.7435;22,0.2708;23,-0.2506;24,5.7662;26,-0.5664;
            28,7.0665;29,0.5425;30,2.6584;31,0.8058;32,3.4491;33,-6.0637;36,-3.9158;37,-7.1591;
            38,7.4466;39,-8.1057;
        ];
  
        return;
        
        min_n = min([scores.one, scores.two, scores.three]);
        scores.one(:,2) = scores.one(:,2) - min_n(2);
        scores.two(:,2) = scores.two(:,2) - min_n(2);
        scores.three(:,2) = scores.three(:,2) - min_n(2);
    
        if strcmp(type_of_effect, 'habituation')
            scores.one(:,2) = scores.one(:,2) * 2.72;
            scores.two(:,2) = scores.two(:,2) * 1.65;
            scores.three(:,2) = scores.three(:,2) * 1.00;
        elseif strcmp(type_of_effect, 'sensitization')
            scores.one(:,2) = scores.one(:,2) * 1.00;
            scores.two(:,2) = scores.two(:,2) * 1.65;
            scores.three(:,2) = scores.three(:,2) * 2.72;
        else
            error('Type of experiment not properly specified');
        end
        
        [n_participants, ~] = size(scores.one);
        
        for k=1:n_participants
            p1 = scores.one(k,1);
            p2 = scores.two(k,1);
            p3 = scores.three(k,1);
            
            if p1 == p2 && p2 == p3                
                if strcmp(type_of_effect, 'habituation')
                    to_remove = scores.three(k,2);
                elseif strcmp(type_of_effect, 'sensitization')
                    to_remove = scores.one(k,2);
                end
                
                scores.one(k,2) = scores.one(k,2) - to_remove;
                scores.two(k,2) = scores.two(k,2) - to_remove;
                scores.three(k,2) = scores.three(k,2) - to_remove;
            else
                error('Participants do not align...')
            end
        end    
        
    elseif strcmp(regression_type, 'headache')
        dataset = [
        1,-0.2574;2,-0.0417;3,-0.6726;4,0.4236;5,1.781;6,-1.0608;7,-0.7657;
        8,0.1279;9,-0.6553;10,-0.2896;11,-0.5122;12,2.1424;13,-0.1803;
        14,1.4491;16,0.1157;17,-0.1649;20,-0.4721;21,1.0486;22,-0.554;23,-0.8912;
        24,-0.4481;25,-0.7581;26,-1.2784;28,0.2989;29,0.0439;30,-0.4732;31,-0.7701;
        32,-0.7037;33,-0.819;34,-0.7987;37,1.1507;38,-0.2806;39,0.8546;40,-0.3823;   
        ];
    
        scores.one = dataset;
        scores.two = dataset;
        scores.three = dataset;
    
        min_n = min(scores.one);
        scores.one(:,2) = scores.one(:,2) - min_n(2);
        scores.two(:,2) = scores.two(:,2) - min_n(2);
        scores.three(:,2) = scores.three(:,2) - min_n(2);
        
        if strcmp(type_of_effect, 'habituation')
            scores.one(:,2) = scores.one(:,2) * 2.72;
            scores.two(:,2) = scores.two(:,2) * 1.65;
            scores.three(:,2) = scores.three(:,2) * 1.00;
        elseif strcmp(type_of_effect, 'sensitization')
            scores.one(:,2) = scores.one(:,2) * 1.00;
            scores.two(:,2) = scores.two(:,2) * 1.65;
            scores.three(:,2) = scores.three(:,2) * 2.72;
        else
            error('Type of experiment not properly specified');
        end
        
        [n_participants, ~] = size(dataset);
        
        for k=1:n_participants
            p1 = scores.one(k,1);
            p2 = scores.two(k,1);
            p3 = scores.three(k,1);
            
            if p1 == p2 && p2 == p3                
                if strcmp(type_of_effect, 'habituation')
                    to_remove = scores.three(k,2);
                elseif strcmp(type_of_effect, 'sensitization')
                    to_remove = scores.one(k,2);
                end
                
                scores.one(k,2) = scores.one(k,2) - to_remove;
                scores.two(k,2) = scores.two(k,2) - to_remove;
                scores.three(k,2) = scores.three(k,2) - to_remove;
            else
                error('Participants do not align...')
            end
        end     
                   
    elseif strcmp(regression_type, 'visual_stress')
        dataset = [
        1,0.3227;2,-0.1086;3,-0.5102;4,1.1336;5,-0.6395;6,-1.2147;
        7,-0.3301;8,0.7524;9,-0.3903;10,-0.7221;11,-0.769;12,-1.063;
        13,-0.8985;14,-1.4672;16,-1.1987;17,0.1542;20,0.5867;21,1.0008;
        22,-0.1169;23,1.7209;24,0.1411;25,0.6221;26,-0.7483;28,0.6739;
        29,-0.0237;30,0.0364;31,0.6996;32,-0.2998;33,-0.65;34,0.0262;
        37,0.7986;38,-0.5883;39,2.3332;40,2.2667;   
        ];

        scores.one = dataset;
        scores.two = dataset;
        scores.three = dataset;
    
        min_n = min(scores.one);
        scores.one(:,2) = scores.one(:,2) - min_n(2);
        scores.two(:,2) = scores.two(:,2) - min_n(2);
        scores.three(:,2) = scores.three(:,2) - min_n(2);
    
        if strcmp(type_of_effect, 'habituation')
            scores.one(:,2) = scores.one(:,2) * 2.72;
            scores.two(:,2) = scores.two(:,2) * 1.65;
            scores.three(:,2) = scores.three(:,2) * 1.00;
        elseif strcmp(type_of_effect, 'sensitization')
            scores.one(:,2) = scores.one(:,2) * 1.00;
            scores.two(:,2) = scores.two(:,2) * 1.65;
            scores.three(:,2) = scores.three(:,2) * 2.72;
        else
            error('Type of experiment not properly specified');
        end
        
        [n_participants, ~] = size(dataset);
        
        for k=1:n_participants
            p1 = scores.one(k,1);
            p2 = scores.two(k,1);
            p3 = scores.three(k,1);
            
            if p1 == p2 && p2 == p3
                if strcmp(type_of_effect, 'habituation')
                    to_remove = scores.three(k,2);
                elseif strcmp(type_of_effect, 'sensitization')
                    to_remove = scores.one(k,2);
                end
                
                scores.one(k,2) = scores.one(k,2) - to_remove;
                scores.two(k,2) = scores.two(k,2) - to_remove;
                scores.three(k,2) = scores.three(k,2) - to_remove;
            else
                error('Participants do not align...')
            end
        end

    elseif strcmp(regression_type, 'discomfort')
        dataset = [
        1,-0.2427;2,0.4398;3,-0.5221;4,1.8399;5,-0.6095;6,0.8092;7,-0.6979;
        8,0.9717;9,-0.8232;10,-0.9152;11,0.3501;12,-0.8418;13,-0.7414;
        14,1.5086;16,1.0678;17,1.5466;20,0.1606;21,0.1343;22,0.6145;23,-1.3703;
        24,2.2964;25,-0.7656;26,-0.5905;28,-0.8957;29,0.3773;30,-0.6245;31,2.1948;
        32,-1.5111;33,1.1882;34,-0.7889;37,-0.5762;38,-1.0582;39,-0.7461;40,0.5811;
        ];
    
        scores.one = dataset;
        scores.two = dataset;
        scores.three = dataset;

        min_n = min(scores.one);
        scores.one(:,2) = scores.one(:,2) - min_n(2);
        scores.two(:,2) = scores.two(:,2) - min_n(2);
        scores.three(:,2) = scores.three(:,2) - min_n(2);
    
    
        if strcmp(type_of_effect, 'habituation')
            scores.one(:,2) = scores.one(:,2) * 2.72;
            scores.two(:,2) = scores.two(:,2) * 1.65;
            scores.three(:,2) = scores.three(:,2) * 1.00;
        elseif strcmp(type_of_effect, 'sensitization')
            scores.one(:,2) = scores.one(:,2) * 1.00;
            scores.two(:,2) = scores.two(:,2) * 1.65;
            scores.three(:,2) = scores.three(:,2) * 2.72;
        else
            error('Type of experiment not properly specified');
        end
        
        [n_participants, ~] = size(dataset);
        
        for k=1:n_participants
            p1 = scores.one(k,1);
            p2 = scores.two(k,1);
            p3 = scores.three(k,1);
            
            if p1 == p2 && p2 == p3
                if strcmp(type_of_effect, 'habituation')
                    to_remove = scores.three(k,2);
                elseif strcmp(type_of_effect, 'sensitization')
                    to_remove = scores.one(k,2);
                end
                
                scores.one(k,2) = scores.one(k,2) - to_remove;
                scores.two(k,2) = scores.two(k,2) - to_remove;
                scores.three(k,2) = scores.three(k,2) - to_remove;
            else
                error('Participants do not align...')
            end
        end

    elseif contains(regression_type, 'headache-scores-entire')
        scores.one = [
        1,-0.227;2,-0.052;3,-0.721;4,0.531;5,1.72;6,-1.176;7,-0.797;
        8,0.199;9,-0.692;10,-0.358;11,-0.585;12,2.041;13,-0.266;14,1.31;
        16,0.002;17,-0.15;20,-0.416;21,1.144;22,-0.565;23,-0.728;24,-0.435;
        25,-0.699;26,-1.35;28,0.363;29,0.042;30,-0.47;31,-0.704;32,-0.732;
        33,-0.881;34,-0.796;37,1.227;38,-0.337;39,1.077;40,-0.167
        ];
    
    elseif contains(regression_type, 'discomfort-scores-entire')
        scores.one = [
        1,-0.216;2,0.439;3,-0.495;4,1.848;5,-0.746;6,0.84;7,-0.659;8,0.987;
        9,-0.793;10,-0.919;11,0.359;12,-1.015;13,-0.758;14,1.368;16,1.022;
        17,1.562;20,0.21;21,0.098;22,0.647;23,-1.258;26,-0.531;28,-0.894;
        29,0.374;30,-0.593;31,2.267;32,-1.475;33,1.221;34,-0.736;38,-1.059;
        39,-0.728;  
        ];
       
    elseif contains(regression_type, 'vs-scores-entire')
        scores.one = [
        1,0.323;2,-0.109;3,-0.51;4,1.134;5,-0.639;6,-1.215;7,-0.33;8,0.752;
        9,-0.39;10,-0.722;11,-0.769;12,-1.063;13,-0.899;14,-1.467;16,-1.199;
        17,0.154;20,0.587;21,1.001;22,-0.117;23,1.721;26,-0.748;28,0.674;
        29,-0.024;30,0.036;31,0.7;32,-0.3;33,-0.65;34,0.026;38,-0.588;39,2.333;
        ];
    end
end
%% load post-processed fildtrip data
function [labels, ft_regression_data, participant_order] = ...
    load_data_for_plotting(main_path, n_participants, filename, partition, peak_electrode, analysis_window)


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
            
            e_idx = find(contains(data.label,peak_electrode));
            
            if partition.is_partition
               if partition.partition_number == 1
                   med = data.p1_med(e_idx,:,:);
                   thick = data.p1_thick(e_idx,:,:);
                   thin = data.p1_thin(e_idx,:,:);
               elseif partition.partition_number == 2
                   med = data.p2_med(e_idx,:,:);
                   thick = data.p2_thick(e_idx,:,:);
                   thin = data.p2_thin(e_idx,:,:);
               elseif partition.partition_number == 3
                   med = data.p3_med(e_idx,:,:);
                   thick = data.p3_thick(e_idx,:,:);
                   thin = data.p3_thin(e_idx,:,:);
               end
            elseif ~partition.is_partition
                med = data.med(e_idx,:,:);
                thin = data.thin(e_idx,:,:);
                thick = data.thick(e_idx,:,:);
            end
            
            ft.med_t = med;
            ft.thin_t = thin;
            ft.thick_t = thick;
            ft.time = data.time{1};
            
            ft.med = mean(med,3);
            ft.thick = mean(thick, 3);
            ft.thin = mean(thin, 3);
            ft.avg = ft.med - (ft.thick + ft.thin)/2;
            
            ft = cut_data_using_analysis_window(ft, analysis_window);
            
            
            ft_regression_data{idx_used_for_saving_data} = ft;
            participant_order{idx_used_for_saving_data} = i;
            idx_used_for_saving_data = idx_used_for_saving_data + 1;
            
            labels = data.label;
        end
    end
end

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

%% generate erp plots
function generate_plots(master_dir, main_path, experiment_type, start_peak, ...
    end_peak, peak_electrode, peak_effect, t_value, df, regression_type, ...
    pvalue, cluster_size, save_dir, effect_type, weight_erps, weighting_factor);

    
    plotting_window = [2800, 4000];
    rmpath W:\PhD\MatlabPlugins\spm8;
    addpath W:\PhD\MatlabPlugins\spm12;
    cd(master_dir);

    %% Are we looking at onsets 2-8 or partitions
    % set up the experiment as needed
    if strcmp(experiment_type, 'onsets-2-8-explicit')
        n_participants = 40;

        partition.is_partition = 0;
        partition.partition_number = 0;


        data_file = 'time_domain_mean_intercept_onsets_2_3_4_5_6_7_8_grand-average_baseline_2.mat';
       [data, ~] = load_postprocessed_data(main_path, n_participants, ...
            data_file, partition);
        e_idx = find(contains(data{1}.label,peak_electrode));
        ci = bootstrap_erps(data, e_idx);
        
    elseif strcmp(experiment_type, 'partitions (no factor)') || strcmp(experiment_type, 'partitions-1') ||  strcmp(experiment_type, 'partitions-2-8')
        n_participants = 39;

        partition1.is_partition = 1; 
        partition1.partition_number = 1;
        partition2.is_partition = 1; 
        partition2.partition_number = 2;
        partition3.is_partition = 1; 
        partition3.partition_number = 3;

        data_file = 'partitions_partitioned_onsets_2_3_4_5_6_7_8_grand-average.mat';
        [data1, participant_order_1] = load_postprocessed_data(main_path, n_participants, ...
            data_file, partition1);
        e_idx = find(contains(data1{1}.label,peak_electrode));
        [data2, participant_order_2] = load_postprocessed_data(main_path, n_participants, ...
            data_file, partition2);
        [data3, participant_order_3] = load_postprocessed_data(main_path, n_participants, ...
            data_file, partition3);
        
        data = [data1,data2,data3];
        
        type_of_effect = 'habituation';
        [data1_h, data1_l] = get_partitions_medium_split(data1, participant_order_1,...
            regression_type, 1, type_of_effect, weight_erps, weighting_factor);
        ci1_h = bootstrap_erps(data1_h, e_idx);
        ci1_l = bootstrap_erps(data1_l, e_idx);
        [data2_h, data2_l] = get_partitions_medium_split(data2, participant_order_2,...
            regression_type, 2, type_of_effect, weight_erps, weighting_factor);
        ci2_h = bootstrap_erps(data2_h, e_idx);
        ci2_l = bootstrap_erps(data2_l, e_idx);
        [data3_h, data3_l] = get_partitions_medium_split(data3, participant_order_3,...
            regression_type, 3, type_of_effect, weight_erps, weighting_factor);
        ci3_h = bootstrap_erps(data3_h, e_idx);
        ci3_l = bootstrap_erps(data3_l, e_idx);

        
    elseif strcmp(experiment_type, 'erps-23-45-67') || strcmp(experiment_type, 'erps-23-45-67-no-factor') 
        type_of_effect = 'sensitization';
        data_file23 = 'mean_intercept_onsets_2_3_grand-average.mat';
        data_file45 = 'mean_intercept_onsets_4_5_grand-average.mat';
        data_file67 = 'mean_intercept_onsets_6_7_grand-average.mat';
            
        n_participants = 40;
        partition.is_partition = 0;
        partition.partition_number = 0;

        [data1, participant_order_1] = load_postprocessed_data(main_path, n_participants, ...
            data_file23, partition);
        e_idx = find(contains(data1{1}.label,peak_electrode));
        [data2, participant_order_2] = load_postprocessed_data(main_path, n_participants, ...
            data_file45, partition);
        [data3, participant_order_3] = load_postprocessed_data(main_path, n_participants, ...
            data_file67, partition);

        
        [data1_h, data1_l] = get_partitions_medium_split(data1, participant_order_1,...
            regression_type, 1, type_of_effect, weight_erps, weighting_factor);
        ci1_h = bootstrap_erps(data1_h,e_idx);
        ci1_l = bootstrap_erps(data1_l,e_idx);
        [data2_h, data2_l] = get_partitions_medium_split(data2, participant_order_2,...
            regression_type, 1, type_of_effect, weight_erps, weighting_factor);
        ci2_h = bootstrap_erps(data2_h,e_idx);
        ci2_l = bootstrap_erps(data2_l,e_idx);
        [data3_h, data3_l] = get_partitions_medium_split(data3, participant_order_3,...
            regression_type, 1,  type_of_effect, weight_erps, weighting_factor);
        ci3_h = bootstrap_erps(data3_h,e_idx);
        ci3_l = bootstrap_erps(data3_l,e_idx);
        
        data = [data1, data2, data3];
    end
    %% generate_supplementary information and indices used to plot
    if strcmp(experiment_type, 'partitions-2-8')
        experiment_name = 'Partitions (2:8) ';
    elseif strcmp(experiment_type, 'erps-23-45-67')
        experiment_name = 'Onsets (2,3; 4,5; 6,7)';
    elseif strcmp(experiment_type, 'erps-23-45-67-no-factor')
        experiment_name = 'ERPs 2,3; 4,5; 6,7 (No Factor)';
    elseif strcmp(experiment_type, 'onsets-2-8-explicit')
        experiment_name = 'Onsets 2:8 Mean Intercept';
    else
       experiment_name = experiment_type;
    end
    
    
    regression_type = regexprep(regression_type,'(\<[a-z])','${upper($1)}');
    effect_type = strcat(regexprep(effect_type,'(\<[a-z])','${upper($1)}'), ' Tail');
    mtitle = strcat(effect_type, {' '}, experiment_name, {' ' }, regression_type, {' '}, peak_electrode);
    mtitle = mtitle{1};
    
    start_peak = start_peak*1000;
    end_peak = end_peak*1000;
    cohens_d = round((2*t_value)/sqrt(df),2);
    effect_size = round(sqrt((t_value*t_value)/((t_value*t_value)+df)),2);
    
    time = data{1}.time * 1000;
    peak_effect = peak_effect*1000;
    t_value = round(t_value, 2);
    cluster_size = round(cluster_size, 0);
    
    woody_filt = 0;
    if woody_filt == 1
        t = time(1:257);
        t = tiledlayout(1,2, 'TileSpacing','Compact');
        nexttile
        p1_pgi_h = woody_filtering(data1_h, e_idx);
        p2_pgi_h = woody_filtering(data2_h, e_idx);
        p3_pgi_h = woody_filtering(data3_h, e_idx);
        
        plot(p1_pgi_h, 'r');
        hold;
        plot(p2_pgi_h, 'g');
        plot(p3_pgi_h, 'b');
        legend({'P1', 'P2', 'P3'},'Location','northeast')
        nexttile;
        
        p1_pgi_l = woody_filtering(data1_l, e_idx);
        p2_pgi_l = woody_filtering(data2_l, e_idx);
        p3_pgi_l = woody_filtering(data3_l, e_idx);
        
        plot(p1_pgi_l, 'r');
        hold;
        plot(p2_pgi_l, 'g');
        plot(p3_pgi_l, 'b');
        legend({'P1', 'P2', 'P3'},'Location','northeast')
        nexttile;
    end
    
    
    t = tiledlayout(5,2, 'TileSpacing','Compact');

    if contains(experiment_type, 'onsets-2-8')
       t = tiledlayout(2,1, 'TileSpacing','Compact');
       time = data{1}.time * 1000;
       nexttile
       hold on;
       plot(time, ci.dist_pgi_avg, 'color', 'r', 'LineWidth', 1.35,'DisplayName','PGI')
       plot(time, ci.dist_pgi_high, 'LineWidth', 0.01, 'color', 'r','DisplayName','');
       plot(time, ci.dist_pgi_low, 'LineWidth', 0.01, 'color', 'r','DisplayName','');
       x2 = [time, fliplr(time)];
       inBetween = [ci.dist_pgi_high, fliplr(ci.dist_pgi_low)];
       h = fill(x2, inBetween, 'r');
       set(h,'facealpha',.13)
       xlim(plotting_window);
       ylim([-5, 8])
       grid on
       xline(start_peak, '-');
       xline(end_peak, '-');
       xline(peak_effect, '--r');
       legend({'PGI'},'Location','northeast')

        nexttile
        hold on;
        
        plot(NaN(1), 'g');
        plot(NaN(1), 'b');
        plot(NaN(1), 'r');
        legend({'Thin', 'Medium', 'Thick'},'Location','northeast')
        
        plot(time, ci.dist_thin_avg, 'color', 'g', 'LineWidth', 1.35, 'HandleVisibility','off')
        plot(time, ci.dist_thin_high, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
        plot(time, ci.dist_thin_low, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
        x2 = [time, fliplr(time)];
        inBetween = [ci.dist_thin_high, fliplr(ci.dist_thin_low)];
        h = fill(x2, inBetween, 'g', 'HandleVisibility','off');
        set(h,'facealpha',.13)
       
        plot(time, ci.dist_med_avg, 'color', 'b','LineWidth', 1.35, 'HandleVisibility','off');
        plot(time, ci.dist_med_high, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
        plot(time, ci.dist_med_low, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
        x2 = [time, fliplr(time)];
        inBetween = [ci.dist_med_high, fliplr(ci.dist_med_low)];
        h = fill(x2, inBetween, 'b', 'HandleVisibility','off');
        set(h,'facealpha',.13)
       

        plot(time, ci.dist_thick_avg, 'color', 'r','LineWidth', 1.35, 'HandleVisibility','off');
        plot(time, ci.dist_thick_high, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
        plot(time, ci.dist_thick_low, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
        x2 = [time, fliplr(time)];
        inBetween = [ci.dist_thick_high, fliplr(ci.dist_thick_low)];
        h = fill(x2, inBetween, 'r', 'HandleVisibility','off');
        set(h,'facealpha',.13)
       
       xline(start_peak, '-', 'HandleVisibility','off');
       xline(end_peak, '-', 'HandleVisibility','off');
       xline(peak_effect, '--r', 'HandleVisibility','off');
       xlim(plotting_window);
       ylim([-6, 8])
       grid on
     

       hold off;
    elseif strcmp(experiment_type, 'partitions-2-8') || strcmp(experiment_type, 'erps-23-45-67')
       time = data{1}.time * 1000;

       % PGI HIGH
       nexttile
       hold on;
       plot(NaN(1), 'r');
       plot(NaN(1), 'g');
       plot(NaN(1), 'b');
       if contains(experiment_type, 'partitions-2-8')
        legend({'P1-PGI', 'P2-PGI', 'P3-PGI'},'Location','northeast')
       elseif strcmp(experiment_type, 'erps-23-45-67')
           legend({'Onsets 2:3', 'Onsets 4:5', 'Onsets 6:7'},'Location','northeast')
       end
       
       plot(time, ci1_h.dist_pgi_avg, 'color', 'r', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci1_h.dist_pgi_high, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       plot(time, ci1_h.dist_pgi_low, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci1_h.dist_pgi_high, fliplr(ci1_h.dist_pgi_low)];
       h = fill(x2, inBetween, 'r', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci2_h.dist_pgi_avg, 'color', 'g', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci2_h.dist_pgi_high, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       plot(time, ci2_h.dist_pgi_low, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci2_h.dist_pgi_high, fliplr(ci2_h.dist_pgi_low)];
       h = fill(x2, inBetween, 'g', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci3_h.dist_pgi_avg, 'color', 'b', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci3_h.dist_pgi_high, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       plot(time, ci3_h.dist_pgi_low, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci3_h.dist_pgi_high, fliplr(ci3_h.dist_pgi_low)];
       h = fill(x2, inBetween, 'b', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       xlim(plotting_window);
       title('High Group: Partitions: PGI');
       ylim([-6, 12])
       grid on;
       hold off;
       
       xline(start_peak, '-','HandleVisibility','off');
       xline(end_peak, '-','HandleVisibility','off');
       xline(peak_effect, '--r','HandleVisibility','off');
      
       % PGI LOW
       nexttile
       hold on;
       plot(NaN(1), 'r');
       plot(NaN(1), 'g');
       plot(NaN(1), 'b');
       if contains(experiment_type, 'partitions-2-8')
        legend({'P1-PGI', 'P2-PGI', 'P3-PGI'},'Location','northeast')
       elseif strcmp(experiment_type, 'erps-23-45-67')
           legend({'Onsets 2:3', 'Onsets 4:5', 'Onsets 6:7'},'Location','northeast')
       end
       
       plot(time, ci1_l.dist_pgi_avg, 'color', 'r', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci1_l.dist_pgi_high, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       plot(time, ci1_l.dist_pgi_low, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci1_l.dist_pgi_high, fliplr(ci1_l.dist_pgi_low)];
       h = fill(x2, inBetween, 'r', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci2_l.dist_pgi_avg, 'color', 'g', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci2_l.dist_pgi_high, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       plot(time, ci2_l.dist_pgi_low, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci2_l.dist_pgi_high, fliplr(ci2_l.dist_pgi_low)];
       h = fill(x2, inBetween, 'g', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci3_l.dist_pgi_avg, 'color', 'b', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci3_l.dist_pgi_high, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       plot(time, ci3_l.dist_pgi_low, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci3_l.dist_pgi_high, fliplr(ci3_l.dist_pgi_low)];
       h = fill(x2, inBetween, 'b', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       xlim(plotting_window);
       title('Low Group: Partitions: PGI');
       ylim([-6, 12])
       grid on;
       hold off;
       
       xline(start_peak, '-','HandleVisibility','off');
       xline(end_peak, '-','HandleVisibility','off');
       xline(peak_effect, '--r','HandleVisibility','off');    
       
       % MED HIGH
       nexttile
       hold on;
       plot(NaN(1), 'r');
       plot(NaN(1), 'g');
       plot(NaN(1), 'b');
       if contains(experiment_type, 'partitions-2-8')
            legend({'Med-P1', 'Med-P2', 'Med-P3'},'Location','northeast')
       elseif strcmp(experiment_type, 'erps-23-45-67')
            legend({'Med-2:3', 'Med-4:5', 'Med-6:7'},'Location','northeast')
       end
       
       
       plot(time, ci1_h.dist_med_avg, 'color', 'r', 'LineWidth', 1.35,'HandleVisibility','off');
       %plot(time, ci1_h.dist_med_high, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       %plot(time, ci1_h.dist_med_low, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       %x2 = [time, fliplr(time)];
       %inBetween = [ci1_h.dist_med_high, fliplr(ci1_h.dist_med_low)];
       %h = fill(x2, inBetween, 'r', 'HandleVisibility','off');
       %set(h,'facealpha',.13)
       
       plot(time, ci2_h.dist_med_avg, 'color', 'g', 'LineWidth', 1.35,'HandleVisibility','off');
       %plot(time, ci2_h.dist_med_high, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       %plot(time, ci2_h.dist_med_low, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       %x2 = [time, fliplr(time)];
       %inBetween = [ci2_h.dist_med_high, fliplr(ci2_h.dist_med_low)];
       %h = fill(x2, inBetween, 'g', 'HandleVisibility','off');
       %set(h,'facealpha',.13)
       
       plot(time, ci3_h.dist_med_avg, 'color', 'b', 'LineWidth', 1.35,'HandleVisibility','off');
       %plot(time, ci3_h.dist_med_high, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       %plot(time, ci3_h.dist_med_low, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       %x2 = [time, fliplr(time)];
       %inBetween = [ci3_h.dist_med_high, fliplr(ci3_h.dist_med_low)];
       %h = fill(x2, inBetween, 'b', 'HandleVisibility','off');
       %set(h,'facealpha',.13)
       
       xlim(plotting_window);
       title('High Group: Medium Through the Partitions');
       ylim([-6, 12])
       grid on;
       hold off;
       
       xline(start_peak, '-','HandleVisibility','off');
       xline(end_peak, '-','HandleVisibility','off');
       xline(peak_effect, '--r','HandleVisibility','off');
       
       % MED LOW
       nexttile
       hold on;
       plot(NaN(1), 'r');
       plot(NaN(1), 'g');
       plot(NaN(1), 'b');
       if contains(experiment_type, 'partitions-2-8')
            legend({'Med-P1', 'Med-P2', 'Med-P3'},'Location','northeast')
       elseif strcmp(experiment_type, 'erps-23-45-67')
            legend({'Med-2:3', 'Med-4:5', 'Med-6:7'},'Location','northeast')
       end
       
       plot(time, ci1_l.dist_med_avg, 'color', 'r', 'LineWidth', 1.35,'HandleVisibility','off');
       %plot(time, ci1_l.dist_med_high, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       %plot(time, ci1_l.dist_med_low, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       %x2 = [time, fliplr(time)];
       %inBetween = [ci1_l.dist_med_high, fliplr(ci1_l.dist_med_low)];
       %h = fill(x2, inBetween, 'r', 'HandleVisibility','off');
       %set(h,'facealpha',.13)
       
       plot(time, ci2_l.dist_med_avg, 'color', 'g', 'LineWidth', 1.35,'HandleVisibility','off');
       %plot(time, ci2_l.dist_med_high, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       %plot(time, ci2_l.dist_med_low, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       %x2 = [time, fliplr(time)];
       %inBetween = [ci2_l.dist_med_high, fliplr(ci2_l.dist_med_low)];
       %h = fill(x2, inBetween, 'g', 'HandleVisibility','off');
       %set(h,'facealpha',.13)
       
       plot(time, ci3_l.dist_med_avg, 'color', 'b', 'LineWidth', 1.35,'HandleVisibility','off');
       %plot(time, ci3_l.dist_med_high, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       %plot(time, ci3_l.dist_med_low, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       %x2 = [time, fliplr(time)];
       %inBetween = [ci3_l.dist_med_high, fliplr(ci3_l.dist_med_low)];
       %h = fill(x2, inBetween, 'b', 'HandleVisibility','off');
       %set(h,'facealpha',.13)
       
       xlim(plotting_window);
       title('Low Group: Medium Through the Partitions');
       ylim([-6, 12])
       grid on;
       hold off;
       
       xline(start_peak, '-','HandleVisibility','off');
       xline(end_peak, '-','HandleVisibility','off');
       xline(peak_effect, '--r','HandleVisibility','off');
       
       % P1 HIGH
       nexttile
       hold on;
       plot(NaN(1), 'r');
       plot(NaN(1), 'g');
       plot(NaN(1), 'b');
       legend({'Thin', 'Medium', 'Thick'},'Location','northeast')
       
       plot(time, ci1_h.dist_thin_avg, 'color', 'r', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci1_h.dist_thin_high, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       plot(time, ci1_h.dist_thin_low, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci1_h.dist_thin_high, fliplr(ci1_h.dist_thin_low)];
       h = fill(x2, inBetween, 'r', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci1_h.dist_med_avg, 'color', 'g', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci1_h.dist_med_high, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       plot(time, ci1_h.dist_med_low, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci1_h.dist_med_high, fliplr(ci1_h.dist_med_low)];
       h = fill(x2, inBetween, 'g', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci1_h.dist_thick_avg, 'color', 'b', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci1_h.dist_thick_high, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       plot(time, ci1_h.dist_thick_low, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci1_h.dist_thick_high, fliplr(ci1_h.dist_thick_low)];
       h = fill(x2, inBetween, 'b', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       xlim(plotting_window);
       title('High Group');
       ylim([-6, 12])
       grid on;
       hold off;
       
       xline(start_peak, '-','HandleVisibility','off');
       xline(end_peak, '-','HandleVisibility','off');
       xline(peak_effect, '--r','HandleVisibility','off');
       
       % P1 LOW
       nexttile
       hold on;
       plot(NaN(1), 'r');
       plot(NaN(1), 'g');
       plot(NaN(1), 'b');
       legend({'Thin', 'Medium', 'Thick'},'Location','northeast')
       
       plot(time, ci1_l.dist_thin_avg, 'color', 'r', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci1_l.dist_thin_high, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       plot(time, ci1_l.dist_thin_low, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci1_l.dist_thin_high, fliplr(ci1_l.dist_thin_low)];
       h = fill(x2, inBetween, 'r', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci1_l.dist_med_avg, 'color', 'g', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci1_l.dist_med_high, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       plot(time, ci1_l.dist_med_low, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci1_l.dist_med_high, fliplr(ci1_l.dist_med_low)];
       h = fill(x2, inBetween, 'g', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci1_l.dist_thick_avg, 'color', 'b', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci1_l.dist_thick_high, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       plot(time, ci1_l.dist_thick_low, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci1_l.dist_thick_high, fliplr(ci1_l.dist_thick_low)];
       h = fill(x2, inBetween, 'b', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       xlim(plotting_window);
       title('Low Group');
       ylim([-6, 12])
       grid on;
       hold off;
       
       xline(start_peak, '-','HandleVisibility','off');
       xline(end_peak, '-','HandleVisibility','off');
       xline(peak_effect, '--r','HandleVisibility','off');
       
       % P2 HIGH
       nexttile
       hold on;
       plot(NaN(1), 'r');
       plot(NaN(1), 'g');
       plot(NaN(1), 'b');
       legend({'Thin', 'Medium', 'Thick'},'Location','northeast')
       
       plot(time, ci2_h.dist_thin_avg, 'color', 'r', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci2_h.dist_thin_high, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       plot(time, ci2_h.dist_thin_low, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci2_h.dist_thin_high, fliplr(ci2_h.dist_thin_low)];
       h = fill(x2, inBetween, 'r', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci2_h.dist_med_avg, 'color', 'g', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci2_h.dist_med_high, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       plot(time, ci2_h.dist_med_low, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci2_h.dist_med_high, fliplr(ci2_h.dist_med_low)];
       h = fill(x2, inBetween, 'g', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci2_h.dist_thick_avg, 'color', 'b', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci2_h.dist_thick_high, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       plot(time, ci2_h.dist_thick_low, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci2_h.dist_thick_high, fliplr(ci2_h.dist_thick_low)];
       h = fill(x2, inBetween, 'b', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       xlim(plotting_window);
       title('High Group');
       ylim([-6, 12])
       grid on;
       hold off;
       
       xline(start_peak, '-','HandleVisibility','off');
       xline(end_peak, '-','HandleVisibility','off');
       xline(peak_effect, '--r','HandleVisibility','off');
       
       % P2 LOW
       nexttile
       hold on;
       plot(NaN(1), 'r');
       plot(NaN(1), 'g');
       plot(NaN(1), 'b');
       legend({'Thin', 'Medium', 'Thick'},'Location','northeast')
       
       plot(time, ci2_l.dist_thin_avg, 'color', 'r', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci2_l.dist_thin_high, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       plot(time, ci2_l.dist_thin_low, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci2_l.dist_thin_high, fliplr(ci2_l.dist_thin_low)];
       h = fill(x2, inBetween, 'r', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci2_l.dist_med_avg, 'color', 'g', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci2_l.dist_med_high, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       plot(time, ci2_l.dist_med_low, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci2_l.dist_med_high, fliplr(ci2_l.dist_med_low)];
       h = fill(x2, inBetween, 'g', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci2_l.dist_thick_avg, 'color', 'b', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci2_l.dist_thick_high, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       plot(time, ci2_l.dist_thick_low, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci2_l.dist_thick_high, fliplr(ci2_l.dist_thick_low)];
       h = fill(x2, inBetween, 'b', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       xlim(plotting_window);
       title('Low Group');
       ylim([-6, 12])
       grid on;
       hold off;
       
       xline(start_peak, '-','HandleVisibility','off');
       xline(end_peak, '-','HandleVisibility','off');
       xline(peak_effect, '--r','HandleVisibility','off');
       
       % P3 HIGH 
       nexttile
       hold on;
       plot(NaN(1), 'r');
       plot(NaN(1), 'g');
       plot(NaN(1), 'b');
       legend({'Thin', 'Medium', 'Thick'},'Location','northeast')
       
       plot(time, ci3_h.dist_thin_avg, 'color', 'r', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci3_h.dist_thin_high, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       plot(time, ci3_h.dist_thin_low, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci3_h.dist_thin_high, fliplr(ci3_h.dist_thin_low)];
       h = fill(x2, inBetween, 'r', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci3_h.dist_med_avg, 'color', 'g', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci3_h.dist_med_high, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       plot(time, ci3_h.dist_med_low, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci3_h.dist_med_high, fliplr(ci3_h.dist_med_low)];
       h = fill(x2, inBetween, 'g', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci3_h.dist_thick_avg, 'color', 'b', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci3_h.dist_thick_high, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       plot(time, ci3_h.dist_thick_low, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci3_h.dist_thick_high, fliplr(ci3_h.dist_thick_low)];
       h = fill(x2, inBetween, 'b', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       xlim(plotting_window);
       title('High Group');
       ylim([-6, 12])
       grid on;
       hold off;
       
       xline(start_peak, '-','HandleVisibility','off');
       xline(end_peak, '-','HandleVisibility','off');
       xline(peak_effect, '--r','HandleVisibility','off');
       
       % P3 LOW
       % P2 LOW
       nexttile
       hold on;
       plot(NaN(1), 'r');
       plot(NaN(1), 'g');
       plot(NaN(1), 'b');
       legend({'Thin', 'Medium', 'Thick'},'Location','northeast')
       
       plot(time, ci3_l.dist_thin_avg, 'color', 'r', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci3_l.dist_thin_high, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       plot(time, ci3_l.dist_thin_low, 'LineWidth', 0.01, 'color', 'r','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci3_l.dist_thin_high, fliplr(ci3_l.dist_thin_low)];
       h = fill(x2, inBetween, 'r', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci3_l.dist_med_avg, 'color', 'g', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci3_l.dist_med_high, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       plot(time, ci3_l.dist_med_low, 'LineWidth', 0.01, 'color', 'g','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci3_l.dist_med_high, fliplr(ci3_l.dist_med_low)];
       h = fill(x2, inBetween, 'g', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       plot(time, ci3_l.dist_thick_avg, 'color', 'b', 'LineWidth', 1.35,'HandleVisibility','off');
       plot(time, ci3_l.dist_thick_high, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       plot(time, ci3_l.dist_thick_low, 'LineWidth', 0.01, 'color', 'b','HandleVisibility','off');
       x2 = [time, fliplr(time)];
       inBetween = [ci3_l.dist_thick_high, fliplr(ci3_l.dist_thick_low)];
       h = fill(x2, inBetween, 'b', 'HandleVisibility','off');
       set(h,'facealpha',.13)
       
       xlim(plotting_window);
       title('Low Group');
       ylim([-6, 12])
       grid on;
       hold off;
       
       xline(start_peak, '-','HandleVisibility','off');
       xline(end_peak, '-','HandleVisibility','off');
       xline(peak_effect, '--r','HandleVisibility','off');
       
    end
    
    title(t, mtitle, 'FontSize', 14);
    effect_peak = strcat('Effect Size Correlation: ', {' '}, num2str(effect_size), ...
        ' Cohens D: ', {' '}, num2str(cohens_d), ' T: ', {' '}, num2str(t_value));
    effect_peak = effect_peak{1};
    cluster_level = strcat('P-value:', {' '}, num2str(pvalue), ' Cluster Size', ...
        {' '}, num2str(cluster_size));
    cluster_level = cluster_level{1};
    subtitle(t, {effect_peak, cluster_level}, 'FontSize', 12)
    set(gcf,'Position',[100 100 1250 1250])
    exportgraphics(gcf,save_dir,'Resolution',500);
end
%% calculate partitions splits
function [data_high, data_low, high_ids, low_ids] ...
    = get_partitions_medium_split(data, participant_order, regression_type, ...
    partition, type_of_effect, weight_erps, weighting_factor)

    function [split_data, p_order] = get_participants(data, all_ids, current_ids)
        cnt = 1;
        split_data = {};
        p_order = {};
        for i=1:numel(data)
            participant = data{i};
            id = all_ids(i);

            if ismember(id{1}, current_ids)
                split_data{cnt} = participant;
                p_order{cnt} = id{1};
                cnt = cnt+1;
            end     
        end      
    end

    function data = weight_erps_based_on_score(data, ranks, type, order, weighting_factor)
        upper_weighting = 1 + weighting_factor;
        lower_weighting = 1 - weighting_factor;
        
        n = size(ranks, 1);
        quartile = floor(n/2);
        if strcmp(type, 'high')
           ids_in_quartile = ranks(1:quartile,2);
        else
           ids_in_quartile = ranks(quartile:end,2);
        end
        
        in_quartile = size(order,2);
        for i=1:in_quartile
            participant_order = order{i};
            participant = data{i};
            
            if ismember(participant_order,ids_in_quartile)
               participant.avg = participant.avg * upper_weighting;
               participant.thin = participant.thin * upper_weighting;
               participant.med = participant.med * upper_weighting;
               participant.thick = participant.thick * upper_weighting;
               participant.weighting = upper_weighting;
               data{i} = participant;
            else
               participant.avg = participant.avg * lower_weighting;
               participant.thin = participant.thin * lower_weighting;
               participant.med = participant.med * lower_weighting;
               participant.thick = participant.thick * lower_weighting;
               participant.weighting = lower_weighting;
               data{i} = participant;
            end
        end      
    end

    scores = return_scores(regression_type, type_of_effect);

    if contains(regression_type, 'median-split')
        ratings = scores.one;
    else
        if partition == 1
            ratings = scores.one;
        elseif partition == 2
            ratings = scores.two;
        elseif partition == 3
            ratings = scores.three;
        end
    end

    if strcmp(regression_type, 'headache-masked')
        ps = ratings(:,1);
        scores = ratings(:,2);
        scores(scores>0)=1;
        scores(scores<0)=-1;
        ratings(:,1) = ps;
        ratings(:,2) = scores;
        low_ids = ratings(ratings(:,2)<0);
        high_ids = ratings(ratings(:,2)>0);
    elseif strcmp(regression_type, 'headache-quartile')
        sorted(:,1) = ratings(:,2);
        sorted(:,2) = ratings(:,1);
        sorted = flipud(sortrows(sorted));
        n = ceil(numel(sorted(:,1))/2);
        high = sorted(1:n/2,:);
        high_ids = high(:,2);
        low = sorted(n+n/2+1:end,:);
        low_ids = low(:,2);
    else
        sorted(:,1) = ratings(:,2);
        sorted(:,2) = ratings(:,1);
        sorted = flipud(sortrows(sorted));
        n = ceil(numel(sorted(:,1))/2);
        high = sorted(1:n,:);
        high_ids = high(:,2);
        low = sorted(n+1:end,:);
        low_ids = low(:,2);
    end
    [data_high, high_order] = get_participants(data, participant_order, high_ids);
    [data_low, low_order] = get_participants(data, participant_order, low_ids); 
    
    if weight_erps == 1
       [data_high] = weight_erps_based_on_score(data_high, high, 'high', high_order, weighting_factor);
       [data_low] = weight_erps_based_on_score(data_low, low, 'low', low_order, weighting_factor);
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
           participant_level.weighting = 1;
           
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

% mainly to handle memory constraints
function ft = cut_data_using_analysis_window(ft, analysis_window)
    time = ft.time*1000;
    [~,start_t] = min(abs(time-analysis_window(1)));
    [~,end_t] = min(abs(time-analysis_window(2)));
    
    fields = fieldnames(ft);
    
    for k=1:numel(fields)
        time_series = ft.(fields{k});
        if numel(size(time_series)) == 2 || numel(size(time_series)) == 1
            filtered_series = time_series(:,start_t:end_t);
            ft.(fields{k}) = filtered_series;
        else
            new_time_series = [];
            [~,~,z] = size(time_series);
            for i=1:z
               trial =  time_series(:,start_t:end_t,i);
               new_time_series(:,:,i) = trial;
            end
            ft.(fields{k}) = new_time_series;
        end          
    end
end

%% applies the wavelett decomposition to the data
function dataset = to_frequency_data(data, save_dir, partition, participant_order, type, fname, frequency_type, channel)
    
     if strcmp(frequency_type, 'fourier')
         cfg = [];
         cfg.channel = 'eeg';
         cfg.method = 'wavelet';
         cfg.width = 3;
         cfg.output = 'fourier';
         cfg.pad = 'nextpow2';
         cfg.foi = 5:30;
         cfg.toi = -0.5:0.002:0.5;  
     elseif strcmp(frequency_type, 'pow')
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.width = 3;
        cfg.foi =   1:30;
        cfg.t_ftimwin = ones(length(cfg.foi),1).*0.25;
        cfg.toi          = -0.5:0.002:0.5;
        cfg.channel      = 'all';
     end

    dataset = {};
    for i=1:numel(data)
        participant = data{i};
        
        disp(strcat('Loading/Processing Participant ', int2str(i)));
        participant_number = participant_order{i};
        med_path = strcat(save_dir, int2str(participant_number), '\', 'partition_', int2str(partition), '_', fname{1});
        thin_path = strcat(save_dir, int2str(participant_number), '\', 'partition_', int2str(partition), '_', fname{2});
        thick_path = strcat(save_dir, int2str(participant_number), '\', 'partition_', int2str(partition), '_', fname{3});
        
        if strcmp(type, 'preprocess')
            med.label = participant.label;
            med.elec = participant.elec;
            med.trial = participant.med;
            med.time = update_with_time_info(med.trial, participant.time);
            med.dimord = 'chan_time';

            thick.label = participant.label;
            thick.elec = participant.elec;
            thick.trial = participant.thick;
            thick.time = update_with_time_info(thick.trial, participant.time);
            thick.dimord = 'chan_time';

            thin.label = participant.label;
            thin.elec = participant.elec;
            thin.trial = participant.thin;
            thin.time = update_with_time_info(thin.trial, participant.time);
            thin.dimord = 'chan_time';

            TFRwave_med = ft_freqanalysis(cfg, med);
            TFRwave_med.info = 'medium';
            save(med_path, 'TFRwave_med', '-v7.3')
            clear TFRwave_med;
            
            TFRwave_thick = ft_freqanalysis(cfg, thick);
            TFRwave_thick.info = 'thick';
            save(thick_path, 'TFRwave_thick', '-v7.3')
            clear TFRwave_thick;

            TFRwave_thin = ft_freqanalysis(cfg, thin);
            TFRwave_thin.info = 'thin';
            save(thin_path, 'TFRwave_thin', '-v7.3')
            clear TFRwave_thin;
        elseif strcmp(type, 'load')
            load(med_path);            
            load(thin_path);
            load(thick_path);
            
            if strcmp(frequency_type, 'pow')
                participant.thin = TFRwave_thin;
                participant.thick = TFRwave_thick;
                participant.med = TFRwave_med; 
            elseif strcmp(frequency_type, 'fourier')
                med = calculate_itc(TFRwave_med, channel);
                thick = calculate_itc(TFRwave_thick, channel);
                thin = calculate_itc(TFRwave_thin, channel);
                participant.thin = thin;
                participant.thick = thick;
                participant.med = med; 
            end
            
            participant.participant_number = participant_number;
            dataset{end+1} = participant;
        end
    end
end

%% update with time information
function trial_info = update_with_time_info(trial, time)
    trial_info = {};
    n = numel(trial);
    for i=1:n
        trial_info{i} = time;
    end
end

%% load the frequency data
function itc_data = calculate_itc(data, required_channel)
    
    itc           = [];
    itc.label     = data.label;
    itc.freq      = data.freq;
    itc.time      = data.time;
    itc.dimord    = 'chan_freq_time';
    
    channel_idx = find(contains(itc.label,required_channel));
    
    F = data.fourierspctrm;
    N = size(F,1);           % number of trials
    F = F(2:end,:,:,:);
    
    % compute inter-trial phase coherence (itpc)
    itc.itpc      = F./abs(F);         % divide by amplitude
    itc.itpc      = sum(itc.itpc,1);   % sum angles
    itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
    itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension
    inter_trial_coherence = itc.itpc(channel_idx, :, :); % get the nth channel
    inter_trial_coherence = squeeze(inter_trial_coherence); % remove 1st dim
    
    itc_data.inter_trial_coherence = inter_trial_coherence;
    itc_data.time = itc.time;
    itc_data.freq = data.freq;
end

%% calculate the arithmetic mean based on weightings
function avg_mtx = calculate_aritmetic_mean(data)
    
    n = size(data,2);
    total_weight = 0;
    matricies = [];
    for i=1:n
        matricies(:,:,i) = data{i}.series;
        total_weight = total_weight + data{i}.weighting;
    end

    sum_mtx = sum(matricies,3);
    avg_mtx = sum_mtx/total_weight;
    
end

%% compute spectrograms of the participant data
 function f_data = compute_spectrogram(data, partition, save_path, channel, regressor, participants, frequency_type)

    if strcmp(frequency_type, 'pow')

        [thin, med, thick] = deal({}, {}, {});
        n_participants = size(data,2);
        
        for i=1:n_participants
            participant = data{i};
            thin{end+1} = participant.thin;
            thick{end+1} = participant.thick;
            med{end+1} = participant.med;
        end

        % save the freq plots
        cfg = [];
        thin_avg = ft_freqgrandaverage(cfg, thin{:});
        thick_avg = ft_freqgrandaverage(cfg, thick{:});
        med_avg = ft_freqgrandaverage(cfg, med{:});

        f_data.thin = thin_avg;
        f_data.thick = thick_avg;
        f_data.med = med_avg;
        f_data.time = data{1}.med.time;
        f_data.freq = data{1}.med.freq;
        

    elseif strcmp(frequency_type, 'fourier')

        [thin, med, thick] = deal([], [], []);
        n_participants = size(data,2);
        
        for i=1:n_participants
            participant = data{i};
            thin(end+1,:,:) = participant.thin.inter_trial_coherence;
            thick(end+1,:,:) = participant.thick.inter_trial_coherence;
            med(end+1,:,:) = participant.med.inter_trial_coherence;
        end

        time = participant.med.time;
        frequencies = participant.med.freq;
        avg_thin = squeeze(mean(thin,'omitnan'));
        avg_thick = squeeze(mean(thick,'omitnan'));
        avg_med = squeeze(mean(med,'omitnan'));
        
        f_data.thin = avg_thin;
        f_data.thick = avg_thick;
        f_data.med = avg_med;
        f_data.time = time;
        f_data.freq = frequencies;
        
    end
 end

 %% compute a barchat of frequencies through time
 function compute_frequency_bar_chart(mtx, time, freq, save_path)
    interval = 0.020; % 20ms
    max_time = 0.300; % 300ms
    
    % theta, alpha, beta, gamma
    bands = [[5,7]; [8,15];[16,31];[32, 100]];
    N = size(bands,1);


    for b=1:N
        f = bands(b,:);

        save_path = strcat(save_path, {'_'}, int2str(f(1)), {'-'}, int2str(f(2)), '.png');
        save_path = save_path{1};
        
        time_of_interest = 0; % start a 0ms
        time_axes = [];
        power_axes = [];
        while time_of_interest <= max_time
            
            % get indices that fall between the desired time intervals
            [~, time_idx] = find(time>=time_of_interest & time<=time_of_interest+interval);
            [~, freq_idx] = find(freq>=f(1) & freq<=f(2));
            
            % get the data for plotting
            subset = mtx(freq_idx,time_idx);
            avg_at_t = mean(subset,'all','omitnan');
            
            time_axes(end+1) = time_of_interest;
            power_axes(end+1) = avg_at_t;
            time_of_interest = time_of_interest + interval;
        end
        
        figure;
        plot(time_axes, power_axes);
        xlabel('Time in 20ms intervals');
        ylabel('ITC value at time T (0-1)')
        t = strcat('ITC values at frequencies:', {' '}, int2str(f(1)), {'-'}, int2str(f(2)));
        t = t{1};
        set(get(gca, 'title'), 'string', t)
        set(gca, 'xtick', time_axes);
        ylim([0,1]);
        grid on
        set(gcf,'Position',[100 100 750 750]);
        exportgraphics(gcf,save_path,'Resolution',500);
        close;
    end
 end
 
 %% hardcoded values from the MUA analysis which returns apriori electrode
 function electrode = return_mua_electrode(factor)
    if strcmp(factor, 'no-factor')
       electrode = 'A22'; 
    elseif strcmp(factor, 'headache')
       electrode = 'A22'; 
    elseif strcmp(factor, 'discomfort')
       electrode = 'A15'; 
    elseif strcmp(factor, 'visual_stress')
       electrode = 'A22';
    end
 end
 
 %% Woody Algorithm Implementation
 % x is a matrix where each row is a participants ERP at electrode X over
 % time
 function out = woody_filtering(x, e_idx)
    %Default parameter values
    
    N = numel(x);
    erps = [];
    for i=1:N
        erps(i,:) = x{i}.avg(e_idx,:);
    end
    
    x = erps';
    x = x(1:257,:);
    
    tol= 0.1;
    max_it=100;
    xcorr_mthd='unbiased';
    
    [N,M]=size(x);
    mx=mean(x,2);
    p=zeros(N,1);
    conv=1;
    run=0;
    sig_x=diag(sqrt(x'*x));
    X=xcorr(mx);
    ref=length(X)/2;

    if(mod(ref,2))
        ref=ceil(ref);
    else
        ref=floor(ref);
    end

    if(nargout>1)
        %In this case we output the lag of the trials as well
        lag_data=zeros(1,M);
    end

    while(conv*(run<max_it))
        disp(run);
        z=zeros(N,1);
        w=ones(N,1);
        for i=1:M

            y=x(:,i);
            xy=xcorr(mx,y,xcorr_mthd);
            
            [~,ind]=max(xy);
            if(ind>ref)
                lag=ref-ind-1;
            else
                lag=ref-ind;
            end
            if(lag>0)
                num=w(lag:end)-1;
                z(1:N-lag+1)=( z(1:N-lag+1).*num + y(lag:end))./w(lag:end);
                w(lag:end)=w(lag:end)+1;
            elseif(lag<0)
                num=w(lag*(-1)+1:end)-1;
                z(lag*(-1)+1:end)=( z(lag*(-1)+1:end).*num + y(1:N+lag) )./w(lag*(-1)+1:end);
                w(lag*(-1)+1:end)=w(lag*(-1)+1:end)+1;
            else
                z=z.*(w-1)./w + y./w;
                w=w+1;
            end
            if(exist('lag_data','var'))
                lag_data(i)=lag;
            end
        end
        
        old_mx=mx;
        mx=z;
        p_old=p;
        p=mx'*x./(sqrt(mx'*mx).*sig_x');
        p=sum(p)./M;
        err=abs(p-p_old);
        if(err<tol)
            conv=0;
        end
        run=run+1;
    end
    out=mx;
 end
 
 %% plot sepctrogram
 function plot_spectrogram(data, save_path, partition, plotting_type, ...
     channel, cat_type)
    participants = size(data,2);
    
    thin = data.thin;
    thick = data.thick;
    med = data.med;

    if strcmp(plotting_type, 'pow')
        cfg = [];
        cfg.baseline = 'yes';
        cfg.baseline     = [-0.5 0];
        cfg.baselinetype = 'db';
        cfg.maskstyle    = 'saturation';
        cfg.xlim = [-0.200,0.500];
        cfg.ylim = [0, 30];
        cfg.zlim = [-3, 3];
        cfg.channel = channel;
        
        % save medium
        title = strcat(cat_type, {' '}, 'Partition:', {' '}, int2str(partition), ',', {' '}, 'Grating:', {' '},...
            'Medium,', {' '}, 'Channel:' ,{' '}, channel);
        title = title{1};
        cfg.title = title;
        figure
        ft_singleplotTFR(cfg, med);
        save_dir = strcat(save_path, '\spectrograms\',  cat_type, '_p', int2str(partition),'_medium_freq.png');
        exportgraphics(gcf,save_dir,'Resolution',500);
        close;

        % save thin
        title = strcat(cat_type, {' '}, 'Partition:', {' '}, int2str(partition), ',', {' '}, 'Grating:', {' '},...
            'Thin,', {' '}, 'Channel:' ,{' '}, channel);
        title = title{1};
        cfg.title = title;
        figure
        ft_singleplotTFR(cfg, thin);
        save_dir = strcat(save_path, '\spectrograms\', cat_type, '_p', int2str(partition),'_thin_freq.png');
        exportgraphics(gcf,save_dir,'Resolution',500);
        close;

        % save thick
        title = strcat(cat_type, {' '}, 'Partition:', {' '}, int2str(partition), ',', {' '}, 'Grating:', {' '},...
            'Thick,', {' '}, 'Channel:' ,{' '}, channel);
        title = title{1};
        cfg.title = title;
        figure
        ft_singleplotTFR(cfg, thick);
        save_dir = strcat(save_path, '\spectrograms\', cat_type, '_p', int2str(partition),'_thick_freq.png');
        exportgraphics(gcf,save_dir,'Resolution',500);
        close;
    elseif strcmp(plotting_type, 'fourier')
        time = data.time;
        freq = data.freq;
        
        t = strcat(cat_type, {' '}, 'Partition:', {' '}, int2str(partition), ',', {' '}, 'Grating:', {' '},...
            'Med,', {' '}, 'Channel:' ,{' '}, channel);
        t = t{1};
        imagesc(time, freq, med);
        axis xy;
        set(get(gca, 'title'), 'string', t)
        save_dir = strcat(save_path, '\itc\', cat_type, {'_'}, 'itc_p', int2str(partition),'_med_freq.png');
        save_dir = save_dir{1};
        xlabel('time -200 to 500ms');
        ylabel('Freq');
        set(gcf,'Position',[100 100 750 750]);
        exportgraphics(gcf,save_dir,'Resolution',500);
        close;


        t = strcat(cat_type, {' '}, 'Partition:', {' '}, int2str(partition), ',', {' '}, 'Grating:', {' '},...
            'Thin,', {' '}, 'Channel:' ,{' '}, channel);
        t = t{1};
        imagesc(time, freq, thin);
        axis xy;
        set(get(gca, 'title'), 'string', t)
        save_dir = strcat(save_path, '\itc\', cat_type, {'_'}, 'itc_p', int2str(partition),'_thin_freq.png');
        save_dir = save_dir{1};
        xlabel('time -200 to 500ms');
        ylabel('Freq');
        set(gcf,'Position',[100 100 750 750]);
        exportgraphics(gcf,save_dir,'Resolution',500);
        close;
        
        t = strcat(cat_type, {' '}, 'Partition:', {' '}, int2str(partition), ',', {' '}, 'Grating:', {' '},...
            'Thick,', {' '}, 'Channel:' ,{' '}, channel);
        t = t{1};
        imagesc(time, freq, thick);
        axis xy;
        set(get(gca, 'title'), 'string', t)
        save_dir = strcat(save_path, '\itc\', cat_type, {'_'}, 'itc_p', int2str(partition),'_thick_freq.png');
        save_dir = save_dir{1};
        xlabel('time -200 to 500ms');
        ylabel('Freq');
        set(gcf,'Position',[100 100 750 750]);
        exportgraphics(gcf,save_dir,'Resolution',500);
        close;
    end
 end
 
 %% aggregate the power data across participants
 function new_data = aggregate_freq_data(data, type)
    
    if strcmp(type, 'pow')

        fields = fieldnames(data);
        N = numel(fields);

        [thin_pwrspec, thick_pwrspec, med_pwrspec] = deal([],[],[]);
        for i= 1:N
           data_i = data.(fields{i}); 
           if i == 1
                example_thin = data_i.thin;
                example_thick = data_i.thick;
                example_med = data_i.med;
           end

           thin_pwr = data_i.thin.powspctrm;
           thick_pwr = data_i.thick.powspctrm;
           med_pwr = data_i.med.powspctrm;

           thin_pwrspec(:,:,:,end+1) = thin_pwr;
           thick_pwrspec(:,:,:,end+1) = thick_pwr;
           med_pwrspec(:,:,:,end+1) = med_pwr;

        end

        avg_thick = mean(thick_pwrspec,4);
        avg_thin = mean(thin_pwrspec,4);
        avg_med = mean(med_pwrspec,4);

        example_thin.powspctrm = avg_thin;
        example_thick.powspctrm = avg_thick;
        example_med.powspctrm = avg_med;

        new_data.thin = example_thin;
        new_data.thick = example_thick;
        new_data.med = example_med;
    elseif strcmp(type, 'fourier')
        new_data.thin = data.thin;
        new_data.thick = data.thick;
        new_data.med = data.med;
    end
 end
 
 %% create power values from each participant
 function average_power_values(data, freq, time, electrode)
    N = numel(data);
    electrode = find(contains(data{1}.label,electrode));
    
    participant_data = [];
    for i=1:N
       participant = data{i}; 
       participant_number = participant.participant_number;
       
       data_time = participant.med.time;
       [~,start_time] = min(abs(data_time-time(1)));
       [~,end_time] = min(abs(data_time-time(2)));
       
       data_freq = participant.med.freq;
       [~,start_freq] = min(abs(data_freq-freq(1)));
       [~,end_freq] = min(abs(data_freq-freq(2)));
       
       
       thin = squeeze(participant.thin.powspctrm(electrode, ...
           start_freq:end_freq, ...
           start_time:end_time));
       avg_thin = mean(thin,'all');
       
       thick = squeeze(participant.thick.powspctrm(electrode, ...
           start_freq:end_freq, ...
           start_time:end_time));
       avg_thick = mean(thick,'all');
       
       med = squeeze(participant.med.powspctrm(electrode, ...
           start_freq:end_freq, ...
           start_time:end_time));
       avg_med = mean(med,'all');

       avg_pgi = avg_med - (avg_thin + avg_thick)/2;
       
       participant_data(i, 1) = avg_pgi;
       participant_data(i,2) = participant_number;
    end
 end