%pipeline for preprocessing in paralell computing to speed up process, Author: Tom Jefferis, built using script by Cihan Dogan

clear all;
restoredefaultpath;
addpath('W:\PhD\MatlabPlugins\fieldtrip-20210906'); %path to fieldtrip
addpath('W:\PhD\MatlabPlugins\spm12') %path to spm
ft_defaults;
cd("W:\PhD\PatternGlareCode");

%%Vars needed ot edit for preprocessing
main_path = 'W:\PhD\PatternGlareData\participants\participant_';
analysis_type = 'mean_intercept';
type_of_analysis = 'time_domain'; % frequency_domain or time_domain
% what onsets you want to include
onsets = [
    2, 3, 4, 5, 6, 7, 8
    ];
number_of_onsets = size(onsets);
number_of_onsets = number_of_onsets(1);
baseline_windows = {[-0.2 0], [2.8 3.0], [3.7 3.9]};
filter_freq = [0.1, 60];
n_participants = 1;

%for each of the baseline windows
for baseline_window = baseline_windows
    %converts from SPM to FT
    for index = 1:n_participants
        %% gets the onsets of interest
        [thin, med, thick, description] = get_onsets(1, analysis_type);
        full_description = strcat(analysis_type, '_', description{1});
        full_description = strcat(type_of_analysis, {'_'}, full_description);
        full_description = full_description{1};

        %% works out where to load the data
        participant_main_path = strcat(main_path, int2str(index));

        if exist(participant_main_path, 'dir')
            participant_main_path = strcat(participant_main_path, '\');
            data_structure = 'spmeeg_P';

            if index < 10
                p = strcat('0', int2str(index));
            else
                p = int2str(index);
            end

            data_structure = strcat(data_structure, p);

            if contains(type_of_analysis, 'frequency_domain')
                data_fname = strcat(data_structure, '_075_80Hz_rejected.dat');
                data_structure = strcat(data_structure, '_075_80Hz_rejected.mat');
                filter_freq = [0.1, 80];
                %baseline_window = [-0.5 0];
            elseif strcmp(type_of_analysis, 'time_domain')
                data_fname = strcat(data_structure, '_075_80Hz.dat');
                data_structure = strcat(data_structure, '_075_80Hz.mat');
                filter_freq = [0.1, 30];
                %baseline_window = [-0.2 0];
            end

            file_main_path = strcat(participant_main_path, data_structure);

            if ~isfile(file_main_path)
                continue;
            end

            cd(participant_main_path);

            %% this function updates the trial information so that you only
            % analyse the conditions of interest
            condition_names = label_data(thin, med, thick, participant_main_path, data_structure, analysis_type);

            %% load and convert from SPM > FieldTrip
            load(file_main_path);
            D.data.fname = strcat(participant_main_path, data_fname);
            spm_eeg = meeg(D);
            raw = spm_eeg.ftraw;

            %% setup the FT preprocessing fns
            % filtering and baselining the data
            cfg = [];
            cfg.demean = 'yes';
            cfg.baselinewindow = baseline_window;

            cfg.bpfilter = 'yes';
            cfg.bpfilttype = 'fir';
            cfg.bpfreq = filter_freq;

            data = ft_preprocessing(cfg, raw);

            format = strcat("ftformat_P", p);
            path = strcat(participant_main_path, format, ".mat");
            parsave(path, raw);
        else
            continue
        end

    end

    clearvars raw

    data_files = {}
    temp_saves = {}

    for index = 1:n_participants
        participant_main_path = strcat(main_path, int2str(index));

        if exist(participant_main_path, 'dir')
            participant_main_path = strcat(participant_main_path, '\');

            if index < 10
                p = strcat('0', int2str(index));
            else
                p = int2str(index);
            end

            format = strcat("ftformat_P", p);
            temp_format = strcat("temp_P", p);
            temp_saves{index} = strcat(participant_main_path, temp_format, ".mat");
            data_files{index} = strcat(participant_main_path, format, ".mat");
        else
            continue
        end

    end

    for index = 1:n_participants

        if exist(data_files{index}, 'file')
            %% setup the FT preprocessing fns
            % filtering and baselining the data
            cfg = [];
            cfg.demean = 'yes';
            cfg.baselinewindow = baseline_window;

            cfg.bpfilter = 'yes';
            cfg.bpfilttype = 'fir';
            cfg.bpfreq = filter_freq;
            raw = load(data_files{index});
            raw = raw.x;
            data = ft_preprocessing(cfg, raw.x);
            parsave(temp_saves{index}, data);
        else
            continue
        end

    end

    clearvars data

    parfor i = 1:n_participants

        if exist(temp_saves{i}, 'file')

            data = load(temp_saves{i});
            % Detect artefacts via thresholding -100:100 uV

            cfg = [];
            cfg.continious = 'no';
            cfg.artfctdef.threshold.min = -100;
            cfg.artfctdef.threshold.max = 100;
            cfg.artfctdef.threshold.channel = get_eeg_channels(data.x);
            cfg.artfctdef.threshold.bpfilter = 'no';

            [~, artifact] = ft_artifact_threshold(cfg, data.x);

            % reject the detected artefacts
            cfg = [];
            cfg.artfctdef.reject = 'complete';
            cfg.artfctdef.zvalue.artifact = artifact;
            postprocessed = ft_rejectartifact(cfg, data);

            parsave(temp_saves{i}, postprocessed)
        else
            continue
        end

    end

    for i = 1:n_participants

        if exist(temp_saves{i}, 'file')
            postprocessed = load(temp_saves{i});
            raw = load(data_files{i});
            raw = raw.x;

            participant_main_path = strcat(main_path, int2str(participant), "\");

            if index < 10
                p = strcat('0', int2str(participant));
            else
                p = int2str(participant);
            end

            data_structure = strcat(participant_main_path, 'spmeeg_P', p, '_075_80Hz.mat');
            file = load(file_main_path);

            % update with the proper trial names after artefact rejection
            postprocessed = label_data_with_trials(raw, postprocessed);
            postprocessed = relabel_conditions(postprocessed, file.D);

            % reject based on count of trials per condition
            reject_participant = reject_particiapnt_based_on_bad_trials(postprocessed, raw);

            if reject_participant == 1
                fprintf(strcat('REJECTED PARTICIPANT...', int2str(participant)));
                continue;
            end

            % get the data ready for FT analysis
            [trial_level, grand_averages] = data_ready_for_analysis(postprocessed, analysis_type);

            % saves the grand average data (easier to load rather than
            % trial level. Also saves trial level if needed.
            save_data(grand_averages, participant_main_path, "mean_intercept_onsets", strcat('_grand-average_baseline_', string(h)))
            save_data(trial_level, participant_main_path, "mean_intercept_onsets", strcat('_trial-level_baseline_pre_off_pl', string(h)))
            delete(temp_saves{i});
            delete(data_files{i});
            fprintf(strcat('PROCESSED PARTICIPANT..', int2str(participant)));

        end

    end

end

function parsave(fname, x)
    save(fname, 'x')
end

%%% Helper Function Wrtitten By Cihan Dogan
function save_data(data, participant_main_path, description, data_type)
    path = strcat(participant_main_path, description, data_type, '.mat');
    save(path, 'data');
end

%% return the desired onsets
function [thin, med, thick, description] = get_onsets(onsets, analysis_type)
    % below is purely a reference so we have all codes
    onsets_thin_REF = {'65411'; '65412'; '65413'; '65414'; '65415'; '65416'; '65417'; '65418'; '65419'};
    onsets_medium_REF = {'65401'; '65402'; '65403'; '65404'; '65405'; '65406'; '65407'; '65408'; '65409'};
    onsets_thick_REF = {'65391'; '65392'; '65393'; '65394'; '65395'; '65396'; '65397'; '65398'; '65399'};

    thin = onsets_thin_REF(onsets);
    med = onsets_medium_REF(onsets);
    thick = onsets_thick_REF(onsets);

    shape = size(onsets);
    number_of_onsets = shape(2);

    description = 'onsets';

    for i = 1:number_of_onsets
        onset = int2str(onsets(i));
        description = strcat(description, {'_'}, onset);
    end

    if strcmp(analysis_type, 'partitions')
        str = '';
        cnt = 1;

        for o = onsets
            str = strcat(str, int2str(o));
            str = strcat(str, '_');
        end

        str = str(1:end - 1);
        description = {strcat('partitioned_onsets', {'_'}, str)};
        description = description{1};
    end

end

%% ensure there is atleast 20% of stimulus type per bucket
function reject_participant = reject_particiapnt_based_on_bad_trials(postprocessed, raw)
    trial_info = raw.trialinfo;
    [original_n_occurence, original_conditions] = hist(trial_info, unique(trial_info));

    pp_trial_info = postprocessed.sampleinfo(:, 3);
    [pp_original_n_occurence, pp_original_conditions] = hist(pp_trial_info, unique(pp_trial_info));

    % make sure we have a minimum number of trials
    if ismember(1, pp_original_n_occurence) || ismember(0, pp_original_n_occurence)
        reject_participant = 1;
        return;
    end

    % lost a condition due to postprocessing
    original_size = numel(original_conditions);
    pp_size = numel(pp_original_conditions);

    if original_size ~= pp_size
        reject_participant = 1;
        return;
    end

    % make sure there are 20% trial per condition
    [~, n_conditions] = size(pp_original_n_occurence);

    for i = 1:n_conditions
        original_condition = original_conditions(i);
        new_condition = pp_original_conditions(i);

        if original_condition == new_condition
            original_count = original_n_occurence(i);
            new_count = pp_original_n_occurence(i);

            if (new_count / original_count) < 0.20
                reject_participant = 1;
                return;
            end

        end

    end

    reject_participant = 0;

end

function [trial_level, grand_averages] = data_ready_for_analysis(postprocessed, data_type)

    postprocessed = remove_electrodes(postprocessed);

    idx_used_for_saving_data = 1;
    trial_names_and_order = postprocessed.trial_order;
    sample_information = postprocessed.sampleinfo;

    % find the condition labels used to match up the data
    if strcmp(data_type, 'partitions')
        p1_thin_idx = find(contains(trial_names_and_order, '_partition_1__thin'));
        p1_thick_idx = find(contains(trial_names_and_order, '_partition_1__thick'));
        p1_med_idx = find(contains(trial_names_and_order, '_partition_1__medium'));
        p2_thin_idx = find(contains(trial_names_and_order, '_partition_2__thin'));
        p2_thick_idx = find(contains(trial_names_and_order, '_partition_2__thick'));
        p2_med_idx = find(contains(trial_names_and_order, '_partition_2__medium'));
        p3_thin_idx = find(contains(trial_names_and_order, '_partition_3__thin'));
        p3_thick_idx = find(contains(trial_names_and_order, '_partition_3__thick'));
        p3_med_idx = find(contains(trial_names_and_order, '_partition_3__medium'));
        % put the trials into the respective buckets
        trials = postprocessed.trial;
        [~, n] = size(trials);

        [p1_thin, p1_thick, p1_med, p2_thin, p2_thick, p2_med, ...
                p3_thin, p3_thick, p3_med] = deal([], [], [], [], [], [], [], [], []);

        for idx = 1:n
            trial = trials{idx};
            condition = sample_information(idx, 3);

            if condition == p1_thin_idx
                p1_thin(:, :, end + 1) = trial;
            elseif condition == p2_thin_idx
                p2_thin(:, :, end + 1) = trial;
            elseif condition == p3_thin_idx
                p3_thin(:, :, end + 1) = trial;
            elseif condition == p1_thick_idx
                p1_thick(:, :, end + 1) = trial;
            elseif condition == p2_thick_idx
                p2_thick(:, :, end + 1) = trial;
            elseif condition == p3_thick_idx
                p3_thick(:, :, end + 1) = trial;
            elseif condition == p1_med_idx
                p1_med(:, :, end + 1) = trial;
            elseif condition == p2_med_idx
                p2_med(:, :, end + 1) = trial;
            elseif condition == p3_med_idx
                p3_med(:, :, end + 1) = trial;
            end

        end

        trial_level.p1_med = convert_to_fieldtrip_format(p1_med);
        trial_level.p2_med = convert_to_fieldtrip_format(p2_med);
        trial_level.p3_med = convert_to_fieldtrip_format(p3_med);
        trial_level.p1_thin = convert_to_fieldtrip_format(p1_thin);
        trial_level.p2_thin = convert_to_fieldtrip_format(p2_thin);
        trial_level.p3_thin = convert_to_fieldtrip_format(p3_thin);
        trial_level.p1_thick = convert_to_fieldtrip_format(p1_thick);
        trial_level.p2_thick = convert_to_fieldtrip_format(p2_thick);
        trial_level.p3_thick = convert_to_fieldtrip_format(p3_thick);
        trial_level.elec = postprocessed.elec;
        trial_level.label = postprocessed.label;
        trial_level.time = postprocessed.time(1, 1);

        % calculate the means
        p1_med = mean(p1_med, 3);
        p2_med = mean(p2_med, 3);
        p3_med = mean(p3_med, 3);
        p1_thin = mean(p1_thin, 3);
        p2_thin = mean(p2_thin, 3);
        p3_thin = mean(p3_thin, 3);
        p1_thick = mean(p1_thick, 3);
        p2_thick = mean(p2_thick, 3);
        p3_thick = mean(p3_thick, 3);

        % setup the data structure for analysis
        grand_averages.p1_pgi = p1_med - (p1_thin + p1_thick) / 2;
        grand_averages.p2_pgi = p2_med - (p2_thin + p2_thick) / 2;
        grand_averages.p3_pgi = p3_med - (p3_thin + p3_thick) / 2;

        grand_averages.p1_med = p1_med;
        grand_averages.p2_med = p2_med;
        grand_averages.p3_med = p3_med;

        grand_averages.p1_thin = p1_thin;
        grand_averages.p2_thin = p2_thin;
        grand_averages.p3_thin = p3_thin;

        grand_averages.p1_thick = p1_thick;
        grand_averages.p2_thick = p2_thick;
        grand_averages.p3_thick = p3_thick;

        grand_averages.trialinfo = [1];
        grand_averages.time = postprocessed.time(1, 1);
        grand_averages.elec = postprocessed.elec;
        grand_averages.dimord = 'chan_time';
        grand_averages.label = postprocessed.label;
    elseif strcmp(data_type, 'mean_intercept')
        thin_idx = find(contains(trial_names_and_order, 'thin'));
        med_idx = find(contains(trial_names_and_order, 'medium'));
        thick_idx = find(contains(trial_names_and_order, 'thick'));

        trials = postprocessed.trial;
        [~, n] = size(trials);

        [thin, medium, thick] = deal([], [], []);

        for idx = 1:n
            trial = trials{idx};
            condition = sample_information(idx, 3);

            if condition == thin_idx
                thin(:, :, end + 1) = trial;
            elseif condition == med_idx
                medium(:, :, end + 1) = trial;
            elseif condition == thick_idx
                thick(:, :, end + 1) = trial;
            end

        end

        trial_level.thin = convert_to_fieldtrip_format(thin);
        trial_level.med = convert_to_fieldtrip_format(medium);
        trial_level.thick = convert_to_fieldtrip_format(thick);
        trial_level.elec = postprocessed.elec;
        trial_level.time = postprocessed.time(1, 1);
        trial_level.label = postprocessed.label;

        % calculate means
        thin = mean(thin, 3);
        thick = mean(thick, 3);
        medium = mean(medium, 3);

        grand_averages.thin = thin;
        grand_averages.thick = thick;
        grand_averages.med = medium;
        grand_averages.trialinfo = [1];
        grand_averages.time = postprocessed.time(1, 1);
        grand_averages.elec = postprocessed.elec;
        grand_averages.dimord = 'chan_time';
        grand_averages.label = postprocessed.label;
    end

end

%% rename the trials with the correct name
function postprocessed = relabel_conditions(postprocessed, D)
    [~, n_trials] = size(D.trials);

    label_names = {};

    cnt = 1;

    for i = 1:n_trials
        trial = D.trials(i).label;

        if ~any(strcmp(label_names, trial))
            label_names{cnt} = trial;
            cnt = cnt + 1;
        end

    end

    postprocessed.trial_order = label_names;
end

%% converts the trials to a fieldtrip readable format
function new_trials = convert_to_fieldtrip_format(trials)
    new_trials = {};

    [~, ~, samples] = size(trials);

    for k = 1:samples
        t = trials(:, :, k);
        new_trials{k} = t;
    end

end

%% updates the EEG data with the onsets we are interested in analysing
function condition_names = label_data(thin, medium, thick, path, fname, analysis_type)
    factor_name = '';
    file = strcat(path, fname);
    load(file); % loads the D object
    D.fname = fname;
    D.path = path;
    n_trials = size(D.trials);
    n_trials = n_trials(2);
    count = 1;
    condition_names = {};

    if ~strcmp(analysis_type, 'partitions')

        for onset = 1:n_trials

            events = D.trials(onset).events;
                [~, rows] = size(events);

                for i = 1:rows
                    condition = events(i).binlabel;

                    if ~strcmp(condition, '""')
                        condition_found = 1;
                        break
                    end

                end

                if sum(contains(condition, thin))
                    condition = strcat(factor_name, '_thin');
                elseif sum(contains(condition, medium))
                    condition = strcat(factor_name, '_medium');
                elseif sum(contains(condition, thick))
                    condition = strcat(factor_name, '_thick');
                else
                    condition = 'N/A';
                end

                condition_names{onset} = condition;
                D.trials(onset).label = condition;
                count = count + 1;
            end

        else
            partition_number = 1;
            max_epoch = D.trials(n_trials).events.epoch;

            for onset = 1:n_trials

                events = D.trials(onset).events;
                    current_epoch = D.trials(onset).events.epoch;

                    [~, rows] = size(events);

                    for i = 1:rows
                        condition = events(i).binlabel;

                        if ~strcmp(condition, '""')
                            condition_found = 1;
                            break
                        end

                    end

                    if ~condition_found == 1
                        error('Condition not found...');
                    end

                    if current_epoch <= (max_epoch / 3)
                        partition_number = 1;
                    elseif (current_epoch > (max_epoch / 3)) && (current_epoch <= (max_epoch / 3) * 2)
                        partition_number = 2;
                    elseif (current_epoch > (max_epoch / 3) * 2) && (current_epoch <= max_epoch)
                        partition_number = 3;
                    end

                    description = strcat('partition_', int2str(partition_number));
                    description = strcat(description, '_');

                    if sum(contains(condition, thin))
                        condition = strcat(description, '_thin');
                        condition = strcat('_', condition);
                        condition = strcat(factor_name, condition);
                    elseif sum(contains(condition, medium))
                        condition = strcat(description, '_medium');
                        condition = strcat('_', condition);
                        condition = strcat(factor_name, condition);
                    elseif sum(contains(condition, thick))
                        condition = strcat(description, '_thick');
                        condition = strcat('_', condition);
                        condition = strcat(factor_name, condition);
                    else
                        condition = 'N/A';
                    end

                    condition_names{onset} = condition;
                    D.trials(onset).label = condition;
                end

            end

            condition_names = unique(cellfun(@num2str, condition_names, 'uni', 0));
            condition_names(ismember(condition_names, 'N/A')) = [];
            save(file, 'D')
        end

        %% update the samples with trial info
        function postprocessed = label_data_with_trials(raw, postprocessed)
            original_info = raw.sampleinfo;
            original_info(:, 3) = raw.trialinfo';
            new_info = postprocessed.sampleinfo;

            [row, ~] = size(new_info);

            for i = 1:row
                start_sample = new_info(i, 1);
                idx = find(original_info(:, 1) == start_sample);
                original_label = original_info(idx, 3);
                new_info(i, 3) = original_label;
            end

            postprocessed.sampleinfo = new_info;
        end
