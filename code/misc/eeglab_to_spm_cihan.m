n_participants = 40;
main_path = 'C:\ProgramFiles\PhD\participant_';
clear matlabbatch

for participant = 1:40
    %% get the correct file path
    clear matlabbatch
    disp(strcat('Procesisng participant..',int2str(participant)));
    participant_main_path = strcat(main_path, int2str(participant));    
    
    
    
    if exist(participant_main_path, 'dir')
        if participant < 10
            p = strcat('P0', int2str(participant));
        else
            p = 'P';
            p = strcat(p,int2str(participant));
        end
        
        cd(participant_main_path);
        
        p = strcat(p, '_075_80Hz.set');
        participant_file_path = strcat(participant_main_path, '\');
        participant_file_path = strcat(participant_file_path, p);
        
        %% convert to SPM ready file
        D = spm_eeg_convert(participant_file_path);
        
        %% clean up the labels
        label = D.chanlabels;
        for i = 1:numel(label)
            lbl = label{i};
            lbl(lbl=='"')=[];
            k = strfind(lbl, '(');
            if ~isempty(k)
                lbl = lbl(1:(k-1));
            end
            label{i} = lbl;
        end
        
        D = chanlabels(D, ':', label);
        
        if isfield(D, 'origchantypes')
            D=D.rmfield('origchantypes');
        end
        
        
        S = [];
        S.D =D;
        S.task = 'defaulttype';
        S.save = 1;
        D = spm_eeg_prep(S);
        
        S = [];
        S.D = D;
        S.task = 'defaulteegsens';
        S.save = 1;
        D = spm_eeg_prep(S);
        
        S = [];
        S.task = 'project3D';
        S.modality = 'EEG';
        S.updatehistory = 0;
        S.D = D;
        
        D = spm_eeg_prep(S);
        
        D = timeonset(D, -0.2);
        
        save(D);
        disp('COMPLETED PARTICIPANT');
        disp(participant);
        
    end
end