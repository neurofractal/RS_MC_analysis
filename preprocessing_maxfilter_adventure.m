%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a Matlab script performs SIMPLE preprocessig for the Maxfilter 
% investigation data. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise_level = {'noise_1','noise_2','noise_3','noise_4'}; % 4 levels of noisiness

maxfilter_setting = {'quat','sss','quat_tsss','mc','mc_tsss'}; % maxfilter parameters

%% Start the loop
for i = 1:length(noise_level)
    for j = 1:length(maxfilter_setting)
        %% Display update
        disp(['Processing the data from ' noise_level{i} ' with Maxfilter parameters ' maxfilter_setting{j}])
        
        %% Get name of rawfile
        if strcmp(noise_level{i},'noise_1')
            filenm = 'rs_asd_rs';
        elseif strcmp(noise_level{i},'noise_2')
            filenm = 'rs_asd_rs_noisy_1';
        elseif strcmp(noise_level{i},'noise_3')
            filenm = 'rs_asd_rs_noisy_2';
        elseif strcmp(noise_level{i},'noise_4')
            filenm = 'rs_asd_rs_noisy_3';
        end
        
        rawfile = (['D:\RS_MC_analysis\' filenm '_aliens_' maxfilter_setting{j} '.fif'])
        
        %% Epoching & Filtering
        % Epoch the whole dataset into one continous dataset and apply
        % the appropriate filters
        cfg = [];
        cfg.headerfile = rawfile;
        cfg.datafile = rawfile;
        cfg.channel = 'MEG';
        cfg.trialdef.triallength = Inf;
        cfg.trialdef.ntrials = 1;
        cfg = ft_definetrial(cfg);
        
        cfg.continuous = 'yes';
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [0.5 250];
        cfg.channel = 'MEG';
        cfg.dftfilter = 'yes';
        cfg.dftfreq = [50];
        alldata = ft_preprocessing(cfg);
        
        % Deal with 50Hz line noise
        cfg = [];
        cfg.bsfilter = 'yes';
        cfg.bsfreq = [49.5 50.5];
        alldata = ft_preprocessing(cfg,alldata);
        
        % Deal with 100Hz line noise
        cfg = [];
        cfg.bsfilter = 'yes';
        cfg.bsfreq = [99.5 100.5];
        alldata = ft_preprocessing(cfg,alldata);
        
        % Epoch your filtered data based on a specific trigger
        cfg = [];
        cfg.headerfile = rawfile;
        cfg.datafile = rawfile;
        cfg.channel = 'MEG';
        cfg.trialdef.eventtype = 'Trigger';
        cfg.trialdef.eventvalue = [16 32];
        cfg.trialdef.prestim = 2.0;         % pre-stimulus interval
        cfg.trialdef.poststim = 2.0;        % post-stimulus interval
        cfg = ft_definetrial(cfg);
        data = ft_redefinetrial(cfg,alldata); %redefines the filtered data
        
        %% Detrend and demean each trial
        cfg = [];
        cfg.demean = 'yes';
        cfg.detrend = 'yes';
        data = ft_preprocessing(cfg,data);
        
        %% Remove Deviant Trials
        if strcmp(noise_level{i},'noise_1')
            trial2remove = 43;
        elseif strcmp(noise_level{i},'noise_2')
            trial2remove = 75;
        elseif strcmp(noise_level{i},'noise_3')
            trial2remove = [];
        elseif strcmp(noise_level{i},'noise_4')
            trial2remove = [14,77,81,85,89,112];
        end
        
        trial_list = [1:1:128];
        trial_list(trial2remove) = [];
        
        cfg = [];
        cfg.trials = trial_list;
        data = ft_redefinetrial(cfg,data);
        clear trial_list trial2remove
        
        %% Redefine into visual and auditory
        indx_visual = find(data.trialinfo == 16);
        indx_auditory = find(data.trialinfo == 32);
        
        cfg = [];
        cfg.trials = indx_visual;
        data_visual = ft_redefinetrial(cfg,data);
        cfg.trials = indx_auditory;
        data_auditory = ft_redefinetrial(cfg,data);
        
        %% SAVE
        disp(['Saving the clean data from ' noise_level{i} ' with Maxfilter parameters ' maxfilter_setting{j}])
        
        cd('D:\RS_MC_analysis\data_clean')
        save([noise_level{i} '_' maxfilter_setting{j} '_visual_data_clean'],'data_visual','-v7.3');
        save([noise_level{i} '_' maxfilter_setting{j} '_auditory_data_clean'],'data_auditory','-v7.3');
        clear data_visual data_auditory data alldata rawfile
    end
end

