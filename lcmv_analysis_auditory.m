%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This scripts perform source reconstruction beamformer analysis in the
% gamma band using an LCMV beamformer. Data is from my alien paradigm
% (auditory clicktrain).
% Individual trial-based statistics are then computed using cluster-based 
% permutation tests based on the Montercarlo method 
% (Maris & Oostenveld, 2007)

% Written by Robert Seymour (ABC) - March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify subject initials

noise_level = {'noise_1','noise_2','noise_3','noise_4'}; % 4 levels of noisiness

maxfilter_setting = {'quat','sss','quat_tsss','mc','mc_tsss'}; % maxfilter parameters

%% Loop for all subjects
for i=1:length(noise_level)
    for j = 1:length(maxfilter_setting)
        disp(['Performing source localisation using the audiotry data from ' noise_level{i} ' with Maxfilter parameters ' maxfilter_setting{j}])
        
        %% Load the required variables you computed earlier
        load(sprintf('D:\\RS_MC_analysis\\data_clean\\%s_%s_auditory_data_clean.mat',noise_level{i},maxfilter_setting{j}))
        
        if strcmp(noise_level{i},'noise_1')
            load('D:\\RS_MC_analysis\\noise_1_forward\\mri_realigned.mat');
            load('D:\\RS_MC_analysis\\noise_1_forward\\sens.mat')
            load('D:\\RS_MC_analysis\\noise_1_forward\\seg.mat')
        elseif strcmp(noise_level{i},'noise_2')
            load('D:\\RS_MC_analysis\\noise_2_forward\\mri_realigned.mat');
            load('D:\\RS_MC_analysis\\noise_2_forward\\sens.mat')
            load('D:\\RS_MC_analysis\\noise_2_forward\\seg.mat')
        elseif strcmp(noise_level{i},'noise_3')
            load('D:\\RS_MC_analysis\\noise_3_forward\\mri_realigned.mat');
            load('D:\\RS_MC_analysis\\noise_3_forward\\sens.mat')
            load('D:\\RS_MC_analysis\\noise_3_forward\\seg.mat')
        elseif strcmp(noise_level{i},'noise_4')
            load('D:\\RS_MC_analysis\\noise_4_forward\\mri_realigned.mat');
            load('D:\\RS_MC_analysis\\noise_4_forward\\sens.mat')
            load('D:\\RS_MC_analysis\\noise_4_forward\\seg.mat')
        end
        
        
        cd('D:\RS_MC_analysis\auditory_source')
        
        chans_included = {'MEGGRAD', '-MEG0322', '-MEG2542','-MEG0111','-MEG0532'};
        
        %% Make sure we have 64 components
        numcomponent = 64;
        disp(sprintf('\n Reducing the data to %d components \n',numcomponent));
        
        cfg = [];
        cfg.method = 'pca';
        cfg.updatesens = 'yes';
        cfg.channel = chans_included;
        comp = ft_componentanalysis(cfg, data_auditory);
        
        cfg = [];
        cfg.updatesens = 'yes';
        cfg.component = comp.label(numcomponent:end);
        data_auditory = ft_rejectcomponent(cfg, comp);
        
        %% Downsample - necessary step if you want your analysis to run this year!
        cfg = [];
        cfg.resamplefs = 250;
        cfg.detrend = 'no';
        data_auditory_200 = ft_resampledata(cfg,data_auditory);
        
        %% BP Filter & Select Gradiometers
        cfg = [];
        cfg.channel = chans_included; %performing the analysis only on the gradiometers
        cfg.bpfilter = 'yes'
        cfg.bpfreq = [35 45];    %band-pass filter in the required range
        data_filtered = ft_preprocessing(cfg,data_auditory_200)
        
        %% Create leadfields in subject's brain warped to MNI space
        %Load template sourcemodel
        load('D:\fieldtrip-20160208\template\sourcemodel\standard_sourcemodel3d8mm.mat');
        template_grid = sourcemodel;
        template_grid = ft_convert_units(template_grid,'mm')
        clear sourcemodel;
        
        % construct the volume conductor model (i.e. head model) for each subject
        cfg        = [];
        cfg.method = 'singleshell';
        headmodel  = ft_prepare_headmodel(cfg, seg);
        
        % create the subject specific grid, using the template grid that has just been created
        cfg                = [];
        cfg.grid.warpmni   = 'yes';
        cfg.grid.template  = template_grid;
        cfg.grid.nonlinear = 'yes'; % use non-linear normalization
        cfg.mri            = mri_realigned;
        cfg.grid.unit      ='mm';
        %cfg.inwardshift = '1.5';
        grid               = ft_prepare_sourcemodel(cfg);
        
        % create leadfield
        cfg = [];
        cfg.channel = chans_included;
        cfg.grid = grid;
        cfg.vol = headmodel;
        cfg.grad = sens;
        cfg.normalize = 'yes';
        
        lf = ft_prepare_leadfield(cfg)
        
        %% Here we redefine trials based on the time-points of interest.
        % Make sure the timepoints are of equivalent length
        cfg = [];
        cfg.toilim = [-1.5 -0.3];
        datapre = ft_redefinetrial(cfg, data_filtered);
        cfg.toilim = [0.3 1.5];
        datapost = ft_redefinetrial(cfg, data_filtered);
        
        % Here we are keeping all parts of the trial for your covariance matrix
        cfg = [];
        cfg.keeptrials = 'yes';
        cfg.channel = chans_included;
        cfg.covariance = 'yes';
        cfg.covariancewindow = [-1.5 1.5]
        avg = ft_timelockanalysis(cfg,data_filtered);
        
        % Time lock analysis for datapre and datapost period
        cfg = [];
        cfg.keeptrials = 'yes';
        cfg.channel = chans_included;
        cfg.covariance='yes';
        avgpre = ft_timelockanalysis(cfg,datapre);
        avgpst = ft_timelockanalysis(cfg,datapost);
        
        %% Do ya beamforming
        % Source reconstruction for the whole trial
        cfg=[];
        cfg.keeptrials = 'yes';
        cfg.channel = chans_included;
        cfg.method='lcmv';
        cfg.grid=lf;
        cfg.vol=headmodel;
        cfg.lcmv.keepfilter='yes';
        sourceavg=ft_sourceanalysis(cfg, avg);
        
        % Now do beamforming for the two time points separately using the same spatial
        % filter computed from the whole trial
        cfg = [];
        cfg.rawtrial = 'yes';
        cfg.channel = chans_included;
        cfg.method='lcmv';
        cfg.grid=lf;
        cfg.grid.filter=sourceavg.avg.filter; %uses the grid from the whole trial average
        cfg.vol=headmodel;
        %Pre-grating
        sourcepreS1 = ft_sourceanalysis(cfg, avgpre);
        %Post-grating
        sourcepstS1=ft_sourceanalysis(cfg, avgpst);
        
        % Make sure your field positions match the template grid
        sourcepreS1.pos=template_grid.pos; % right(?)
        sourcepstS1.pos=template_grid.pos; % right(?)
        
        save([noise_level{i} '_' maxfilter_setting{j} '_auditory_sourcepre'],'sourcepreS1','-v7.3');
        save([noise_level{i} '_' maxfilter_setting{j} '_auditory_sourcepost'],'sourcepstS1','-v7.3');
        
        %Plot the difference - not necessary but useful for illustration purposes
        cfg = [];
        cfg.parameter = 'pow';
        cfg.operation = '((x1-x2)./x2)*100';
        sourceR=ft_math(cfg,sourcepstS1,sourcepreS1);
        sourceR.pow = mean(sourceR.pow,2); %collapse across trials for illustration
        
        mri = ft_read_mri('D:\fieldtrip-20160208\template\anatomy\single_subj_T1.nii');
        
        cfg              = [];
        cfg.voxelcoord   = 'no';
        cfg.parameter    = 'pow';
        cfg.interpmethod = 'nearest';
        sourceI  = ft_sourceinterpolate(cfg, sourceR, mri);
        
        %% Load atlas and create a binary mask in visual cortex
        atlas = ft_read_atlas('D:\fieldtrip-20160208\template\atlas\aal\ROI_MNI_V4.nii');
        atlas = ft_convert_units(atlas,'mm');% assure that atlas and template_grid are expressed in the %same units
        
        cfg              = [];
        cfg.voxelcoord   = 'no';
        cfg.parameter    = 'pow';
        cfg.interpmethod = 'nearest';
        parcel  = ft_sourceinterpolate(cfg, sourceR, atlas);
        
        dummy=atlas;
        for k=1:length(parcel.pow)
            dummy.tissue(find(dummy.tissue==k))=parcel.pow(k);
        end;
        
        sourceI.parcel=dummy.tissue;
        sourceI.coordsys = 'mni';
        
        cfg=[];
        cfg.method = 'ortho';
        cfg.funparameter = 'pow';
        cfg.funcolormap    = 'jet';
        cfg.atlas = atlas;
        ft_sourceplot(cfg,sourceI);
        saveas(gcf,([noise_level{i} ' ' maxfilter_setting{j} '.png']))
        %clear vol sens seg mri_realigned sourceR sourceI sourcepreS1 sourcepstS1 mri grid lf
        clc
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATISTICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noise_level = {'noise_1','noise_2','noise_3','noise_4'}; % 4 levels of noisiness

maxfilter_setting = {'tsss_mc'}; % maxfilter parameters

for i=1:length(noise_level)
    for j = 1:length(maxfilter_setting)
        disp(['Performing statistical analysis using the auditory data from ' noise_level{i} ' with Maxfilter parameters ' maxfilter_setting{j}])
        
        %% Load the required variables you computed earlier
        load(sprintf('D:\\RS_MC_analysis\\auditory_source\\%s_%s_auditory_sourcepost.mat',noise_level{i},maxfilter_setting{j}))
        load(sprintf('D:\\RS_MC_analysis\\auditory_source\\%s_%s_auditory_sourcepre.mat',noise_level{i},maxfilter_setting{j}))
        
        
        cd('D:\RS_MC_analysis\auditory_source\stats')
        
        cfg=[];
        cfg.parameter   = 'pow';
        cfg.correctm    = 'cluster';
        cfg.computecritval = 'yes';
        cfg.dim         = sourcepreS1.dim;
        cfg.method      = 'montecarlo';
        cfg.statistic   = 'ft_statfun_depsamplesT';
        
        cfg.numrandomization = 1000;
        cfg.alpha       = 0.05; % Set alpha level
        cfg.clusteralpha = 0.05;
        cfg.tail        = 0;    % Two sided testing
        
        % Design Matrix
        ntrials=numel(sourcepreS1.trial);
        design = zeros(2,2*ntrials);
        design(1,1:ntrials) = 1;
        design(1,ntrials+1:2*ntrials) = 2;
        design(2,1:ntrials) = [1:ntrials];
        design(2,ntrials+1:2*ntrials) = [1:ntrials];
        
        cfg.design = design;
        cfg.uvar        = 2; % row of design matrix that contains unit variable (in this case: subjects)
        cfg.ivar        = 1; % row of design matrix that contains independent variable (the conditions)
        
        stat = ft_sourcestatistics(cfg,sourcepstS1, sourcepreS1);
        
        save([noise_level{i} '_' maxfilter_setting{j} '_stat'],'stat','-v7.3');
        
        % Show raw source level statistics
        % cfg               = [];
        % cfg.method        = 'ortho';
        % cfg.funparameter  = 'stat';
        % cfg.location = 'max';
        % %cfg.maskparameter = 'mask';%turn on to show mask
        % cfg.funcolormap = 'jet';
        % ft_sourceplot(cfg,stat);
        
        %% Interpolate onto SPM T1 brain
        mri = ft_read_mri('D:\fieldtrip-20160208\template\anatomy\single_subj_T1.nii');
        
        cfg              = [];
        cfg.voxelcoord   = 'no';
        cfg.parameter    = 'stat';
        cfg.interpmethod = 'nearest';
        statint  = ft_sourceinterpolate(cfg, stat, mri); %your FT variable corresponding to the subject specific nifti
        % cfg.parameter    = 'mask';
        % maskint  = ft_sourceinterpolate(cfg, stat,mri);
        % statint.mask = maskint.mask;
        
        %% Plot interpolated results
        % cfg               = [];
        % cfg.funparameter  = 'stat';
        % cfg.maskparameter = 'mask';
        % cfg.funcolorlim = 'maxabs';
        % cfg.location = 'max';
        % cfg.funcolormap = 'jet';
        % ft_sourceplot(cfg,statint);
        
        %% Export to nifti
        % Run this if you want to export the clusters rather than the raw stats
        % statint.stat(isnan(statint.stat)) = 0;
        % statint.stat = (statint.stat(:).*statint.mask(:)) %masks
        
        % Use ft_sourcewrite to export to nifti
        cfg = [];
        cfg.filetype = 'nifti';
        cfg.filename = [noise_level{i} '_' maxfilter_setting{j} '_stat_nifti'];
        cfg.parameter = 'stat';
        ft_sourcewrite(cfg,statint);
        
        %% Export to connectome workbench
        
        system(['D:\Software\workbench\bin_windows64\wb_command -volume-to-surface-mapping D:\RS_MC_analysis\auditory_source\stats\' noise_level{i} '_' maxfilter_setting{j} '_stat_nifti.nii ' 'D:\Software\workbench\bin_windows64\Conte69_atlas-v2.LR.32k_fs_LR.wb\32k_ConteAtlas_v2\Conte69.L.midthickness.32k_fs_LR.surf.gii D:\RS_MC_analysis\auditory_source\stats\gifti\' noise_level{i} '_' maxfilter_setting{j} '_LEFT.shape.gii -trilinear'])
        system(['D:\Software\workbench\bin_windows64\wb_command -volume-to-surface-mapping D:\RS_MC_analysis\auditory_source\stats\' noise_level{i} '_' maxfilter_setting{j} '_stat_nifti.nii ' 'D:\Software\workbench\bin_windows64\Conte69_atlas-v2.LR.32k_fs_LR.wb\32k_ConteAtlas_v2\Conte69.R.midthickness.32k_fs_LR.surf.gii D:\RS_MC_analysis\auditory_source\stats\gifti\' noise_level{i} '_' maxfilter_setting{j} '_RIGHT.shape.gii -trilinear'])
    end
end

%% Analyse statistics
noise_level = {'noise_1','noise_3','noise_2','noise_4'}; % 4 levels of noisiness
maxfilter_setting = {'quat','sss','quat_tsss','mc','tsss_mc'}; % maxfilter parameters

stat_out = [];

for i=1:length(noise_level)
    for j = 1:length(maxfilter_setting)
        disp(['Getting juicy stats from ' noise_level{i} ' with Maxfilter parameters ' maxfilter_setting{j}])
        
        %% Load the required variables you computed earlier
        load(sprintf('D:\\RS_MC_analysis\\auditory_source\\stats\\%s_%s_stat.mat',noise_level{i},maxfilter_setting{j}))
        
        stat_out(i,j) = stat.posclusters(1).prob;
    end
end
