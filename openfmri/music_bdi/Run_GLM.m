data_dir = 'C:\original_scan_data\music_mdd\ds171_R1.0.0';

% Find subjects
subject_dirs = spm_select('FPList',data_dir,'dir','^sub-.*');
subject_dirs = cellstr(subject_dirs);

% Output directory
base_dir = pwd;

n = length(subject_dirs);

TR = 3; % Seconds

%% Build onsets (concatenate sessions)

% Original design
TONES = 1;
RESPONSE = 2;
NVE_MUSIC = 3;
PVE_MUSIC = 4;
NVE_NONMUSIC = 5;
PVE_NONMUSIC = 6;
names1 = {'tones','responses','nve-music','pve-music','nve-nonmusic','pve-nonmusic'};

% Simplified design
TASK     = 1;
RESP     = 2;
names2 = {'task','response'};

for i = 1:n
    func_dir = fullfile(subject_dirs{i},'func');
    
    onsets    = cell(1,length(names1));
    durations = cell(1,length(names1));
    
    t = 0; % Time offset for run (seconds)
    
    for r = 1:5        
        % Load events file for this run (onsets,durations,names)
        filename = sprintf('^sub-.*_run-%d_events.tsv$',r);
        P        = spm_select('FPList',func_dir,filename);
        events   = import_events_file(P);
        
        % Identify trials of each type
        is_tone         = strcmp(events(:,3),'tones');
        is_response     = strcmp(events(:,3),'response');
        is_nve_music    = strcmp(events(:,3),'negative_music');
        is_pve_music    = strcmp(events(:,3),'positive_music');
        is_nve_nonmusic = strcmp(events(:,3),'negative_nonmusic');
        is_pve_nonmusic = strcmp(events(:,3),'positive_nonmusic');
        
        % Extract onsets and durations
        ons = cell2mat(events(:,1));
        dur = cell2mat(events(:,2));
        
        onsets{TONES}    = [onsets{TONES};    ons(is_tone) + t];
        durations{TONES} = [durations{TONES}; dur(is_tone)];

        onsets{RESPONSE}    = [onsets{RESPONSE};    ons(is_response) + t];
        durations{RESPONSE} = [durations{RESPONSE}; dur(is_response)];
        
        onsets{NVE_MUSIC}    = [onsets{NVE_MUSIC};    ons(is_nve_music) + t];
        durations{NVE_MUSIC} = [durations{NVE_MUSIC}; dur(is_nve_music)];
        
        onsets{PVE_MUSIC}    = [onsets{PVE_MUSIC};    ons(is_pve_music) + t];
        durations{PVE_MUSIC} = [durations{PVE_MUSIC}; dur(is_pve_music)];
        
        onsets{NVE_NONMUSIC}    = [onsets{NVE_NONMUSIC};    ons(is_nve_nonmusic) + t];
        durations{NVE_NONMUSIC} = [durations{NVE_NONMUSIC}; dur(is_nve_nonmusic)];
        
        onsets{PVE_NONMUSIC}    = [onsets{PVE_NONMUSIC};    ons(is_pve_nonmusic) + t];
        durations{PVE_NONMUSIC} = [durations{PVE_NONMUSIC}; dur(is_pve_nonmusic)];
        
        % Increase time offset by length of run
        filename = sprintf('^wr.*_run-%d_bold.nii$',r);
        P  = spm_select('ExtFPList',func_dir,filename);
        nv = size(P,1);
        t  = t + (nv * TR); 
    end
    
    % Count trials of each type
    nt = [length(onsets{NVE_MUSIC}); 
          length(onsets{PVE_MUSIC}); 
          length(onsets{NVE_NONMUSIC}); 
          length(onsets{PVE_NONMUSIC})];
    
    % Simplify design: task, is_pve, is_music, response
    onsets2    = cell(1,length(names2));
    durations2 = cell(1,length(names2));
    pmod       = [];
    
    onsets2{TASK}    = [onsets{NVE_MUSIC}; onsets{PVE_MUSIC}; onsets{NVE_NONMUSIC}; onsets{PVE_NONMUSIC}];
    durations2{TASK} = [durations{NVE_MUSIC}; durations{PVE_MUSIC}; durations{NVE_NONMUSIC}; durations{PVE_NONMUSIC}];
    
    onsets2{RESP}    = onsets{RESPONSE};
    durations2{RESP} = durations{RESPONSE};
    
    % Parametric regressors (main effects)
    pmod(1).name{1}  = 'is-pve';
    pmod(1).param{1} = [repmat(-1,1,nt(1)) ...
                        repmat( 1,1,nt(2)) ...
                        repmat(-1,1,nt(3)) ...
                        repmat( 1,1,nt(4))]; %#ok<REPMAT>
    pmod(1).poly{1}  = 1;
    
    pmod(1).name{2}  = 'is-music';
    pmod(1).param{2} = [repmat( 1,1,nt(1)) ...
                        repmat( 1,1,nt(2)) ...
                        repmat(-1,1,nt(3)) ...
                        repmat(-1,1,nt(4))]; %#ok<REPMAT>
    pmod(1).poly{2}  = 1;
    
    % Mean correct main effects
    pmod(1).param{1} = pmod(1).param{1} - mean(pmod(1).param{1});
    pmod(1).param{2} = pmod(1).param{2} - mean(pmod(1).param{2});
    
    % Interaction
    pmod(1).name{3}  = 'interaction';
    pmod(1).param{3} = (pmod(1).param{1} .* pmod(1).param{2});
    pmod(1).poly{3}  = 1;    
    
    pmod(1).orth = 0;
    
    % Save
    names     = names2;
    onsets    = onsets2;
    durations = durations2;
    
    [~,subject_name] = fileparts(subject_dirs{i});
    
    out_dir   = fullfile('onsets');
    if ~exist(out_dir,'file')
        mkdir(out_dir);
    end
    
    save(fullfile(out_dir,[subject_name '_onsets.mat']),...
        'names','onsets','durations','pmod');
end
%% Specify first level GLM
for i = 1:n
    fprintf('Subject %d\n',i);
    
    % Load data
    anat_dir = fullfile(subject_dirs{i},'anat');
    func_dir = fullfile(subject_dirs{i},'func');    
    [~,subject_name] = fileparts(subject_dirs{i});
    
    % Make output dir
    glm_dir = fullfile(base_dir,'GLM',subject_name);
    if ~exist(glm_dir,'file')
        mkdir(glm_dir);
    end
    
    % Count the functionals in each run
    nv = zeros(1,5);
    for r = 1:5
        pattern = sprintf('^swr.*_run-%d_bold.nii$',r);
        nv(r) = size(spm_select('ExtFPList',func_dir,pattern),1);
    end
    
    % Collate functionals
    P = [cellstr(spm_select('ExtFPList',func_dir,'^swr.*_run-1_bold.nii$'));
         cellstr(spm_select('ExtFPList',func_dir,'^swr.*_run-2_bold.nii$'));
         cellstr(spm_select('ExtFPList',func_dir,'^swr.*_run-3_bold.nii$'));
         cellstr(spm_select('ExtFPList',func_dir,'^swr.*_run-4_bold.nii$'));
         cellstr(spm_select('ExtFPList',func_dir,'^swr.*_run-5_bold.nii$'))];
     
    % Create concatenated movement file
    R = [];
    for r = 1:5
        % Movement
        pattern = sprintf('^rp_.*_run-%d_bold.txt$',r);
        rp = importdata(spm_select('FPList',func_dir,pattern));        
        R = [R; rp];
    end
    
    % Save movement file
    rp_file = fullfile(func_dir,'rp_concat.txt');
    save(rp_file,'R','-ascii');
     
    % Onsets file
    onsets = fullfile(base_dir,'onsets',[subject_name '_onsets.mat']);

    % Specify batch
    load('glm_batch.mat');
    matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(glm_dir);    
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = P;
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = cellstr(onsets);
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = cellstr(rp_file);    
    
    if exist(fullfile(glm_dir,'SPM.mat'),'file')
        delete(fullfile(glm_dir,'SPM.mat'));
    end
    
    % Run GLM specification
    spm_jobman('run',matlabbatch);
    
    % Adjust for concatenation
    spm_fmri_concatenate(fullfile(glm_dir,'SPM.mat'), nv);
    
    % Run GLM estimation
    clear matlabbatch;
    matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr(fullfile(glm_dir,'SPM.mat'));
    matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;    
    spm_jobman('run',matlabbatch);    
    
end
%% Add contrasts
load('contrasts_batch.mat');
for i = 1:n
    fprintf('Subject %d\n',i);
    
    [~,subject_name] = fileparts(subject_dirs{i});
    
    spm_mat = fullfile(base_dir,'GLM',subject_name,'SPM.mat');
    matlabbatch{1}.spm.stats.con.spmmat = cellstr(spm_mat);
    
    spm_jobman('run',matlabbatch);
end
%% Second level (collapsed over group)
load('one_sample_ttest_batch.mat');

out_dir = fullfile(base_dir,'second_level');
if ~exist(out_dir,'file')
    mkdir(out_dir);
end

names = {'task','positive','negative','music','notmusic',...
         'interaction_pve','interaction_nve','response'};

for j = 1:length(names)
    
    con_name = sprintf('con_000%d.nii', j+1);    
    
    % Assemble con images ordered by subject
    P = {};
    for k = 1:20        
        P{k} = fullfile(base_dir,'GLM',sprintf('sub-control%02d',k),con_name);
    end
    for k = 1:19 
        P{k+20} = fullfile(base_dir,'GLM',sprintf('sub-mdd%02d',k),con_name);
    end
    
    d = fullfile(out_dir,names{j});
    if ~exist(d,'file'), mkdir(d); end
    
    if exist(fullfile(d,'SPM.mat'),'file')
        delete(fullfile(d,'SPM.mat'));
    end
    
    matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(d);
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cellstr(P)';
    
    spm_jobman('run',matlabbatch);
end
%% Second level (BDI difference)
load('one_sample_ttest_batch.mat');

out_dir = fullfile(base_dir,'second_level_bdi');
if ~exist(out_dir,'file')
    mkdir(out_dir);
end

names = {'task','positive','negative','music','notmusic',...
         'interaction_pve','interaction_nve','response'};

% Load BDI and z-score
behav   = importdata(fullfile(base_dir,'behavioural_data_from_paper','S1File.csv'));
bdi_col = find(strcmp(behav.colheaders, 'BDI_Tot'));
bdi     = behav.data(:,bdi_col);
bdi     = (bdi - mean(bdi)) ./ std(bdi);

for j = 1:length(names)
    
    con_name = sprintf('con_000%d.nii', j+1);    
    
    % Assemble con images ordered by subject
    P = {};
    for k = 1:20        
        P{k} = fullfile(base_dir,'GLM',sprintf('sub-control%02d',k),con_name);
    end
    for k = 1:19 
        P{k+20} = fullfile(base_dir,'GLM',sprintf('sub-mdd%02d',k),con_name);
    end
       
    % Prepare directory
    d = fullfile(out_dir,names{j});
    if ~exist(d,'file'), mkdir(d); end    
    
    if exist(fullfile(d,'SPM.mat'),'file')
        delete(fullfile(d,'SPM.mat'));
    end
    
    matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(d);
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cellstr(P)';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = bdi;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'BDI';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
    
    % Add second contrast for BDI
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'BDI';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';    
        
    spm_jobman('run',matlabbatch);
end
%% Second level (BDA difference)
load('one_sample_ttest_batch.mat');

out_dir = fullfile(base_dir,'second_level_bai');
if ~exist(out_dir,'file')
    mkdir(out_dir);
end

names = {'task','positive','negative','music','notmusic',...
         'interaction_pve','interaction_nve','response'};

% Load BDI and z-score
behav   = importdata(fullfile(base_dir,'behavioural_data_from_paper','S1File.csv'));
bdi_col = find(strcmp(behav.colheaders, 'BAI_Tot'));
bdi     = behav.data(:,bdi_col);
bdi     = (bdi - mean(bdi)) ./ std(bdi);

for j = 1:length(names)
    
    con_name = sprintf('con_000%d.nii', j+1);    
    
    % Assemble con images ordered by subject
    P = {};
    for k = 1:20        
        P{k} = fullfile(base_dir,'GLM',sprintf('sub-control%02d',k),con_name);
    end
    for k = 1:19 
        P{k+20} = fullfile(base_dir,'GLM',sprintf('sub-mdd%02d',k),con_name);
    end
       
    % Prepare directory
    d = fullfile(out_dir,names{j});
    if ~exist(d,'file'), mkdir(d); end    
    
    if exist(fullfile(d,'SPM.mat'),'file')
        delete(fullfile(d,'SPM.mat'));
    end
    
    matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(d);
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cellstr(P)';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = bdi;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'BAI';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
    
    % Add second contrast for BDI
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'BAI';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';    
        
    spm_jobman('run',matlabbatch);
end