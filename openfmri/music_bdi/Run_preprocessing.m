% Path to the downloaded and unzipped data
data_dir = 'C:\original_scan_data\music_mdd\ds171_R1.0.0';

% Find subjects
subject_dirs = spm_select('FPList',data_dir,'dir','^sub-.*');
subject_dirs = cellstr(subject_dirs);

% Output directory
base_dir = pwd;

n = length(subject_dirs);

TR = 3; % Seconds

%% Untar images
for i = 1:n
    fprintf('Untarring subject %d of %d\n',i,n);
    
    anat_dir = fullfile(subject_dirs{i},'anat');
    func_dir = fullfile(subject_dirs{i},'func');
    
    struct = cellstr(spm_select('FPList',anat_dir,'.*_T1w.nii.gz$'));
    funcA  = cellstr(spm_select('FPList',func_dir,'.*task-nonmusic.*gz$'));
    funcB  = cellstr(spm_select('FPList',func_dir,'.*task-music.*gz$'));
    
    gunzip(struct{1}, anat_dir);
    
    for j = 1:length(funcA)
        gunzip(funcA{j},func_dir);
    end
    for j = 1:length(funcB)
        gunzip(funcB{j},func_dir);
    end
end
%% Preprocess images (realign, segment, coregister, normalise, smooth)
load('preprocess_batch.mat');
for i = 1:n
    fprintf('Preprocessing subject %d of %d\n',i,n);
    
    % Subject's directories
    anat_dir = fullfile(subject_dirs{i},'anat');
    func_dir = fullfile(subject_dirs{i},'func');
    
    % Functionals
    P1 = cellstr(spm_select('ExtFPList',func_dir,'^s.*_run-1_bold.nii$'));
    P2 = cellstr(spm_select('ExtFPList',func_dir,'^s.*_run-2_bold.nii$'));
    P3 = cellstr(spm_select('ExtFPList',func_dir,'^s.*_run-3_bold.nii$'));
    P4 = cellstr(spm_select('ExtFPList',func_dir,'^s.*_run-4_bold.nii$'));
    P5 = cellstr(spm_select('ExtFPList',func_dir,'^s.*_run-5_bold.nii$'));
    
    matlabbatch{1}.spm.spatial.realign.estwrite.data = {P1 P2 P3 P4 P5}';    
    
    % Anatomical
    anat = spm_select('FPList',anat_dir,'^s.*_T1w.nii$');
    matlabbatch{2}.spm.spatial.preproc.channel.vols = cellstr(anat);
    
    save('delme.mat','matlabbatch');
    spm_jobman('run',matlabbatch);
end
%% Make average structural
load('create_mean_structural.mat');
P = spm_select('FPListRec',data_dir,'wbrain.nii');
matlabbatch{1}.spm.util.imcalc.input = cellstr(P);
matlabbatch{1}.spm.util.imcalc.outdir = cellstr(base_dir);
spm_jobman('run',matlabbatch);