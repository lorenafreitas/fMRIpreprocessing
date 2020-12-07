% JOB_GLM: This function performs a GLM Analysis on a subject from
% the Building the Path fMRI dataset for the Reality Filtering task.
%____________________________________________________________________________
% Copyright (C) 2016 Lorena Freitas; MIP:Lab. 


% $Id: job_GLM_RF.m 11 2016-06-08 16:26:24F Lorena $

function [matlabbatch] = job_GLM_RF(b)


%% ------------------------------------------------------------------------
% PARAMETER SETTING
%--------------------------------------------------------------------------

% MOTION
RPARAM24  = 1; % 1 => Reduced Voltera expansion of realignment parameters using SVD
DESPIKING = 1; % 1 => Includes scan nulling regressors
RAD = 50; % Head size to convert motion to mm

% REGRESSORS
TRIANGLE_REGRESSORS=0; % Model brain response as triangle? 1 = yes; 0 = no;
BLOCK_DESIGN = 0; % 1 = block design; 0 = event related, unless TRIANGLE_REGRESSORS == 1;


% REPETITION TIME
TR = 0.72;
date = '01-Apr-2017'; % date of the saved analysis

% Run for both tasks: RF1 and RF2
task = 'RF1';
taskFolder = 'RealityFiltering1';


fprintf(['\n\n========================================================================\n',...
    'Running ' task ' for ' b.curSubj '!\n', ...
    '========================================================================\n']);


% ------------------------------------------------------------------------
% TASK ONSETS
% Onsets are lists of onset times *IN SECONDS* 
% For an event related design, durations = 0. 
% ------------------------------------------------------------------------
[onsetD1, onsetT1, durations, b] = loadBehaviouralData(b, task);

% ------------------------------------------------------------------------
% SELECT FUNCTIONAL IMAGES
%--------------------------------------------------------------------------
prefix = 'ws6ubold';
files  = spm_select('FPListRec',fullfile(deblank(b.dataDir), taskFolder, 'unwarped'),['^' prefix '.*\.nii$']);

% -------------------------------------------------------------------------
% SELECT OUTPUT DIRECTORY
%--------------------------------------------------------------------------

analyze_dir = [b.dataDir filesep taskFolder filesep 'GLM_' task '_' date filesep];
b.currentTaskPath = [b.dataDir filesep taskFolder filesep 'unwarped'];

% -------------------------------------------------------------------------
% LOAD THE TIMING CORRESPONDING TO THE TASKS, FOR REGRESSORS
%--------------------------------------------------------------------------

cond={'D1', 'T1'};
onsets = {onsetD1, onsetT1};
dur = {durations{1}, durations{2}}; % durations of each type of stimulus

% -------------------------------------------------------------------------
% Define Contrasts
%--------------------------------------------------------------------------

contrast_name={'Distractors', 'Targets', 'D1-T1', 'T1-D1'};
contrast_vect={[1 0 ],[0 1 ],[1 -1],[-1 1]};


%%  STATISTICAL PROCESSING
% -------------------------------------------------------------------------
step_count = 1;
mkdir(analyze_dir);

fs=10; % Sampling of the regressor
hrf = spm_hrf(1/fs);

% SPECIFICATIONS
matlabbatch{step_count}.spm.stats.fmri_spec.dir{1,1} = analyze_dir;
matlabbatch{step_count}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{step_count}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{step_count}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{step_count}.spm.stats.fmri_spec.volt = 1;
matlabbatch{step_count}.spm.stats.fmri_spec.global = 'None';
matlabbatch{step_count}.spm.stats.fmri_spec.mask = {''};
matlabbatch{step_count}.spm.stats.fmri_spec.cvi = 'AR(1)';

% -------------------------------------------------------------------------
%  DEALING WITH MOTION
%  -------------------------------------------------------------------------
% Load motion parameter file
rp_file=spm_select('FPList', deblank([b.currentTaskPath]), '^rp_b');
rp_param=load(rp_file);
clearvars rp_file;


% FrameWise Displacement
% -------------------------------------------------------------------------
if RAD % USE ROTATIONS
    % Convert rotations (radians) into translations (mm)
    rp_param_rad = [rp_param(:,1:3) rp_param(:,4:6)*(2*RAD*pi/360)];
    
    % Compute the Framewise Displacement
    FD = [0;sum(abs(diff(rp_param_rad)),2)];
    
else % DO NOT USE ROTATIONS
    FD = [0;sqrt(sum(diff(rp_param(:,1:3)).^2,2))];
end


% If adding "scan nulling regressors":
if DESPIKING
    FDspikes           = FD > 0.5;
    [f2r, percentF2R]  = frames2remove((FDspikes)');
    b.pctFramesRemoved = percentF2R;
    nSpikeRegressors   = sum(f2r);
    spikeRegressors    = zeros(length(FD), nSpikeRegressors);
    spikeCount         = 0;
    for isr = 1:length(FD)
        if f2r(isr) == 1
            spikeRegressors(isr,spikeCount+1) = 1;
            spikeCount                      = spikeCount + 1;
        end
    end
end

b.subjectInfo.nSpikeRegressors = spikeCount;

% If motion parameters file exists, load it, if not, generate it.
% -------------------------------------------------------------------------
if ~exist(char(strcat(b.currentTaskPath, '/rp_sess_new.txt')),'file')
    
    if RPARAM24
        rp_param = [zeros(1,6) ; rp_param];
        rp_param = [rp_param(2:end,:) rp_param(2:end,:).^2 rp_param(1:end-1,:) rp_param(1:end-1,:).^2];
        
        % Save session motion parameters file
        rp_sess = rp_param;%[rp_param(sess(s,1):sess(s,2),:)];
        rp_sess = (rp_sess-repmat(mean(rp_sess),size(rp_sess,1),1))./repmat(std(rp_sess),size(rp_sess,1),1);
        
        % Reduce 24 parameters using SVD
        [U,S,V]=svd(rp_sess,0);
        S=diag(S); % S is a diagonal matrix containing the sqroots of eigenvalues from U or V in descending order
        Svar=100*S.^2/sum(S.^2); % Percentage variance of the single values
        tmp_asvar=find(cumsum(Svar)>99); % Find the components that explain 99% of the motion-related variance
        NPC=min(6,tmp_asvar(1)); % Select either the first 6 SVD or the first components that explain the first 99% (whatever is smaller)
        rp_sess=U(:,1:NPC)*diag(S(1:NPC));
    else
        rp_sess = rp_param;
    end
    
    % If despiking, include despiking regressors
    if DESPIKING
        rp_sess = [rp_sess, spikeRegressors];
    end
    rp_sess_file = fullfile(deblank([b.currentTaskPath]),['rp_sess_new.txt']);
    save(rp_sess_file,'-ascii','-tabs','rp_sess');
    clearvars rp_sess;
else
    rp_sess_file = fullfile(deblank([b.currentTaskPath]),['rp_sess_new.txt']);
end



matlabbatch{step_count}.spm.stats.fmri_spec.sess.scans = cellstr(deblank(files));

if TRIANGLE_REGRESSORS
    % BUILD REGRESSORS
    fs=10; % Sampling of the regressor
    NTR=size(files,1);
    reg = zeros(ceil(NTR*TR*fs),numel(cond));
    for c = 1 : numel(cond)
        tmp_onset = onsets{c}(:);
        for o = 1 : numel(tmp_onset)
            trialPos = round(fs*tmp_onset(o)+1:fs*(tmp_onset(o)+dur{c}(o)));
            triangle = fliplr(0:1/length(trialPos):1);
            reg(trialPos,c)=triangle(1:end-1);
        end
    end
    reg = conv2(reg,hrf);
    
    regressor = zeros(NTR,numel(cond));
    regressor = reg(round(fs*TR/2:TR*fs:TR*fs*NTR),:);
    dregressor = [zeros(1,numel(cond));abs(diff(regressor))];
    for i=1:size(dregressor,2)
        [~,ind_dreg]=findpeaks(dregressor(:,i),'minpeakheight',0.9*max(dregressor(:,i)));
        tmp=zeros(size(dregressor,1),1);
        tmp(ind_dreg)=1;
        dregressor(:,i)=tmp;
    end
    clearvars reg c tmp_onset o n i tmp ind_dreg;
end

for c = 1 : numel(cond)
    
    if TRIANGLE_REGRESSORS
        matlabbatch{step_count}.spm.stats.fmri_spec.sess.regress(c).name = char(cond{c});
        matlabbatch{step_count}.spm.stats.fmri_spec.sess.regress(c).val = regressor(:,c);
    else
        matlabbatch{step_count}.spm.stats.fmri_spec.sess.cond(c).name = char(cond{c});
        matlabbatch{step_count}.spm.stats.fmri_spec.sess.cond(c).onset = onsets{c}(1,:);
        if BLOCK_DESIGN %if design is "block type". If not, it's considered event related.
            matlabbatch{step_count}.spm.stats.fmri_spec.sess.cond(c).duration = dur{c}(1,:);
        else %event related
            matlabbatch{step_count}.spm.stats.fmri_spec.sess.cond(c).duration = zeros(1,length(dur{c}(1,:)));
        end
    end
    
end

matlabbatch{step_count}.spm.stats.fmri_spec.sess.multi_reg{1,1} = rp_sess_file;
matlabbatch{step_count}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{step_count}.spm.stats.fmri_spec.mthresh = 0.8;
   

%% MODEL ESTIMATION
% -------------------------------------------------------------------------

step_count = step_count + 1;
matlabbatch{step_count}.spm.stats.fmri_est.spmmat{1,1} = fullfile(deblank(analyze_dir),'SPM.mat');
matlabbatch{step_count}.spm.stats.fmri_est.method.Classical = 1;

% Contrasts
step_count = step_count + 1;
matlabbatch{step_count}.spm.stats.con.spmmat{1,1} = fullfile(deblank(analyze_dir),'SPM.mat');

%Effects of interest (F-Contrast)
matlabbatch{step_count}.spm.stats.con.consess{1}.fcon.name = 'Effects of interest';
matlabbatch{step_count}.spm.stats.con.consess{1}.fcon.weights = [eye(length(cond)) zeros(length(cond),1)];
matlabbatch{step_count}.spm.stats.con.consess{1}.fcon.name = 'All';
matlabbatch{step_count}.spm.stats.con.consess{1}.fcon.weights = ones(1,length(cond));

% T-Contrasts
for c = 1 : numel(contrast_vect)
    matlabbatch{step_count}.spm.stats.con.consess{c+1}.tcon.name = contrast_name{c};
    matlabbatch{step_count}.spm.stats.con.consess{c+1}.tcon.convec = contrast_vect{c};
    matlabbatch{step_count}.spm.stats.con.consess{c+1}.tcon.sessrep = 'replsc'; % Replicate contrast over sessions and scale to take into account the variable number of session
end
matlabbatch{step_count}.spm.stats.con.delete = 1; % Delete previous contrasts
clearvars c;



%% RUN ANALYSIS
% -------------------------------------------------------------------------
save(fullfile(deblank(analyze_dir),['SPM_analysis_' date '_' task '.mat']),'matlabbatch');
spm_jobman('run',matlabbatch);


end