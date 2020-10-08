% -------------------------------------------------------
%
%    ppCreateDartelTemplate - creating group template using dartel
% 
%    Ver. 1.0.0
%
%    Created:         Daniela Zoeller      (4.11.2015)
%    Last modified:   Daniela Zoeller      (4.11.2015)
%
%    Medical Image Processing Lab
%    EPFL - UniGe
%
% ------------------------------------------------------
%
% Helper script beeing called from ppMain:
% Loads GM and WM images of preprocessed nii images and creating template
% out of it
%


% clear all
% close all
% clc


%% set up path
codeBasePath = fullfile(filesep,'Users','daniela','Code');
addpath(genpath(fullfile(codeBasePath,'source','22q11_preprocessingfMRI','01_preprocessing','preprocessing_spm12')));
addpath(genpath(fullfile(codeBasePath,'source','general')));
addpath(genpath(fullfile(codeBasePath,'packages','spm12')));

% generate full path of this file
thisLocation=which('preprocess12.m');
% max-compatibility version of fileparts (old releases don't have ~)
if verLessThan('matlab', '7.13.0')
    [jobsParentPath, quux1, quux2, quux3] = fileparts(thisLocation);
else
    [jobsParentPath, quux1, quux2] = fileparts(thisLocation);
end

%% list of all subject folders containing data to preprocess
% dataBasePath = fullfile(filesep,'Users','daniela','Data2','FunImg');
% % [folderNames, folders] = getSubjectFolderNames(fullfile(dataBasePath));
% [folderNames, folders] = getSubjectFolderNames(fullfile(dataBasePath,'checked'));
% 
% [folderNames2, folders2] = getSubjectFolderNames(fullfile(dataBasePath,'cutTemp'));
% 
% folders = [folders; folders2];
% folderNames = [folderNames;folderNames2];


%% get list of rc1* files and rc2* files
rc1FileList = cell(length(folders),1);
rc2FileList = cell(length(folders),1);
for iDir = 1:length(folders)
    tmp = dir(fullfile(folders{iDir},'anat','Segmented','rc1*'));
    rc1FileList{iDir} = fullfile(folders{iDir},'anat','Segmented',tmp.name);
    tmp = dir(fullfile(folders{iDir},'anat','Segmented','rc2*'));
    rc2FileList{iDir} = fullfile(folders{iDir},'anat','Segmented',tmp.name);
end


%% run dartel job
% % List of open inputs
% % Run DARTEL (create Templates): Images - cfg_files
% % Run DARTEL (create Templates): Images - cfg_files

jobs = {fullfile(jobsParentPath,'jobs','createDartelTemplate_job.mat')};
inputs = cell(2, 1);
inputs{1} = rc1FileList; % Run DARTEL (create Templates): Images - cfg_files
inputs{2} = rc2FileList; % Run DARTEL (create Templates): Images - cfg_files
spm('defaults', 'FMRI');
spm_jobman('serial', jobs, '', inputs{:});

%% save executed job and arrange data
for iDir = 1:length(folders)
    dartelFolder = fullfile(folders{iDir},'anat','DARTEL');
    mkdir(dartelFolder);
    if iDir==1; % save job and templates
        copyfile(jobs{:}, dartelFolder); 
        movefile(fullfile(folders{iDir},'anat','Segmented','Template*'),dartelFolder);
    end
    movefile(fullfile(folders{iDir},'anat','Segmented','u_*'),dartelFolder);
end

%% compute Dartel to MNI
% dartel2MNImat = fullfile(folders{1},'anat','DARTEL','Template_6_2mni.mat');
% if exist(dartel2MNImat,'file'); rmdir(dartel2MNImat,'s'); end
% 
% % List of open inputs
% % Normalise to MNI Space: DARTEL Template - cfg_files
% jobs = {fullfile(codeBasePath,'preprocessing_pipeline/functions/jobs/dartelToMNI_job.m')};
% 
% inputs = cell(1);
% inputs{1} = {fullfile(folders{1},'anat','DARTEL','Template_6.nii')}; % Normalise to MNI Space: DARTEL Template - cfg_files
% spm('defaults', 'FMRI');
% spm_jobman('serial', jobs, '', inputs{:});



