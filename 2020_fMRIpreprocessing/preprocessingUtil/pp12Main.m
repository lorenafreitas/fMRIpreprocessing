% -------------------------------------------------------
%
%    pp12Main - loading nifti data and executing first preprocessing steps
%
%    Ver. 1.0.0
%
%    Created:         Daniela Zoeller      (30.10.2015)
%    Based on: preprocROItimeCourseExample.m (Jonas Richiardi and Giulia Preti)
%    Last modified:   Daniela Zoeller      (15.03.2016)
%
% ------------------------------------------------------
%
% should be executed with matlab2015 and spm12
%
% images should already be converted to NIFTI and stored in the
% dataBasePath under "anat" and "func/run_0001"
%

clear all
close all
clc


%% set up path
codeBasePath = fullfile(filesep,'Users','daniela','Code');
addpath(genpath(fullfile(codeBasePath,'source','22q11_preprocessingfMRI','functions')));
addpath(genpath(fullfile(codeBasePath,'source','general')));
addpath(genpath(fullfile(codeBasePath,'packages','spm12')));

%% list of all subject folders
% dataBasePath = fullfile(filesep,'Users','daniela','Data','FunImg');
dataBasePath = fullfile(filesep,'Volumes','WD My Passport','Data','FunImg');
opts.considerFolders = {'cutTempIncluded','checked','cutTemp'};
[~, folders, folderNames] = readScans(opts,dataBasePath);

%% atlas file
AALfile = fullfile(filesep,'Users','daniela','Data','AtlasData','AAL90_correctLR.nii');

%% run preprocessing
for iDir = 1:length(folders)
    fprintf(['-------------------------\n' folderNames{iDir} '\n-------------------------\n']);
    structPath = fullfile(folders{iDir},'anat');
    functPath = fullfile(folders{iDir},'func','run_0001'); % only first run of functional will be utilized, if quality is bad: rearrange data and run again
    preprocess12(functPath,structPath,'procChain',{'reset','QC','realign','smooth','coregister','segmentDartel','label'},...
      'QCcoef',1.5,'smoothFWHM',[6 6 6],'atlasFile',AALfile,'atlasType','AAL');
end


%% create Dartel templates
% ppCreateDartelTemplate % path containing rc* files should be adapted!





