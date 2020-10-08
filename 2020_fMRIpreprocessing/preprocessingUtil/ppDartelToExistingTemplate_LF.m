% -------------------------------------------------------
%
%    ppDartelToExistingTemplate_LF - create flowfields for subjs not in
%    template
%
%    Ver. 1.0.0
%
%    Created:         Daniela Zoeller      (4.11.2015)
%    Last modified:   Lorena Freitas       (7.08.2018)
%
%    Medical Image Processing Lab
%    EPFL - UniGe
%
% ------------------------------------------------------
%
% Helper script beeing called from ppWarpToMNIViaDartel_LF:
% Creates flowfield files (u_*) for subjects not included in template
%


%clear all
%close all
%clc


function out = ppDartelToExistingTemplate_LF(b, templateName)


if nargin < 1
    %% Template name for which to execute
    templateName='Template_SDBOLD_20190313';
    b = initialize_vars('504_mind_c', 'controls', 'pre');
end


% generate full path of this file
thisLocation=which('preprocess12.m');
% max-compatibility version of fileparts (old releases don't have ~)
if verLessThan('matlab', '7.13.0')
    [jobsParentPath, quux1, quux2, quux3] = fileparts(thisLocation);
else
    [jobsParentPath, quux1, quux2] = fileparts(thisLocation);
end



%% list of all subject folders
dataBasePath = b.dataDir;
templateDir  = '/Volumes/EPFL_Lorena/BtP/Data/DARTELtemplate/';


%% get list of rc1* files and rc2* files (only subjects that have no u_* file yet)

%deformField = spm_select('FPListRec',b.structData,['^u.*\' templateName]);

%if exist(deformField,'file')
%    fprintf([deformField ': flow field already exists, skipping\n']);
%end
rc1File{1} = spm_select('FPListRec',[b.structData '/'],'^rc1.*');
rc2File{1} = spm_select('FPListRec',[b.structData '/'],'^rc2.*');


%% get list of existing template files

templatesChar = spm_select('FPList',templateDir,[templateName '.*\.nii']);

% Use only templates 1-6, since Template_0.nii is generated only after
% rigid alignment.
for iTemp = 2:7
    templates{iTemp-1,1} = fullfile(deblank(templatesChar(iTemp,:)));
end


%% run dartel job
% % List of open inputs
% % Run DARTEL (create Templates): Images - cfg_files
% % Run DARTEL (create Templates): Images - cfg_files

jobs = {fullfile(jobsParentPath,'jobs','dartelToExistingTemplate_job.mat')};
inputs = cell(8, 1);
inputs{1} = rc1File; % Run DARTEL (create Templates): Images - cfg_files
inputs{2} = rc2File; % Run DARTEL (create Templates): Images - cfg_files
%inputs(3:8)= templates;
for iTemplates = 1:length(templates)
    inputs{iTemplates+2} = cellstr(templates{iTemplates});
end
spm('defaults', 'FMRI');
spm_jobman('serial', jobs, '', inputs{:});

% %% save executed job and arrange data
% for iDir = 1:length(foldersToProcess)
%     dartelFolder = fullfile(foldersToProcess{iDir},'anat','DARTEL');
%     mkdir(dartelFolder);
%     fNames=dir(fullfile(foldersToProcess{iDir},'anat','Segmented','u_rc1*'));
%     movefile(fullfile(foldersToProcess{iDir},'anat','Segmented',fNames.name),fullfile(dartelFolder,[fNames.name(1:end-4) '_' templateName '.nii']));
% end
end




