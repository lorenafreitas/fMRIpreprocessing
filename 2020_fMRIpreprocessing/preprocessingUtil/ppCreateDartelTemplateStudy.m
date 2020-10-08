% -------------------------------------------------------
%
%    ppCreateDartelTemplateStudy - creating group template using dartel
%
%    Created:         Daniela Zoeller     (04.11.2015)
%    Last modified:   Lorena Freitas      (30.07.2018)
%
% ------------------------------------------------------
%
% script to create a study-specific Template using Dartel

% 20190123: all preterms
% 20190313: preterms + controls


clear all
close all
clc


%% Set up Path 

% Generate full path of this file
thisLocation=which('preprocess12.m');

% Max-compatibility version of fileparts (old releases don't have ~)
if verLessThan('matlab', '7.13.0')
    [jobsParentPath, quux1, quux2, quux3] = fileparts(thisLocation);
else
    [jobsParentPath, quux1, quux2] = fileparts(thisLocation);
end

%% Selects files created during the segmentation step, in order to generate DARTEL template
dataBasePath = fullfile(filesep,'media','miplab-nas','Data2','Nikolina_VAV_RF');
date = datestr(date,'YYYYmmDD');


rc1FileListChar = spm_select('FPListRec',dataBasePath,'^rc1.*'); 
rc1FileList = cell(size(rc1FileListChar,1),1);for i = 1:size(rc1FileListChar,1), rc1FileList{i} = deblank(rc1FileListChar(i,:)); end
rc2FileListChar = spm_select('FPListRec',dataBasePath,'^rc2.*');
rc2FileList = cell(size(rc2FileListChar,1),1);for i = 1:size(rc2FileListChar,1), rc2FileList{i} = deblank(rc2FileListChar(i,:)); end

%% Run DARTEL job
% % List of open inputs
jobs = {fullfile(jobsParentPath,'jobs','createDartelTemplateSpecifyName_job.mat')};
inputs = cell(2, 1);
inputs{1} = rc1FileList; 
inputs{2} = rc2FileList; 
inputs{3} = ['Template_VAV_BOLD_' date]; 
spm('defaults', 'FMRI');
spm_jobman('serial', jobs, '', inputs{:});

%% Save executed job and arrange data
uFiles = spm_select('FPListRec',dataBasePath,['^u_.*' date '*']);
for iDir = 1:size(uFiles,1)
    
    [filepath,name,ext] = fileparts(uFiles(iDir,:));
    dartelFolder = [filepath '/../DARTEL/'];
    if ~exist(dartelFolder, 'dir'), mkdir(dartelFolder); end
    if exist(deblank(uFiles(iDir,:)), 'file'), movefile(deblank(uFiles(iDir,:)),dartelFolder); end
    
    
    if iDir==1 && exist(fullfile(filepath,'Template*'), 'file') % save job and templates
        copyfile(jobs{:}, dartelFolder); 
        movefile(fullfile(filepath,'Template*'),dartelFolder);
    end
  
end




