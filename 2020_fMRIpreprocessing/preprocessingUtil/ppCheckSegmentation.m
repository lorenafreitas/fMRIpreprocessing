

close all;
clear;
clc;


%% set up folders
addpath(genpath(fullfile(filesep,'Users','daniela','Code','source','22q11_preprocessingfMRI','01_preprocessing','preprocessing_spm12',filesep)));
addpath(genpath(fullfile(filesep,'Users','daniela','Code','packages','spm12')));
addpath(genpath(fullfile(filesep,'Users','daniela','Code','source','general')));

%% list of all subject folders containing data to preprocess
dataBasePath = fullfile(filesep,'Users','daniela','Data','FunImg','checked',filesep);
% load([dataBasePath 'folderNamesList20151028.mat']); % list of subjects with double points
% subjBasePath = fullfile(dataBasePath,'excluded');
subjBasePath = fullfile(dataBasePath);%,'cutTemp');

%% list of all subject folders containing data to check
[folderName, folders] = getSubjectFolderNames(subjBasePath);

%% visual check of realigned and coregistered data
checkDir = fullfile(dataBasePath,'segOK');%,'cutTemp');
exclDir = fullfile(dataBasePath);%,'excluded');


%%
for idx = 1:length(folderName)
    fprintf(['-------------------------\n' folderName{idx} '\n-------------------------\n']);
    subjDir = fullfile(subjBasePath,folderName{idx});
    filesOk = 0;
    
    if ~exist(subjDir,'dir')
        continue
    end
    
    % paths:
    structPath = fullfile(subjDir,'anat');
    % files:
    segFiles = dir(fullfile(structPath,'Segmented','c*.nii'));
    for iSegFile = 1:length(segFiles)
        images{iSegFile}=fullfile(structPath,'Segmented',segFiles(iSegFile).name);
    end

    if size(segFiles,1)
        spm_check_registration(char(images));
        pause;

        answer = input('Images ok? (Press N, if not)','s');
        if ~strcmp(answer,'n') && ~strcmp(answer,'N')
            filesOk = 1;
        end
    end
    
    if filesOk && ~strcmp(checkDir,subjBasePath)
        % when all ok, move everything
%         movefile(subjDir,fullfile(checkDir,folderName{idx})); % if ok, move to checked
    elseif ~filesOk && ~strcmp(exclDir,subjBasePath)
        % when no runs are ok, exclude subject
        movefile(subjDir,fullfile(exclDir,folderName{idx})); % if not ok, move to excluded
    end
end
