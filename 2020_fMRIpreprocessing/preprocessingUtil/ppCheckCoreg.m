

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
subjBasePath = fullfile(dataBasePath);

%% list of all subject folders containing data to check
[folderNames, folders] = getSubjectFolderNames(subjBasePath);

%% visual check of realigned and coregistered data
exclDir = fullfile(dataBasePath,'excluded');
movDir = fullfile(dataBasePath,'excluded','moved');
cutDir = fullfile(dataBasePath,'excluded','cut');
checkDir = dataBasePath;
% exclDir = fullfile(dataBasePath,'excluded');

for idx = 1:length(folderNames)
%     if ~strcmp(folderNames{idx},'5287_3_3T')
%         continue
%     end
    
    fprintf(['-------------------------\n' folderNames{idx} '\n-------------------------\n']);
    subjDir = fullfile(subjBasePath,folderNames{idx});
    filesOk = 0;
    
    if ~exist(subjDir,'dir')
        continue
    end
    
    iRun=1;
    % paths:
    rePath = fullfile(subjDir,'func',sprintf('run_%0.4d', iRun),'realigned');
    structPath = fullfile(subjDir,'anat');
    % files:
    meanFile = dir(fullfile(rePath,'mean*'));
    tmp = fullfile(structPath,{'s*.nii';'t1*.nii';'ms*.nii'});
    structFile = []; iFile = 1;
    while ~size(structFile,1) && iFile <= length(tmp)
        structFile = dir(tmp{iFile});
        iFile = iFile + 1;
    end

%         normFile = dir(fullfile(structPath,'Normalized','w*.nii'));
    reParams = dir(fullfile(rePath,'rp*.txt'));

    if size(meanFile,1) && size(structFile,1) && size(reParams,1) % check file existence
        if ppCheckMovement(fullfile(rePath,reParams.name)) % check movement limit values
            spm_check_registration(char(fullfile(structPath,structFile.name),...
                fullfile(rePath,meanFile.name)),...
                {folderNames{idx},folderNames{idx}}); % visual check
            pause;

            answer = input('Images ok? (Press N, if not)','s');
            if ~strcmp(answer,'n') && ~strcmp(answer,'N')
                filesOk(iRun) = 1;
            end
        else
            filesOk(iRun)=-1; %subject moved
        end
    
    
        if all(filesOk==1) && ~strcmp(checkDir,subjBasePath)
            % when all ok, move everything
            movefile(fullfile(subjDir,fullfile(checkDir,folderNames{idx}))); % if ok, move to checked
        elseif any(filesOk==1) && ~strcmp(checkDir,subjBasePath)
            if exist(fullfile(checkDir,folderNames{idx}),'dir') || mkdir(fullfile(checkDir,folderNames{idx}))
                % when only some runs are ok, move the ones that are ok
                movefile(fullfile(subjDir,'anat'),fullfile(checkDir,folderNames{idx},'anat')); % move anat to checked
                movefile(fullfile(subjDir,'jobs'),fullfile(checkDir,folderNames{idx},'jobs'));
                for iRun = find(filesOk)
                    runPath = fullfile(subjDir,'func',sprintf('run_%0.4d', iRun));
                    if exist(fullfile(checkDir,folderNames{idx},'func'),'dir') ...
                            || mkdir(fullfile(checkDir,folderNames{idx},'func'))
                        movefile(runPath,fullfile(checkDir,folderNames{idx},...
                            'func',sprintf('run_%0.4d', iRun))); % move good run
                    end
                end
            end
            rmdir(subjDir,'s');
        elseif ~any(filesOk==1) && ~strcmp(cutDir,subjBasePath)
            if any(filesOk==-1) %subject moved
                movefile(subjDir,fullfile(movDir,folderNames{idx})); % if not ok, move to excluded
            else
                % when no runs are ok, exclude subject
                movefile(subjDir,fullfile(cutDir,folderNames{idx})); % if not ok, move to excluded
            end
        end
    end
end
