% -------------------------------------------------------
%
%    importDICOM - script to convert dicom files using the function
%    dicom2nifti and spm
%
%    Ver. 1.0.0
%
%    Created:         Daniela Zoller      (16.10.2015)
%    Last modified:   Daniela Zoller      (16.10.2015)
%
%    Medical Image Processing Lab
%    EPFL - UniGe
%
% ------------------------------------------------------
%

clc
close all
clear

%% set up folders
% rmpath(genpath(fullfile(filesep,'Users','daniela','Documents','Code','Third Party','spm12')));
addpath(genpath(fullfile(filesep,'Users','daniela','Code','packages','spm12')));
addpath(genpath(fullfile(filesep,'Users','daniela','Code','packages','spmtools')));

codeBasePath = fullfile(filesep,'Users','daniela','Code');
addpath(genpath(fullfile(codeBasePath,'source','general')));

%% get directory with all dicom files
% dicomDir = uigetdir(pwd,'DICOM file directory');
% dicomDir = '/Volumes/Neuroimaging_Projetcs/Daniela_resting/FunRaw';
dicomDir = '/Users/daniela/Data/FunRaw';
[folderName, folderList] = getSubjectFolderNames(dicomDir);

%% determine directory to save the converted files in
% niftiDir = uigetdir(pwd,'converted NIFTI file directory');
niftiDir = '/Users/daniela/Data/FunImg';

%% subjects to skip
subjToSkip = {};%'5157_4_3T','6542_1_3T','6608_1_3T','6623_1_3T','6624_1_3T',...
%     '6650_1_3T','6679_1_3T','6978_1_3T'};


%% convert files of every subject in the folderNameList
for idx = 1:length(folderName)
    fprintf(['-------------------------\n' folderName{idx} '\n-------------------------\n']);
    % skip defined subjects and already checked/excluded subjects
    if ismember(folderName{idx},subjToSkip) ...
            || exist(fullfile(niftiDir,'checked',folderName{idx}),'dir') ...
            || exist(fullfile(niftiDir,'cutTemp',folderName{idx}),'dir') ...
            || exist(fullfile(niftiDir,'excluded',folderName{idx}),'dir')
        continue
    end
    
    dicomSubjDir = fullfile(dicomDir,folderName{idx});
    niftiSubjDir = fullfile(niftiDir,folderName{idx});
    
    % check if there are already converted files
    if exist(fullfile(niftiSubjDir,'anat'),'dir')
        anat_exist = 0;
        file = dir( fullfile(niftiSubjDir,'anat','*.nii'));
        if size(file,1)
            anat_exist = 1;
            fun_exist = 0;
            for iRun = 1:3
                if exist(fullfile(niftiSubjDir,'func',sprintf('run_%0.4d', iRun)),'dir')
                    file = dir(fullfile(niftiSubjDir,'func',sprintf('run_%0.4d', iRun),'*.nii'));
                    if size(file,1)
                        fun_exist = 1;
                    else
                        fun_exist = 0;
                    end
                end
            end
            if fun_exist
                warning('%s : files already converted',folderName{idx});
                continue;
            end
        end
    end
    dicom2nifti('dicom_dir', dicomSubjDir,...
        'subject_dir', niftiSubjDir);
end


