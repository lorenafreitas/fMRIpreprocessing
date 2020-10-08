% 10.5.2016: modified to use SPM deformations instead of "create warped"
% 21.9.2016: volterra 24 covariates, study specific template

function out = ppWarpC1ToMNIViaDartel_LF(b, group)


% Unit test
if nargin <1
    b = initialize_vars('01_vav_p', 'preterm');
end




%% Template name

templateName=['Template_VAV_BOLD_' date '.nii'];
template{1} = ['/Volumes/EPFL_Lorena/BtP/Data/DARTELtemplate/' templateName(1:end-4) '_6.nii';]
% generate full path of this file
thisLocation=which('preprocess12.m');
% max-compatibility version of fileparts (old releases don't have ~)
if verLessThan('matlab', '7.13.0')
    [jobsParentPath, ~,~ , ~] = fileparts(thisLocation);
else
    [jobsParentPath, ~, ~] = fileparts(thisLocation);
end
jobs = {fullfile(jobsParentPath,'jobs','deformToMNIWithDartelLF_job.mat')};



% Find deform field
dataBasePath = b.dataDir;
deformField{1} = spm_select('FPListRec',b.structData,['^u.*' date '*']);


%fprintf('Warpink task: %s to MNI via DARTEL...\n',tasks{it});
filesToWarpChar = spm_select('FPListRec',[b.structData '/'],'^rc1.*');

for iftw = 1:size(filesToWarpChar,1)
    filesToWarp{iftw} = deblank(filesToWarpChar(iftw,:));
end

if exist(b.dataDir, 'dir')
    outDir = fullfile(b.structData, 'warpedMNIviaDartel');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    else
        warning(['warpedMNIviaDartel folder already exists for subject ' b.curSubj]);
    end
    % warp files to Dartel space
    % --------------------------
    if ~strcmp(deformField{1}(end-3:end),'.nii')
        warning('subject not preprocessed yet, please finalize processing! Skipping...\n');
        out = 0;
    end
    
    tmp=strfind(filesToWarp{1},'/');tmp=tmp(end);
    if exist(fullfile(outDir,['w' filesToWarp{1}(tmp+1:end)]),'file');
        warning('subject already warped, skipping...\n');
        out = 0;
    end
    
    inputs{1} = template; % Deformations: Flow fields - cfg_files
    inputs{2} = deformField; % Deformations: Images - cfg_files
    inputs{3} = filesToWarp'; % Deformations: Output directory - cfg_files
    spm('defaults', 'FMRI');
    spm_jobman('run', jobs, inputs{:});
    out = 1;
else
    warning('Data folder %s does not exist!', b.dataDir);
end


end