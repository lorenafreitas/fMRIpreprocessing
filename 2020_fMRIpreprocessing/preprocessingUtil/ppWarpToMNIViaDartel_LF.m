function out = ppWarpToMNIViaDartel_LF(b, group)



% tasks to warp
tasks = {'Rest'}; %RealityFiltering2

% Unit test
if nargin <1
    b = initialize_vars('018_mind_p', 'controls', 'pre');
end




%% Template name
dartelDIR = '/Volumes/EPFL_Lorena/BtP/Data/DARTELtemplate/';
templateName='Template_SDBOLD_20190313.nii';
template{1} = [dartelDIR templateName(1:end-4) '_6.nii'];

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
tmp = spm_select('FPListRec',b.structData,['^u*']);

% If deformField has been previously calculated, recalculate it to avoid having transformation errors!
if ~isempty(tmp)
    % First, re-segment structural data to generate a new c1* file
    preprocess12(b.dataDir,b.structData,'procChain',{'segmentDartel'});
end


% If deformField does not exist, create it!
fprintf('Creating deform field to template %s for subject: %s...\n',templateName, b.curSubj);
ppDartelToExistingTemplate_LF(b, templateName(1:end-4));
fprintf('\nDone!\n');
tmp = spm_select('FPListRec',b.structData,['^u*']);


for iDeformField = 1 :size(tmp,1)
    dfFilesList(iDeformField) = dir(deblank(tmp(iDeformField,:)));
end
[~,idx] = sort([dfFilesList.datenum], 'descend');


deformField{1} = deblank(tmp(idx(1),:)); % use most recently calculated deformfield


for it = 1:length(tasks)
    fprintf('Warping task: %s to MNI via DARTEL...\n',tasks{it});
    
    if strcmp(tasks{it}, 'Rest')
        prefix = '^CovVolt.*'; % select prefix of files to be warped
    else
        prefix = '^s6ubold.*';
    end
    
    filesToWarpChar = spm_select('FPListRec',[b.dataDir tasks{it} '/'],prefix);
    
    
    for iftw = 1:size(filesToWarpChar,1)
        filesToWarp{iftw} = deblank(filesToWarpChar(iftw,:));
    end
    
    if exist(b.dataDir, 'dir')
        outDir = fullfile(b.dataDir, tasks{it}, 'warpedMNIviaDartel');
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        else
            warning(['warpedMNIviaDartel folder already exists for task ' tasks{it} ' for subject ' b.curSubj]);
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
end