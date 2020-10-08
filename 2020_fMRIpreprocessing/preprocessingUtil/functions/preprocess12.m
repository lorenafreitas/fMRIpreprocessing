% -------------------------------------------------------
%
%    preprocess12 - function to execute initial preprocessing steps
%
%    Created:         Daniela Zoeller      (30.10.2015)
%    Based on: preprocess.m (Jonas Richiardi and Giulia Preti)
%    Last modified:   Lorena Freitas       (27.08.2018)
%
% ------------------------------------------------------
%
% Modification preprocess.m:
%   30.10.2015: - run newSegment with Dartel optionf
%               - no deformation task and registration to MNI
%   15.03.2016: - modification of functions and jobs for SPM12
%   24.08.2018: - create Voxel Displacement Map (VDM) based on FieldMap (L.F.)
%   27.08.2018: - realign and unwarp (L.F.)




function preprocess12(functPath, structPath, varargin)
% Preprocessing pipeline for fMRI data (bias field correction, functional
% realignment, functional-structural coregistration, segmentation,
% normalisation, atlasing)
%
% IN:
%   functPath: full path to raw fatlasFileunctional data
%   structPath: full path to raw structural data
%   varargin:
%       1)  the processing chain to use, e.g.
%       'procChain',{'realign','QC','coregister','newSegment','label'}
%       2) other arguments for quality control
%       'QCcoef', 1.5
%       3) high-pass filtering
%       'highpass', 0
%
% REQUIREMENTS
% - SPM8
% - IBASPM toolbox
%
% INSTALLATION
% - SPM8 must be on your matlab path
%
% ACKNOWLEDGEMENTS
% This file borrows liberally from the following sources
% - Rik Henson's example script for SPM 5 (MRC Cognition and Brain Sciences Unit)
% http://en.wikibooks.org/wiki/SPM-Example_batch_script
% - Various functions in the IBASPM toolbox, from Yasser Aleman-Gomez and
% colleages at the Cuba Neurosciences Center
%
% v1.0 Jonas Richiardi / Medical Image Processing Laboratory
% - initial release for SPM5
% v2.0 Manuel W??thrich
% - SPM8 support
% - use spm jobmanager
% - support for newSegment and unifiedSegmentation
% - port IBASPM auto_labelling to SPM8
% - basic QC capability
% v2.0.1 Jonas Richiardi
% - cross-platformish
% - path settings fixed
% v2.0.1b Jonas Richiardi
% - removed ~ ignore code and changed to 4-output forme for filepart for
% backwards compatibility with pre-7.9 versions e.g. as in Greedy SCC cluster
% - corrected path setting code
% v2.0.2 June 2011 Jonas Richiardi
% - supports smoothing
% - supports coregistration direction setting (FtoS or StoF)
% v2.0.3 June 2011 JR
% - supports custom atlases
% v2.1 Oct 2011 JR
% - fixes a critical bug introduced in v2.0 - wrong interpolation
% parameters were used in new_Auto_Labelling.m
% v2.1.1 Dec 2011 JR
% - stricter dir-and-file checking in pp_loadVolumes
% v2.1.2 Dec 2011 JR
% - add version tagging to tasksDone structure to help traceability
% - bugfix and better existence type checking in pp_loadVolume
% - support for Matlab 7.13 (fileparts syntax)

DEBUGMODE=false; % set this to true to help diagnosing problems

%% check and set necessary paths
if exist('spm.m','file')~=2
    error('SPM seems not to be on your matlab path.');
else
    spmLoc=spm('Dir');
    myPath=path();
    % add 'config' SPM dir to the path
    spmConfigLoc=[spmLoc filesep 'config'];
    if exist(spmConfigLoc,'dir')~=7
        error(['SPM installation seems to be missing a config dir at '...
            spmConfigLoc ', quitting.']);
    else
        addpath(spmConfigLoc);
    end
end

if (DEBUGMODE==true)
    path
end

%% initialize and load data
disp('** STARTING PROCESSING');
Tstart = clock;
spm('defaults', 'fMRI');
pp12_loadVolumes;
%if
if isfield(tasksDone, 'realignUnwarp') && tasksDone.realignUnwarp
    alignFolder = unwarpFolder;
    alignmentPrefix = 'u';
else
    alignmentPrefix = 'r';
end


% generate full path of this file
thisLocation=which('preprocess12.m');
% max-compatibility version of fileparts (old releases don't have ~)
if verLessThan('matlab', '7.13.0')
    [jobsParentPath, quux1, quux2, quux3] = fileparts(thisLocation);
else
    [jobsParentPath, quux1, quux2] = fileparts(thisLocation);
end


%% reset of origin
if tasksTodo.reset
    % -------------------- reset of functional images -------------------------------------
    
    
    cd(functPath);
    
    [P,sts] = spm_select('List',functPath,['(^f?^bold).*\.*.(nii|hdr)']);
    if ~sts, return; else P = cellstr(P); end
    spm_progress_bar('Init',numel(P),'Resetting orientations',...
        'Images Complete');
    for i=1:numel(P)
        V    = spm_vol(P{i});
        M    = V.mat;
        vox  = sqrt(sum(M(1:3,1:3).^2));
        if det(M(1:3,1:3))<0, vox(1) = -vox(1); end
        orig = (V.dim(1:3)+1)/2;
        off  = -vox.*orig;
        M    = [vox(1) 0      0      off(1)
            0      vox(2) 0      off(2)
            0      0      vox(3) off(3)
            0      0      0      1];
        spm_get_space(P{i},M);
        spm_progress_bar('Set',i);
    end
    spm_progress_bar('Clear');
end

%---------------- update tasksTodo -------------------------
tasksDone.reset = 1;
save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');


%% realign

if  tasksTodo.realign
    % -------------------- realign -------------------------------------
    if exist(alignFolder,'dir'); rmdir(alignFolder,'s'); end;
    mkdir(alignFolder);
    jobs = {fullfile(jobsParentPath,'jobs','align_job.mat')};
    spm_jobman('initcfg');
    spm_jobman('serial', jobs, '', functFiles);
    
    % ----------------- move files ----------------------
    movefile(fullfile(functPath,'rbold*'),alignFolder);
    movefile(fullfile(functPath,'rp*'),alignFolder);
    movefile(fullfile(functPath,'mean*'),alignFolder);
    
    copyfile(jobs{:}, jobFolder);
    
    %---------------- update tasksTodo -------------------------
    tasksDone.realign = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
    
    
end

%% realign and Unwarp

if  tasksTodo.realignUnwarp
    
    
    % ------------------ create VDM ---------------------------
    
    vdmFile = spm_select('FPList',[functPath '/../'],'^vdm5');
    
    if ~exist(vdmFile)
        jobs = {fullfile(jobsParentPath,'jobs','createVDM_jobServer.mat')};
        
        % --------------- updateInputs and run -------------------
        % Here, I'm dealing with the chance that a subject may have got out of the
        % scanner during the recording, and that a second or even third fieldmap
        % may have been acquired when they entered the scanner again. When this is
        % the case, I want to use the last fieldmap acquired *before* the current
        % functional sequence being analysed.
        
        % Fieldmap Phase
        % ----------------------------------------
        inputs{1} = cellstr(spm_select('FPList',[functPath '/../'],'^gre_field_mapping_phase.*.nii'));   % phase data
        for iPhaseFile = 1 :size(inputs{1},1)
            dPhaseFilesList(iPhaseFile) = dir(char(inputs{1}(iPhaseFile,:)));
        end
        load([structPath(1:end-3) 'dcmHeaders.mat']);
        for i = 1 : length(dPhaseFilesList)
            c = strsplit(dPhaseFilesList(i).name, '.');
            tmp = getfield(h, c{1});
            seriesNumberPhase(i) = getfield(tmp, 'SeriesNumber');
        end
        
        [orderedDateNumsP,idx] = sort(seriesNumberPhase, 'ascend');
        firstBoldFile = dir(spm_select('FPList', functPath, '^bold', '00001.nii'));
        c = strsplit(firstBoldFile.name, '_');
        thisSequence = [c{1} '_' c{2} '_' c{3} ];
        tmp = getfield(h, thisSequence);
        seriesNumberThisSequence = getfield(tmp, 'SeriesNumber');
        
        % Out of the fieldmaps made *BEFORE* the first functional file, I want the
        % last one created.
        
        seriesDiff = seriesNumberThisSequence - seriesNumberPhase;
        if sum(seriesDiff>0) > 0 % if a fieldmap was collected before this sequence, use this one
            orderedDateNumsP(seriesDiff<0) = [];
            idx(seriesDiff<0) = []; % if no fieldmaps were collected before this sequence, use the first one that was collected
        else
            idx(2:end) = [];
        end
        
        inputs{1} = cellstr(deblank(inputs{1}(idx(end),:)));
        
        % Fieldmap magnitude
        % ----------------------------------------
        inputs{2} = cellstr(spm_select('FPList',[functPath '/../'],'^gre_field_mapping_.*00001.nii')); % magnitude data
        for iMagnitudeFile = 1 :size(inputs{2},1)
            dMagnitudeFilesList(iMagnitudeFile) = dir(char(inputs{2}(iMagnitudeFile,:)));
        end
        
        
        for i = 1 : length(dMagnitudeFilesList)
            c = strsplit(dMagnitudeFilesList(i).name, '.');
            c = strsplit(c{1}, '_');
            c1 = [c{1} '_' c{2} '_' c{3}];
            if isfield(h,c1)
                tmp = getfield(h, c1);
            else
                c = [c{1} '_' c{2} '_' c{3} '_' c{4}  ];
                tmp = getfield(h, c);    
            end
            
            seriesNumberMagnitude(i) = getfield(tmp, 'SeriesNumber');
        end
        
        
        [orderedDateNumsM,idxM] = sort(seriesNumberMagnitude, 'ascend');
        
        
        seriesDiff = seriesNumberThisSequence - seriesNumberMagnitude;
        if sum(seriesDiff>0) > 0 % if a fieldmap was collected before this sequence, use this one
            orderedDateNumsM(seriesDiff<0) = [];
            idxM(seriesDiff<0) = []; % if no fieldmaps were collected before this sequence, use the first one that was collected
        else
            idxM(2:end) = [];
        end
        
        inputs{2} = cellstr(deblank(inputs{2}(idxM(end),:)));
        %         % --------------- updateInputs and run -------------------
        %         % fieldmap phase
        %         inputs{1} = cellstr(spm_select('FPList',[functPath '/../'],'^gre_field_mapping_phase.*.nii'));   % phase data
        %         for iPhaseFile = 1 :size(inputs{1},1)
        %             dPhaseFilesList(iPhaseFile) = dir(char(inputs{1}(iPhaseFile,:)));
        %         end
        %         [~,idx] = sort([dPhaseFilesList.datenum], 'ascend');
        %
        %         inputs{1} = cellstr(deblank(inputs{1}(idx(1),:)));
        %
        %
        %        % fieldmap magnitude
        %        inputs{2} = cellstr(spm_select('FPList',[functPath '/../'],'^gre_field_mapping_.*00001.nii')); % magnitude data
        %         inputs{2} = cellstr(inputs{2}(1));
        
        tmp_funcData = spm_select('FPList',functPath,'^bold.*.nii');
        inputs{3} = cellstr([spmLoc '/toolbox/FieldMap/T1.nii']);
        inputs{4} = cellstr(deblank(tmp_funcData(1,:)));
        spm_jobman('serial', jobs, '', inputs{:});
        copyfile(jobs{:}, jobFolder);
        clear tpm_phase_file tmp_funcData
    end
    
    % --------------- realign and unwarp -------------------
    jobs = {fullfile(jobsParentPath,'jobs','alignUnwarp_job.mat')};
    
    if exist(unwarpFolder,'dir'); rmdir(unwarpFolder,'s'); end;
    mkdir(unwarpFolder);
    inputs = [];
    inputs{1} = functFiles;   % functional files
    inputs{2} = cellstr(spm_select('FPList',[functPath '/../'],'^vdm')); % voxel displacement map
    
    spm_jobman('serial', jobs, '', inputs{:});
    
    % ------------------- move files -------------------------
    movefile(fullfile(functPath,'ubold*'),unwarpFolder);
    movefile(fullfile(functPath,'rp*'), unwarpFolder);
    movefile(fullfile(functPath,'mean*'),unwarpFolder);
    %movefile(fullfile(functPath,'wfmag*'),unwarpFolder); %the forward warped field map magnitude image used for co-registration If using a non-EPI field map, the VDM is used to forward warp the magnitude image
    %which is then coregistered to the EPI. The forward warped image is saved with the filename
    %wfmag_NAME-OF-FIRST-INPUT-IMAGE.img
    
    copyfile(jobs{:}, jobFolder);
    
    
    %---------------- update tasksTodo -------------------------
    copyfile(jobs{:}, jobFolder);
    tasksDone.realignUnwarp = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
    
    alignFolder = unwarpFolder;
    alignmentPrefix = 'u';
    
end




%% QC
if tasksTodo.QC
    
    % -------------------- load ----------------------------------------
    
    transData = load(fullfile(alignFolder, ['rp_' functFilenames(1,1:end-3) 'txt']));
    
    functSize = size(spm_read_vols(spm_vol(functFiles{1})));
    intens = zeros(size(functFiles,1),1);
    intensTop100 = zeros(size(functFiles,1),1);
    angle = zeros(size(transData,1),1);
    transl = zeros(size(transData,1),1);
    maxTransl = zeros(size(functFiles,1),1);
    corners = cell(8,1);
    
    for k = 0:7
        mask = dec2bin(k,3) == '1';
        corners{k+1} = (mask .* functSize)';
    end
    
    for i = 1:size(functFiles,1)
        transMatrix = spm_matrix(transData(i,:));
        quat = dcm2qua(transMatrix(1:3,1:3));
        angle(i) = 2*acos(quat(1))/(2*pi)*360;
        transl(i) = norm(transMatrix(1:3,4));
        for k = 1:8
            dist = norm(corners{k} - transMatrix(1:3,:) * [corners{k};1]);
            if dist > maxTransl(i), maxTransl(i) = dist; end
        end
        
        image = spm_read_vols(spm_vol(functFiles{i}));
        intens(i) = mean(image(:));
        array = sort(image(:),'descend');
        intensTop100(i) = mean(array(1:100));
    end
    
    intens = (intens-mean(intens))/std(intens);
    intensTop100 = (intensTop100-mean(intensTop100))/std(intensTop100);
    
    % ---------- compute outliers --------------------------------------
    x = (1:size(intens,1))';
    c = polyfit(x,intens,2);
    f = polyval(c,x);
    
    % ---------  on mean intensity -------------------------------------
    sortIntens = sort(intens);
    lQ = sortIntens(round(0.25*size(intens,1)));
    uQ = sortIntens(round(0.75*size(intens,1)));
    IQR = uQ - lQ;
    min = f + lQ - QCcoef*IQR;
    max = f + uQ + QCcoef*IQR;
    
    indexInt = [find(intens > max); find(intens < min)];
    
    %     figure;
    %     plot(intens,'b'); hold on;
    %     plot(max,'g');
    %     plot(min,'g');
    %     plot(indexInt,intens(indexInt),'ro');hold off;
    
    %    ---------- compute outliers --------------------------------------
    %     x = (1:size(intensTop100,1))';
    %     c = polyfit(x,intensTop100,2);
    %     f = polyval(c,x);
    
    %----------- on mean intensity of top 100 voxels ------------------
    %     sortIntensTop100 = sort(intensTop100);
    %     lQ = sortIntensTop100(round(0.25*size(intensTop100,1)));
    %     uQ = sortIntensTop100(round(0.75*size(intensTop100,1)));
    %     IQR = uQ - lQ;
    %     min = f + lQ - 1.5*IQR;
    %     max = f + uQ + 1.5*IQR;
    %
    %     indexTop100 = [find(intensTop100 > max); find(intensTop100 < min)];
    %
    %
    %
    %
    %
    %     figure;
    %     plot(intensTop100,'b'); hold on;
    %     plot(max,'g');
    %     plot(min,'g');
    %     plot(indexTop100,intensTop100(indexTop100),'ro');hold off;
    %
    % ------------ apply ---------------------------------------------
    mask = zeros(size(intens,1),1);
    mask(indexInt) = 1;
    
    if highpass
        mkdir(QCFolder);
        alignFilenames = spm_select('List',alignFolder,['^' alignmentPrefix '.*\.' volExt '$']);
        alignFiles=cell(size(alignFilenames,1),1);
        QCFiles=cell(size(alignFilenames,1),1);
        for f=1:size(alignFilenames,1)
            alignFiles{f}=fullfile(alignFolder,alignFilenames(f,:));
            QCFiles{f}=fullfile(QCFolder,alignFilenames(f,:));
        end
        
        for i = size(functFiles,1):-1:1
            if mask(i) == 0, copyfile(alignFiles{i},QCFolder);
            else QCFiles(i) = []; end
        end
        
        highpassFilter(char(QCFiles), 2, 0);
        
        for f=1:size(QCFiles,1)
            delete(QCFiles{f});
        end
    end
    
    
    
    % -------------- save artifacts ----------------------------------
    artifacts = struct('rotation', {angle}, 'translation', {transl},'maxTranslation',{maxTransl},...
        'intensity', {intens}, 'intensityTop100', {intensTop100},'mask',mask);
    save(fullfile(alignFolder, 'artifacts'), 'artifacts');
    
    %---------------- update tasksTodo -------------------------
    tasksDone.QC = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% coregister: estimate
if  tasksTodo.coregister
    
    
    jobs = {fullfile(jobsParentPath,'jobs','coreg_job.mat')};
    inputs = cell(3, 1);
    coregDirection;
    %---------------- load meanfile -------------------------
    functFilename = spm_select('List',alignFolder,['^mean.*\.*']);
    meanFile=fullfile(alignFolder,functFilename(1,:));
    %
    if strcmp(coregDirection,'StoF')
        % coreg struct to func (change headers of struct file)
        
        inputs{1} = cellstr(meanFile);
        inputs{2} =  cellstr(structFile); % Coreg: Estimate: Source Image - cfg_files
        inputs{3} = {''};
    else
        % coreg funct to struct (change headers of func files)
        inputs{1} =  cellstr(structFile); % Coreg: Estimate: target Image
        inputs{2} = cellstr(meanFile);
        tmp_others=cellstr(spm_select('List',alignFolder,['^' alignmentPrefix 'bold.*\.' volExt '$']));
        tmp_others=cellfun(@(x) fullfile(alignFolder,x),tmp_others,'UniformOutput',false); % prepend dir
        inputs{3} = tmp_others; % others (all other funcs)
    end
    spm_jobman('serial', jobs, '', inputs{:});
    
    %---------------- update tasksTodo -------------------------
    copyfile(jobs{:}, jobFolder);
    tasksDone.coregister = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% coregister: estimate and and reslice
if  tasksTodo.coregReslice
    
    jobs = {fullfile(jobsParentPath,'jobs','coregReslice_job.mat')};
    inputs = cell(4, 1);
    if strcmp(coregDirection,'FtoS')
        inputs{1} = cellstr(structFile);  % ref (fixed) -> struct
        inputs{2} = {fullfile(alignFolder, ['mean' functFilenames(1,:)])}; % source (moved) -> func
        tmp_others=cellstr(spm_select('List',alignFolder,['^' alignmentPrefix 'bold.*\.' volExt '$']));
        tmp_others=cellfun(@(x) fullfile(alignFolder,x),tmp_others,'UniformOutput',false); % prepend dir
        inputs{3} = tmp_others; % others (all other funcs)
        inputs{4} = 'r'; % prefix
    elseif strcmp(coregDirection,'FtoF')
        inputs{1} = cellstr(refMeanFunc);  % ref (fixed) -> meanfunc
        inputs{2} = {fullfile(alignFolder, ['mean' functFilenames(1,:)])}; % source (moved) -> func
        tmp_others=cellstr(spm_select('List',alignFolder,['^' alignmentPrefix 'bold.*\.' volExt '$']));
        tmp_others=cellfun(@(x) fullfile(alignFolder,x),tmp_others,'UniformOutput',false); % prepend dir
        inputs{3} = tmp_others; % others (all other funcs)
        inputs{4} = 'r'; % prefix
    else
        error(['unknown coregDirection: ' coregDirection]);
    end
    spm_jobman('serial', jobs, '', inputs{:});
    
    %---------------- update tasksTodo -------------------------
    copyfile(jobs{:}, jobFolder);
    tasksDone.coregReslice = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% smooth

if  tasksTodo.smooth
    jobs = {fullfile(jobsParentPath,'jobs','smooth_job.mat')};
    inputs = cell(3, 1);
    % select all funcs
    tmp_data=cellstr(spm_select('List',alignFolder,['^' alignmentPrefix '.*\.' volExt '$']));
    if strcmp(alignmentPrefix, 'u')
        tmp_data{end+1}=['meanu' functFilenames(1,:)];
    else
        tmp_data{end+1}=['mean' functFilenames(1,:)];
    end
    tmp_data=cellfun(@(x) fullfile(alignFolder,x),tmp_data,'UniformOutput',false); % prepend dir
    inputs{1} = tmp_data;   % functional data
    inputs{2} = smoothFWHM; % smoothing kernel specs
    inputs{3} = ['s' num2str(smoothFWHM(1))];
    spm_jobman('serial', jobs, '', inputs{:});
    
    %---------------- update tasksTodo -------------------------
    copyfile(jobs{:}, jobFolder);
    tasksDone.smooth = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end


%% Segmentation and normalization ("New segmentation" in SPM8)
if tasksTodo.segment == SEGMENT
    % --------- do segmentation ---------------------------------
    fprintf('%s','* Segmentation... ');
    mkdir(segFolder);
    
    jobs = {fullfile(jobsParentPath,'jobs','segmentWithWarped_job.mat')};
    spm_jobman('serial', jobs, '', cellstr(structFile));
    
    % arrange files
    movefile(fullfile(structPath,'i*'), segFolder);
    movefile(fullfile(structPath,'y*'), segFolder);
    movefile(fullfile(structPath,'c*'), segFolder);
    movefile(fullfile(structPath,'*seg8.mat'), segFolder);
    
    mkdir(mniFolder);
    movefile(fullfile(structPath,'w*'),mniFolder);
    if size(dir(fullfile(structPath,'mw*')),1)
        movefile(fullfile(structPath,'mw*'),mniFolder);
    end
    
    copyfile(jobs{:}, jobFolder);
    
    fprintf('%s\n','Segmentation done.');
    
    % ------------ update tasksTodo -------------------------------
    tasksDone.segment = SEGMENT;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% Segmentation + DARTEL
if tasksTodo.segment == SEGMENTDARTEL
    % --------- do segmentation ---------------------------------
    
    if ~exist(fullfile(segFolder,'rc1t1_mprage_sag_p2_iso.nii'), 'file')
        mkdir(segFolder);
    else
        warning('rc1* file already exists. Recalculating...');
    end
    
    
    jobs = {fullfile(jobsParentPath,'jobs','segmentDartel_jobLF.mat')};
    spm_jobman('serial', jobs, '', cellstr(structFile));
    
    % arrange files
    movefile(fullfile(structPath,'i*'), segFolder);
    movefile(fullfile(structPath,'y*'), segFolder);
    movefile(fullfile(structPath,'c*'), segFolder);
    movefile(fullfile(structPath,'rc*'), segFolder);
    movefile(fullfile(structPath,'*seg8.mat'), segFolder);
    
    copyfile(jobs{:}, jobFolder);
    
    % ------------ update tasksTodo -------------------------------
    tasksDone.segment = SEGMENTDARTEL;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end





%% IBASPM Labelling
if 	tasksTodo.label
    fprintf('%s','* Atlasing... ');
    list = spm_select('List',segFolder,['^c.*\.' volExt '$']);
    segmentFiles = cell(size(list,1),1);
    for i = 1:size(list,1)
        segmentFiles{i}= fullfile(segFolder, list(i,:));
    end
    
    %if (tasksDone.segment == SEGMENT) || (tasksDone.segment == SEGMENTDARTEL), deform = fullfile(segFolder, ['iy_' structFilename]);
    %elseif tasksDone.segment == OLDSEGMENT, deform = fullfile(normFolder, [ structFilename(1:end-4) '_sn.mat']); end
    deform = fullfile(segFolder, ['iy_' structFilename]);
    
    if atlasType(1)=='G' %Greicius
        %repeat the labeling for each atlas map
        for aaa=1:size(atlasFiles,2)
            new_Auto_Labelling(segmentFiles, atlasFiles{aaa}.fname, deform, atlasFolder, []);
        end
    end
    
    if atlasType(1)=='A' %AAL
        new_Auto_Labelling(segmentFiles, atlasFile, deform, atlasFolder, []);
    end
    
    
    fprintf('%s\n','Atlasing done.');
    
    % I comment it so I don't update the labeling job and I can redo it
    % with different atlases...Giulia
    %---------------- update tasksTodo -----------------------
    tasksDone.label = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% end
Ttotal=etime(clock, Tstart);
disp(['** DONE PROCESSING. Total time: ' num2str(Ttotal/60,'%3.1f') ' min.']);

