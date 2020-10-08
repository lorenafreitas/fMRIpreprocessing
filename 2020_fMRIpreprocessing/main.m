% This batch script analyses the Building the Path to Resilience dataset
%__________________________________________________________________________
% Copyright (C) 2016 MIP:Lab, Lorena Freitas

% Lorena Freitas
% $Id: main.m 11 2020-22-09 16:26:24F Lorena $

clear all;

% Get and set paths
%--------------------------------------------------------------------------
scriptdir = pwd;
addpath(genpath('~/spm12/')); % change to path where spm12 is installed
% also add path where this code is :) 

% Specify variables
%--------------------------------------------------------------------------

controls  = {'501_vav_c', '502_vav_c', '503_vav_c',...
    % INCLUDE HERE ALL SUBJECTS :) 
    '033_mind_p','303_mind_p','Controls'};


preterms  = {'01_vav_p', '02_vav_p', '03_vav_p',...
    'Preterm'};
    

groups      = {preterms, controls}; 
interv      = { 'pre'};
batch_functions = {'preproc_pipeline'}; % options: 'ppWarpC1ToMNIViaDartel_LF', 'ppWarpToMNIViaDartel_LF', 'job_getFD';  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize error log
%--------------------------------------------------------------------------
errorlog = {}; ctr=1;

subjectsTime = zeros(1, (sum(cellfun('length',groups))-length(groups))); % #subjs in all groups

% loop over groups
for g=1:length(groups)
    fprintf('\nWorking on group %s...\n',groups{g}{end});
    
    % loop over subjects
    for i=1:length(groups{g})-1 % -1: last item is the group name
        fprintf(['\n\n========================================================================\n',...
            'Working on subject %s\n', ...
            '%d of %d\n',...
            '========================================================================\n'],groups{g}{i},i,(length(groups{g})-1) ); tic;
        
        % loop over batch functions
        for m=1:length(batch_functions)
            fprintf('\nFunction: %s\n',batch_functions{m});
            
            
            for i_int = 1:length(interv)
                % get subject-specific variables
                fprintf('Running %s -intervention analysis \n', interv{i_int});
                b = initialize_vars(groups{g}{i}, groups{g}{end}, interv{i_int});
                
                if exist(b.dataDir, 'dir')
                    
                    if (strcmp(b.curSubj(end), 'c') && i_int==2), continue; end
                    %run matlabbatch job
                    %--------------------------------------------------------------
                    try
                        %run current job, passing along subject-specific inputs
                        
                        batch_output = eval(strcat(batch_functions{m},'(b)'));
                        
                    catch err % if there's an error, take notes & move on
                        errorlog{ctr,1} = groups{g}{i};
                        errorlog{ctr,2} = b.interv;
                        errorlog{ctr,3} = batch_functions{m};
                        errorlog{ctr,4} = err;
                        ctr = ctr + 1;
                        cd(scriptdir);
                        continue;
                    end
                    subjectsTime(i) = toc;
                    cd(scriptdir);
                else
                    warning('Attention: no %s -intervention data found for subject %s', interv{i_int}, b.curSubj);
                end
            end
        end
        
    end
end

% Print how long it took to run the analysis
% ------------------------------------------
fprintf('The analysis of %d subjects took %d minutes and %f seconds\n', ...
    length(subjectsTime),floor(sum(subjectsTime)/60), ...
    rem(sum(subjectsTime),60));

% If any errors occur, log them at the end
% ------------------------------------------
if ~isempty(errorlog)
    disp(errorlog) 
else
    disp('No errors detected.');
end
