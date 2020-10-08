
% This script performs the preprocessing pipeline for one subject from the
% Building the Path to Resilience project

function OK = preproc_pipeline(subj)

% Unit test, in case this function is called with no "subj" parameters
if nargin == 0
    subj = initialize_vars('01_vav_p', 'preterm', 'pre');
end
initialPath     = subj.dataDir;
structPath      = subj.structData;
tasks           = subj.tasks;


% Loop through all possible tasks
for thisTask = 1:length(tasks)
    
    try
        functionalPath = [initialPath tasks{thisTask} '/'];
        %preprocess12(functionalPath,structPath,'procChain',{'realignUnwarp','QC','coregister','smooth','segmentDartel'}, 'smoothFWHM',[6 6 6],'coregDirection','FtoS');
        %preprocess12(functionalPath,structPath,'procChain',{'realignUnwarp', 'QC','coregister','smooth'}, 'smoothFWHM',[6 6 6],'coregDirection','FtoS');
        
    catch ME
        msg = ['Preprocessing of sequence ' tasks{thisTask} ' failed for subject ', subj.curSubj];
        causeException = MException('MATLAB:myCode:preprocessing',msg);
        ME = addCause(ME,causeException);
        rethrow(ME)
        continue;
    end
end
OK = 1;
end
