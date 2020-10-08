% This function initialises the subject variables for the batch analysis of
% the Building the Path data.
%__________________________________________________________________________
% Copyright (C) 2016 MIP:Lab, Lorena Freitas

% Lorena Freitas
% $Id: initialize_vars.m 11 2016-12-01 16:26:24F Lorena $


function [b] = initialize_vars(subject,group, intervention)

% SPM info
%--------------------------------------------------------------------------
b.spmDir = fileparts(which('spm')); % path to SPM installation


% Directory information
%--------------------------------------------------------------------------
dataDir = '/media/miplab-nas2/Data2/Nikolina_VAV_RF/';


% Subject information
%--------------------------------------------------------------------------

str = strsplit(deblank(subject), '_');
switch str{3}
    case 'p'
        b.group = 'preterm';
        if nargin < 3
            subfolder = '/pre/';
        else
            subfolder = [filesep intervention filesep];
            b.interv = intervention;
        end
    case 'c'
        subfolder = '/';
        b.interv = '';
        b.group = 'controls';
end


b.curSubj = subject;
b.dataDir = strcat(dataDir, b.curSubj, subfolder); % make data directory subject-specific
b.structData = strcat(dataDir, b.curSubj, subfolder, 't1/');

% test if subject has t1 for the post intervention:
if length(dir(b.structData)) < 3
    b.structData = strcat(dataDir, b.curSubj, '/pre/', 't1/');
end



% Tasks information
%--------------------------------------------------------------------------
b.tasks = {'RealityFiltering1', 'RealityFiltering2'};
b.atlased = 0;
b.interpolate = 0;
b.gaveResponses = NaN; % this property will be set in loadBehaviouralData.m


% Template information
%--------------------------------------------------------------------------
b.templateDate = '';

end
