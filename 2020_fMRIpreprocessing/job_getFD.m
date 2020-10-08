% JOB_GETFD: This function calculates de Framewise Displacement for a
% subject and saves it in the folder given as "FDdir".
% 
% For more information, read Power et al., NeuroImage, 2014
%____________________________________________________________________________
% Copyright (C) 2016 MIP:Lab

% Lorena Freitas
% $Id: job_getFD.m 11 2017-09-08 12:04:37F Lorena $

function FD = job_getFD(b)

if nargin == 0
    b = initialize_vars('01_mind_p', 'preterm', 'pre');
end

% Load Parameters
% ____________________________________________
RAD             = 50; % Use rotations
FDdir           = '/Volumes/EPFL_Lorena/BtP/Data/FramewiseDisplacement/';
initialPath     = b.dataDir;
tasks           = b.tasks;


if isempty(b.interv)
    intervLabel     ='';
else
    intervLabel = ['_' b.interv];
end

% Create directory to save FD information
if ~exist(FDdir, 'dir')
    mkdir(FDdir);
end

for thisTask = 1:length(tasks)
    
    FD_file = strcat('FD_', tasks{thisTask}, '_', b.curSubj, intervLabel, '.mat');
    
    % If there exists a file with the FD values for this subject, load it
    % ___________________________________________________________________
    %if exist( char(strcat(FDdir, FD_file)),'file')
    %    load(char(strcat(FDdir,FD_file)));
    %else
        % Load motion parameter file
        % ____________________________________________
        thisTaskPath = [initialPath tasks{thisTask} '/unwarped/'];
        rp_file=fullfile(thisTaskPath,spm_select('List',thisTaskPath,'^rp_bold.*\.*'));
        rp_param=load(rp_file);
        clearvars rp_file;
        
        % Calculate FD
        % The code below is from Power et al., 2014
        % For more information check https://www.jonathanpower.net/2014-ni-motion-2.html
        % ____________________________________________
        if RAD % USE ROTATIONS
            % Convert rotations (rad) into translations (mm)
            rp_param_rad = [rp_param(:,1:3) rp_param(:,4:6)*(2*RAD*pi/360)]; % from Power's code
            % Compute the Framewise Displacement
            FD = [0;sum(abs(diff(rp_param_rad)),2)];
        else % DO NOT USE ROTATIONS
            FD = [0;sqrt(sum(diff(rp_param(:,1:3)).^2,2))];
        end
        
        % Save FrameWise Displacement for this subject
        % ____________________________________________
        save(char(strcat(FDdir, FD_file)), 'FD');
        
        
    %end
    
    
end



