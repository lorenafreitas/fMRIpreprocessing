% -------------------------------------------------------
%
%    ppCheckMovement - checking if subject has moved too much
%
%    Ver. 1.0.0
%
%    Created:         Daniela Zoller      (21.10.2015)
%    Last modified:   Daniela Zoller      (21.10.2015)
%
%
% ------------------------------------------------------
%
% movementOk = ppCkeckMovement(realignedDataFile)
% 
% reads transformation matrix stored after realignment and verifies if 
% the values are less than 3mm
%
% Input: realignedFataFile: full path to text file containing tranformation
%                           data (translation + rotation)
% Output: movementOk: - 1 if movement parameters are fine
%                     - 0 if parameters are too large
%

function movementOk = ppCheckMovement(realignedDataFile)

movementOk = 1;

%% import data
transfMatr = importdata(realignedDataFile);

%% check max translation 3mm
maxTransl = max(max(transfMatr(:,1:3)));
fprintf('max translation : %.2f mm\n',maxTransl);
if maxTransl > 3
    fprintf('More than 3 mm translation\n');
    movementOk = 0;
    return;
elseif maxTransl > 2.5
    warning(['Be careful, maximum translation: ',sprintf('%.2f mm',maxTransl)]);
end


%% check max rotation 3?
% transform to degree
angles = rad2deg(transfMatr(:,4:6));
maxAngle = max(max(angles));
fprintf('max rotation : %.2f degree\n',maxAngle);
if maxAngle > 3
    fprintf('More than 3 deg rotation\n');
    movementOk = 0;
    return;
elseif maxAngle > 2.5
   warning(['Be careful, maximum rotation: ',sprintf('%.2f degree',maxAngle)]);
end 





