function allFDValues = generateFD_groupLevel
% This function joins the FDs for all subjects, to create a plot of all of them

FDdir = '/INSERT HERE THE PARTH TO WHERE YOU WANT TO SAVE THE FILE allFDs.mat/FramewiseDisplacement/';
allFDsFile = 'allFDs.mat';
M_all = nan(1,30); I_all = nan(1,30); V_all = nan(1,30);
tasks           = { 'RealityFiltering1', 'RealityFiltering2', 'Rest'};




% If there exists a file with the FD values for this subject, load it
if exist( char(strcat(FDdir, allFDsFile)),'file')
    load(char(strcat(FDdir, allFDsFile)));
    warning('AllFDsFile already exists. This will overwrite it!');
    delete(char(strcat(FDdir, allFDsFile)));
else
    
    allFDdata =  dir([FDdir  '*.mat']);
    allFDValues = nan(1000,length(allFDdata));
    
    for i = 1:length(allFDdata)
        thisFD = load([ FDdir allFDdata(i).name]);
        padding = 1000-length(thisFD.FD); % if run is shorter than 310, pad it so it fits in the matrix
        thisFD = padarray(thisFD.FD, [padding, 0], nan, 'post');
        allFDValues(:,i) = thisFD(:);
        
       
        [M,V] = regexp(sprintf('%i',[0 diff((thisFD>1)')==0]),'1+','match');
        [M,I] = max(cellfun('length',M));
        M_all(i) = M + 1;  % Holds the max length.
        I_all(i) = V(I) - 1;  % Holds the starting index of first run of length M.
        V_all(i) = thisFD(I);  % The value of the run.
    end
    
    % Save FrameWise Displacement for all subjects
    % ____________________________________________
    save(char(strcat(FDdir, allFDsFile)), 'allFDValues');
    
end

% Plot FD vs Subject
% ____________________________________________
figure;
imagesc((allFDValues>1)');
title('Framewise displacement in the Language IRM dataset');
xlabel('Volumes');
ylabel('Subjects');
legend('test', 'test');

bla = allFDValues>1;



end