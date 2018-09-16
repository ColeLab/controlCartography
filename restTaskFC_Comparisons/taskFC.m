%% A function to estimate resting-state functional connectivity with Pearson's correlation coefficient 
% Carrisa Cocuzza, August 2018

%% Function specifications:
% Need to be within the directory your files will be pulled from and you'd like output to be saved to (use cd to change directory) 
% 
% This function assumes pre-processing, parcellation, nuisance regression, and task-regression has already taken place (referenced study used FIR task regression)
% This function assumes task-activity design matrices have already been generated (via task regression)
%
% INPUTS: 
%   1) taskDataFile = a string; filename of task residuals (ex: taskDataFile = 'taskData_MyStudy.mat';)
%       - two formats are accepted:
%       - 1) a cell array the size of n subjects (ex: N=50, cell size = 1x50); each cell contains a region x TR matrix of task residuals 
%       - 2) a 3D matrix of task residuals, of size: [regions x TRs x subjects]
%
%   2) taskDesignFile = a string; filename of task design matrices (ex: taskDesignFile = 'taskDesign_MyStudy.mat')
%       - two formats are accepted:
%       - 1) a cell array the size of n subjects (ex: N=50, cell size = 1x50); each cell contains a TR x task number matrix of binarized task design locations 
%       - 2) a 3D matrix of task design locations, of size: [TR x task number x subjects]
%
%   3) runTag = a string; this will be used to tag output filenames
%
% OUTPUTS: 
%   all outputs are saved out files; please load them back into MATLAB for further analyses; output files include:
%       - taskFCArray_<runTag>.mat = task-state FC estimates, sample size (N) cell array, each cell of size [regions x regions x task number] 
%       - taskFCpVals_<runTag>.mat = non-corrected p-values associated with FC estimates, of size [regions x regions x task number x subjects] 
%
% NOTES:
%       - If transforming data in later analyses (e.g., averaging across subjects), use the fisher z-transform via atanh --> the operation --> then tanh
%       - Regions are not ordered in this version; see netStruct_CA.m for the CA partition ordering variables
%
% EXAMPLE (using files saved in the GitHub repository):
%       - taskDataFile = 'firTaskData.mat'; taskDesignFile = 'firTaskDesign.mat'; runTag = 'TestRun'; 
%       - taskFC(taskDataFile,taskDesignFile,runTag); 
%
% USAGE NOTES:
%       - 50 subjects should take ~ 5-10 minutes on modern operating systems 
%% FUNCTION 
function taskFC(taskDataFile,taskDesignFile,runTag)

%% LOAD & CONFIGURE VARIABLES 
loadData = load(taskDataFile); loadDesign = load(taskDesignFile); 
taskCell = struct2cell(loadData); testCell = taskCell{1,1};
designCell = struct2cell(loadDesign); 

% note: assumes that formatting of taskDataFile and taskDesignFile are consistent (e.g., if one is a cell array, so is the other) 
if iscell(testCell)
    taskData = taskCell{:}; numSubjs = numel(taskData); sampleSubj = taskData{1}; [numRegions,numTRorig] = size(sampleSubj); 
    taskDesign = designCell{:}; sampleDesign = taskDesign{1}; [numTR,numTasks] = size(sampleDesign); 
elseif ~iscell(testCell)
    taskData = cell2mat(taskCell); [numRegions,numTRorig,numSubjs] = size(taskData); 
    taskDesign = cell2mat(designCell); [numTR,numTasks,~] = size(taskDesign); 
end 

if numTR~=numTRorig
    error('The number of TRs in the task-design file and task-data file are not equal, please check.'); 
end

%% Use task design matrix to index task data 
taskMatrix4D = cell(1,numSubjs);

for subjNum = 1:numSubjs
    
    if numSubjs>1
        thisSubjTaskTS = taskDesign{subjNum};
    elseif numSubjs==1
        thisSubjTaskTS = taskDesign; 
    end
    taskTSidx = thisSubjTaskTS >= 0.5; taskResids = taskData{subjNum}; % only use activations above 0.5
    taskResids3D = repmat(taskResids,1,1,numTasks); taskIdx3D = zeros(numRegions,numTR,numTasks);
    
    for taskNum = 1:numTasks
        thisTaskCol = taskTSidx(:,taskNum); thisTaskPageIdx = repmat(thisTaskCol',numRegions,1); taskIdx3D(:,:,taskNum) = thisTaskPageIdx;
    end
    
    taskArrayNaN = taskIdx3D.*taskResids3D; taskArrayNaN(taskArrayNaN==0) = NaN;
    taskMatrix4D{subjNum} = taskArrayNaN; %sample-size based cell array with [region x TR x task number] matrices
end

%% Pearson R calculations for task FC (looped over subjects)
taskFCArray = cell(1,numSubjs); netMatTemp = NaN(numRegions,numRegions,numTasks);
taskFCpVals = cell(1,numSubjs); pMatTemp = NaN(numRegions,numRegions,numTasks); 

for subjNum = 1:numSubjs
    thisSubject = taskMatrix4D{subjNum};
    
    for taskNum = 1:numTasks
        thisTaskTS = thisSubject(:,:,taskNum); [netMatEachTask,pVal,~,~] = corrcoef(thisTaskTS','rows','complete');
        netMatEachTask(logical(eye(size(netMatEachTask)))) = NaN; netMatTemp(:,:,taskNum) = netMatEachTask;
        pMatTemp(:,:,taskNum) = pVal; 
    end
    taskFCArray{subjNum} = netMatTemp; taskFCpVals{subjNum} = pMatTemp; 
    
end

%% Save out results 
save(['taskFCArray_' runTag '.mat'],'taskFCArray','-v7.3'); save(['taskFCpVals_' runTag '.mat'],'taskFCpVals','-v7.3');

end 




