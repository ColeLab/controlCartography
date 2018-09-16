%% A function to estimate resting-state functional connectivity with Pearson's correlation coefficient 
% Carrisa Cocuzza, August 2018

%% Function specifications:
% You need to be within the directory your files will be pulled from and you'd like output to be saved to (use cd to change directory) 
% 
% This function assumes pre-processing, parcellation, and nuisance regression has already taken place
%
% INPUTS: 
%   1) restDataFile = a string; filename of resting-state residuals (ex: restDataFile = 'restData_MyStudy.mat';)
%       - two formats are accepted:
%       - 1) a cell array the size of n subjects (ex: N=50, cell size = 1x50); each cell contains a region x TR matrix of resting-state residuals 
%       - 2) a 3D matrix of resting-state residuals, of size: [regions x TRs x subjects]
%
%   2) runTag = a string; this will be used to tag output filenames
%
% OUTPUTS: 
%   all outputs are saved out files; please load them back into MATLAB for further analyses; output files include:
%       - restFCArray_<runTag>.mat = resting-state FC estimates, of size [regions x regions x subjects] 
%       - restFCpVals_<runTag>.mat = p-values associated with FC estimates, of size [regions x regions x subjects] 
%
% NOTES:
%       - If transforming data in later analyses (e.g., averaging across subjects), use the fisher z-transform via atanh --> the operation --> then tanh
%       - Regions are not ordered in this version; see netStruct_CA.m for the CA partition ordering variables
%
% EXAMPLE (using files saved in the GitHub repository):
%       - restDataFile = 'firRestData.mat'; runTag = 'TestRun'; 
%       - restFC(restDataFile,runTag); 
%% rsFC function
function restFC(restDataFile,runTag)

%% LOAD & CONFIGURE VARIABLES 
loadData = load(restDataFile); restCell = struct2cell(loadData); testCell = restCell{1,1}; 
if iscell(testCell)
    restData = restCell{:}; numSubjs = numel(restData); sampleSubj = restData{1};
    [numRegions,~] = size(sampleSubj);
elseif ~iscell(testCell)
    restData = cell2mat(restCell); [numRegions,~,numSubjs] = size(restData); 
end 

%% PERFORM FC ESTIMATION WITH PEARSON'S CORRELATION COEFFICIENT 
restFCpVals = NaN(numRegions,numRegions,numSubjs); restFCArray = NaN(numRegions,numRegions,numSubjs); 

for subjNum = 1:numSubjs
    if iscell(testCell)
        thisSubjsData = restData{subjNum}; 
    elseif ~iscell(testCell) 
        thisSubjsData = restData(:,:,subjNum);
    end 
    [corrMatSubj,pVal] = corrcoef(thisSubjsData','rows','complete');
    restFCpVals(:,:,subjNum) = pVal; corrMatSubj(logical(eye(size(corrMatSubj)))) = NaN; % NaN out diagonal
    restFCArray(:,:,subjNum) = corrMatSubj;
end 

%% SAVE OUTPUTS 
save(['restFCArray_' runTag '.mat'],'restFCArray'); 
save(['restFCpVals_' runTag '.mat'],'restFCpVals'); 

end 





