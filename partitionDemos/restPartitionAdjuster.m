%% A function to generate empirically-adjusted resting-state partition
% Requires netStruct.m and associated files 
% 
% OVERVIEW: 
%   - Step 1: order rest FC values according to CA partition
%   - Step 2: find "preferences" in your dataset's resting-state FC values ("preference" = regional max)
%   - Step 3: re-order FC values according to those preferences and save out this adjusted partition 
%   - Step 4: save out rest and task FC matrices, sorted according to the new partition 
%
% INPUTS: 
%   - restFCArray: resting-state FC matrices; 2 forms accepted:
%       1) Cell-array, size = sample size (each cell = 1 subject's data); each cell contains 3D FC matrices of size: [nodes x nodes] = [360 x 360]
%       2) 3D matrix of size: [nodes x nodes x subjects] 
%
%   - taskFCArray: task-evoked FC matrices; 2 forms accepted:
%       1) Cell-array, size = sample size (each cell = 1 subject's data); each cell contains 3D FC matrices of size: [nodes x nodes x tasks]
%       2) 4D matrix of size: [nodes x nodes x tasks x subjects] 
% 
% OUTPUTS: 
%   - taskFCArray_Ordered: a 4D matrix of taskFC values, sorted acccording to the empirically-adjusted CA partition. Size: [nodes x nodes x tasks x subjects] 
%   - restFCArray_Ordered: a 3D matrix of restFC values, sorted according to the empirically adjusted CA partition. Size: [nodes x nodes x subjects] 
%   - nodeOrderNew: a vector that can be used to index CA-sorted regional values 
%   - boundariesNew: matrix of region boundaries for partition scheme that is size: [12 networks x 3], where the 3 columns represent:
%       1) Beginning node number (index of regional start boundary) for each network
%       2) End node number (index of regional end boundary) for each network 
%       3) Network size for each network (e.g., the number of nodes in between the start and end boundaries, per network) 
%       * EXAMPLE (if partition had 178 regions in 5 networks): boundariesNew = [1 6 6; 7 60 54; 61 99 39; 100 155 56; 156 178 23];
%       * Can be used as the input for other functions, ex: bvcAlgorithm.m 
%
% *****IMPORTANT NOTES:******
%   - Assumes regions composing FC matrices = Glasser parcels (or 360 x 360 adjacency matrices; 180 regions per hemisphere; see Glasser, MF., Coalson, TS., 
%     Robinson, EC., Hacker, CD., Harwell, J., Yacoub, E., ? Van Essen, DC. (2016). A multi-modal parcellation of human cerebral cortex.  Nature, 536(7615), 171?178.)
%   - Uses helper script netStruct.m for CA partition info 

%% FUNCTION CALL
function [taskFCArray_Ordered, restFCArray_Ordered, nodeOrderNew, boundariesNew] = restPartitionAdjuster(restFCArray,taskFCArray)

%% Load and pre-allocate
netStruct; 

if iscell(taskFCArray) 
    numSubjects = numel(taskFCArray); sampleMat = taskFCArray{1}; [numNodes,~,numTasks] = size(sampleMat); clear sampleMat;
elseif ~iscell(taskFCArray) 
    [numNodes,~,numTasks,numSubjects] = size(taskFCArray); 
end 

numNetworks = size(boundariesCA,1); 

%% restFC: sort according to CA partition, then find empirical preferences 
restPreferences = NaN(numNodes,numSubjects); 
for subjNum = 1:numSubjects
    if iscell(restFCArray)
        thisSubj = restFCArray{subjNum}; 
    elseif ~iscell(restFCArray) 
        thisSubj = restFCArray(:,:,subjNum); 
    end 
    subjSortedCA = thisSubj(nodeOrder,nodeOrder); 
    for nodeNum = 1:numNodes
        thisNodeVec = subjSortedCA(nodeNum,:); tempVec = NaN(1,numNetworks);
        for netNum = 1:numNetworks
            startNode = boundariesCA(netNum,1); endNode = boundariesCA(netNum,2); tempVec(1,netNum) = tanh(nanmean(atanh(thisNodeVec(startNode:endNode)))); 
        end 
        %[~,prefIdx] = find(tempVec==max(tempVec)); restPreferences(nodeNum,subjNum) = prefIdx; 
        [~,prefIdx] = sort(tempVec,'descend'); restPreferences(nodeNum,subjNum) = prefIdx(1); 
    end 
end 

%% find consensus: if 50% or more of subjects have the preference 
restMode = mode(restPreferences,2); percentAgree = NaN(numNodes,1); 
for nodeNum = 1:numNodes
    thisMode = restMode(nodeNum); thisVec = restPreferences(nodeNum,:); percentAgree(nodeNum,1) = (numel(find(thisVec==thisMode)))/numSubjects;
end 

restConsensus = NaN(numNodes,1); 
for netNum = 1:numNetworks
    startNode = boundariesCA(netNum,1); endNode = boundariesCA(netNum,2);
    for nodeNum = startNode:endNode
        thisNode = percentAgree(nodeNum); 
        if thisNode<0.5 
            restConsensus(nodeNum,1) = netNum; 
        else 
            restConsensus(nodeNum,1) = restMode(nodeNum); 
        end 
    end 
end 

%% generate boundariesNew variable for use in other functions 
[nodeIndicesNew,nodeOrderNew] = sort(restConsensus,'ascend');

boundariesNew = NaN(numNetworks,3); 
for netNum = 1:numNetworks
    finder = find(nodeIndicesNew==netNum); vals = [min(finder) max(finder)]; valVec = min(finder):max(finder); sampSize = numel(valVec);
    boundariesNew(netNum,1:2) = vals; boundariesNew(netNum,3) = sampSize;
end

%% new restFC partition: adjust above by preferences
restFCArray_Ordered = NaN(numNodes,numNodes,numSubjects); 
for subjNum = 1:numSubjects
    if iscell(restFCArray)
        thisSubj = restFCArray{subjNum}; 
    elseif ~iscell(restFCArray) 
        thisSubj = restFCArray(:,:,subjNum); 
    end 
    restCA = thisSubj(nodeOrder,nodeOrder); restAdjusted = restCA(nodeOrderNew,nodeOrderNew); 
    restFCArray_Ordered(:,:,subjNum) = restAdjusted; 
end 

%% apply to taskFC matrices 
taskFCArray_Ordered = NaN(numNodes,numNodes,numTasks,numSubjects); 
for subjNum = 1:numSubjects
    if iscell(taskFCArray)
        thisSubj = taskFCArray{subjNum}; 
    elseif ~iscell(taskFCArray) 
        thisSubj = taskFCArray(:,:,:,subjNum); 
    end 
    for taskNum = 1:numTasks
        thisTask = thisSubj(:,:,taskNum); taskCA = thisTask(nodeOrder,nodeOrder); taskAdjusted = taskCA(nodeOrderNew,nodeOrderNew); 
        taskFCArray_Ordered(:,:,taskNum,subjNum) = taskAdjusted; 
    end 
end 


end 
