%% A function to calculate network partition deviation (NPD) and network partition adherence (NPA) 
% Carrisa Cocuzza, August 2018
% 
% INPUTS:
%   - taskFCArray: task-evoked FC matrices; 2 forms accepted:
%       1) Cell-array, size = sample size (each cell = 1 subject's data); each cell contains 3D FC matrices of size: [nodes x nodes x tasks]
%       2) 4D matrix of size: [nodes x nodes x tasks x subjects] 
%
%   - boundariesNew = matrix of region boundaries for your partition scheme that is size: [networks x 3], where the 3 columns represent:
%       1) Beginning node number (index of regional start boundary) for each network
%       2) End node number (index of regional end boundary) for each network 
%       3) Network size for each network (e.g., the number of nodes in between the start and end boundaries, per network) 
%       * EXAMPLE (if partition had 178 regions in 5 networks): boundariesNew = [1 6 6; 7 60 54; 61 99 39; 100 155 56; 156 178 23];
%
%   - conditionTag: a string; will be tagged onto ouput filenames, ex: conditionTag = 'myStudy'
%
%   - visualizeResults: if set to 1, a graph will appear at the end, and save out using the above file tag; NOTE:
%       * netStruct.m is called here, this is for network colors of the CA partition; if using a different 
%       coloring/partition scheme, turn this off and visualize separately 
% 
% OUTPUTS:
%   - netAffinities: node x network matrix of affinities (affinity profiles, region-wise)
%   - clusteredAffinities: network x network matrix of affinities (affinity profiles, network-wise)
%   - adherenceRFs: network x 1 vector of NPA values (frequency of tasks network adheres to partition)
%   - deviationRFs: network x 1 vector of NPD values (frequency of tasks network deviations from the partition)
%   - deviationVecNodes: node x 1 vector of NPD values (frequency of tasks region adheres to partition)
%   - affinityVecNodes: node x 1 vector of NPA values (frequency of tasks region adheres to partition)
%   - nodePrefIdxs: node x task vector showing which network is the region's preference for that given task 
%
% *****IMPORTANT NOTES:******
%   - This assumes nodes in the FC matrices (taskFCArray input) have already been ordered according to your preferred partition; this should
%   correspond with the info in boundariesNew
%
%   - see netStruct.m (https://github.com/ColeLab/controlCartography) to use the CA partition, and potentially adjusting by empirical 
%   resting-state data (if available in your study); assumes the Glasser 2016 regional parcellation scheme (360 nodes across both hemispheres);
%   but number of nodes is determined in-script here, so a different parcellation-->partition scheme can be used as well 
%
%   - this also assumes that FC values are not normalized beforhand (fisher-z transform used in-script to handle regular R values across subjects) 
% 
%% FUNCTION CALL
function [netAffinities,clusteredAffinities,adherenceRFs,deviationRFs, deviationVecNodes,...
    affinityVecNodes,nodePrefIdxs] = npdAlgorithm(taskFCArray,boundariesNew,conditionTag,visualizeResults)

%% Pre-allocate 

if iscell(taskFCArray) 
    numSubjects = numel(taskFCArray); sampleMat = taskFCArray{1}; [numNodes,~,numTasks] = size(sampleMat); clear sampleMat;
elseif ~iscell(taskFCArray) 
    [numNodes,~,numTasks,numSubjects] = size(taskFCArray); 
end 

numNetworks = size(boundariesNew,1); 

%% Mean FC Over Subjects --> [node x node x task] FC matrices

if iscell(taskFCArray)
    bigFC = NaN(numNodes,numNodes,numTasks,numSubjects);
    for subjNum = 1:numSubjects
        thisSubj = taskFCArray{subjNum}; bigFC(:,:,:,subjNum) = thisSubj; % make the cell a 4D matrix
    end
    meanFC = tanh(nanmean(atanh(bigFC),4));
elseif ~iscell(taskFCArray)
    meanFC = tanh(nanmean(atanh(taskFCArray),4));
end

%% Cluster FC vals down to network-level (weighted avg) --> [node x network x task FC matrices]

clusteredFCArray = NaN(numNodes,numNetworks,numTasks); tempNet = NaN(1,numNetworks); tempTask = NaN(numNodes,numNetworks);
for taskNum = 1:numTasks
    thisTask= meanFC(:,:,taskNum); thisTask(logical(eye(size(thisTask)))) = NaN;
    for nodeNum = 1:numNodes
        thisNodeVec = thisTask(nodeNum,:);
        for netNum = 1:numNetworks
            startNode = boundariesNew(netNum,1); endNode = boundariesNew(netNum,2); thisCluster = thisNodeVec(startNode:endNode);
            clusterMean = tanh(nanmean(atanh(thisCluster))); tempNet(1,netNum) = clusterMean;
        end
        tempTask(nodeNum,:) = tempNet;
    end
    clusteredFCArray(:,:,taskNum) = tempTask;
end

%% Find max FC values & indices (e.g., connecting-node number) of net means from above --> [node x task preference matrix]
tempNodeVecPrefs = NaN(numNodes,1); nodePrefValsMeanFirst = NaN(numNodes,numTasks); 
tempNodeVecIdxs = NaN(numNodes,1); nodePrefIdxs = NaN(numNodes,numTasks); % values here will be index of connecting net # for each task

for taskNum = 1:numTasks
    thisTask = clusteredFCArray(:,:,taskNum);
    for nodeNum = 1:numNodes
        thisNodeVec = thisTask(nodeNum,:); [maxVal,maxIdx] = max(thisNodeVec);
        tempNodeVecPrefs(nodeNum,1) = maxVal; tempNodeVecIdxs(nodeNum,1) = maxIdx;
    end
    nodePrefValsMeanFirst(:,taskNum) = tempNodeVecPrefs; nodePrefIdxs(:,taskNum) = tempNodeVecIdxs; 
end

%% Iterating over all networks: test if preference index is a member of that network --> [node x network x task binarized affinity matrix]
maxMembershipsTF = NaN(numNodes,numTasks,numNetworks); maxMembershipsMean = NaN(numNodes,numTasks,numNetworks);
for nodeNum = 1:numNodes
    thisNodeVec = nodePrefIdxs(nodeNum,:);
    for netNum = 1:numNetworks
        memberTF = thisNodeVec==netNum; maxMembershipsTF(nodeNum,:,netNum) = memberTF; % kept as binary
        maxMembershipsMean(nodeNum,:,netNum) = memberTF*netNum; % 1 x netNum = netIdx
    end
end

%% Find affinity scores (relative frequencies or RFs) for all connecting networks --> [node x network relative freq NPA values]
netAffinities = NaN(numNodes,numNetworks);
for netNum = 1:numNetworks
    thisNet = maxMembershipsTF(:,:,netNum);
    for nodeNum = 1:numNodes
        thisNodeVec = thisNet(nodeNum,:); thisNodeSum = sum(thisNodeVec); % Tally
        thisNodeRF = thisNodeSum/numTasks; % RF
        netAffinities(nodeNum,netNum) = thisNodeRF; % Each row should add to 1
    end
end

%% Cluster ROIs into NOIs  --> [network x network relative freq NPA values]
clusteredAffinities = NaN(numNetworks,numNetworks);
for netNum = 1:numNetworks
    startNode = boundariesNew(netNum,1); endNode = boundariesNew(netNum,2); thisCluster = netAffinities(startNode:endNode,:);
    clusterMeanRFs = nanmean(thisCluster); % should still add to 1
    clusteredAffinities(netNum,:) = clusterMeanRFs;
end

% Check that they all add to 100%
for netNum = 1:numNetworks
    thisNet = clusteredAffinities(netNum,:); fullPercentTest = sum(thisNet)*100; disp(['NOI # ' num2str(netNum) ' adds to ' num2str(fullPercentTest) '%']);
end

%%  Pulling out adherence to pre-defined network partition --> [network x 1] NPA vector 
adherenceRFs = diag(clusteredAffinities); 

%% Calculating NPD as 1-NPA --> [network x 1] NPD vector 
deviationRFs = 1-adherenceRFs; 

%% Visualized adherence/fluctuation dichotomy (they are complements) 
if visualizeResults==1
    netStruct; % Used for network colors 
    stackedMat = cat(1,deviationRFs',adherenceRFs');
    stackFig = figure(1); get(gcf); stackFig.Units = 'centimeters'; stackFig.Position = [1,1,27,17];
    stackFig.PaperUnits = 'normalized'; xVec = 1:12; bb = bar(stackedMat','stacked'); legend('Location','best');
    set(bb(1),'displayname','Fluctuation'); set(bb(2),'displayname','Adherence');
    for netNum = 1:numNetworks
        yLocAdhere = 0.95; yLocFluct = 0.2;
        text(netNum,yLocAdhere,['~' num2str(round(stackedMat(2,netNum)*100)) '%'],'color',...
            'k','horizontalalignment','center','fontsize',9,'FontName','Avenir');
        text(netNum,yLocFluct,['~' num2str(round(stackedMat(1,netNum)*100)) '%'],'color',...
            [1 1 0.9],'horizontalalignment','center','fontsize',9,'FontName','Avenir');
    end
    
    ylim([0 1]);xlabel('Network (Version 4 Partitioner)','fontsize',14); ylabel('Network Affinities (out of 100%)','fontsize',14);
    title(['Network Partition Adherence & Fluctuation; ' num2str(numTasks) ' ' conditionTag ' Tasks, (mean of n=' num2str(numSubjects) ')'],'fontsize',13);
    set(gca,'box','off','color','none','xtick',xVec,'xticklabel',axesLabels,'xticklabelrotation',90,'FontName','Avenir','fontsize',12);
    figName = ['npa_npa_barGraph_meanFirst_' conditionTag]; print(stackFig, figName, '-dpng');
end

%% Region-wise NPA and NPD vectors --> [nodes x 1]; useful for workbench 
affinityVecNodes = NaN(numNodes,1);
for netNum = 1:numNetworks
    startNode = boundariesNew(netNum,1); endNode = boundariesNew(netNum,2); thisMat = netAffinities(startNode:endNode,netNum); 
    affinityVecNodes(startNode:endNode,1) = thisMat; 
end 
deviationVecNodes = 1-affinityVecNodes; 

%% Save important variables 
save(['clusteredFCArray_' conditionTag '.mat'],'clusteredFCArray'); 
save(['nodePrefIdxs_' conditionTag '.mat'],'nodePrefIdxs');
save(['maxMembershipsTF_' conditionTag '.mat'],'maxMembershipsTF');
save(['netAffinities_' conditionTag '.mat'],'netAffinities');
save(['clusteredAffinities_' conditionTag '.mat'],'clusteredAffinities'); 
save(['adherenceRFs_' conditionTag '.mat'],'adherenceRFs');
save(['deviationRFs' conditionTag '.mat'],'deviationRFs');
save(['deviationVecNodes_' conditionTag '.mat'],'deviationVecNodes');
save(['affinityVecNodes_' conditionTag '.mat'],'affinityVecNodes'); 

end 



