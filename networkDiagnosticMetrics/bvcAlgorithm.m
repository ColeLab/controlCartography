function [bvcNodes, bvcNodesNormed, bvcNets, bvcNetsNormed, netBVCsemsNormed,bvcNodesSubjects,...
    bvcNetsSubjectsNormed] = bvcAlgorithm(taskFCArray,boundariesNew,conditionTag,visualizeResults)

%% A function to calculate BVC; saves out a few vectors/matrices for different purposes
% Carrisa Cocuzza, December 2017

% INPUTS:
%   - taskFCArray: task-evoked FC matrices; 2 forms accepted:
%       1) Cell-array, size = sample size (each cell = 1 subject's data); each cell contains 3D FC matrices of size: [nodes x nodes x tasks]
%       2) 4D matrix of size: [nodes x nodes x tasks x subjects] 
%       - This assumes the nodes have already been sorted by preferred partition (described by boundariesNew variable, see next)
%
%   - boundariesNew: a matrix of region boundaries for your partition scheme that is size: [networks x 3], where the 3 columns represent:
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
% OUTPUTS (use ~ to silence):
%   - bvcNodes: node x 1 vector summarizing each node's average BVC value (across subjects) 
%   - bvcNodesNormed: same as above, normalized
%   - bvcNets: net x 1 vector summarizing each brain network's average BVC value (averaged across node's comprising each network, and across subjects)
%   - bvcNetsNormed: same as above, normalized 
%   - netBVCsemsNormed: standard error of the mean for above 
%   - bvcNodesSubjects: node x subject matrix of regional BVC scores; i.e., bvcNodes before averaging across subjects 
%   - bvcNetsSubjectsNormed: same as above, normalized 
%
% *****IMPORTANT NOTES:******
%   - This assumes nodes in the FC matrices (taskFCArray input) have already been ordered according to your preferred partition; this should
%   correspond with the info in boundariesNew
%
%   - see netStructV4.m (https://github.com/ColeLab/controlCartography) to use the CA partition, and potentially adjusting by empirical 
%   resting-state data (if available in your study); assumes the Glasser 2016 regional parcellation scheme (360 nodes across both hemispheres);
%   but number of nodes is determined in-script here, so a different parcellation-->partition scheme can be used as well 
%
%   - this also assumes that FC values are not normalized beforhand (fisher-z transform used in-script to handle regular R values across subjects) 
% 
%% Pre-allocate 

if iscell(taskFCArray) 
    numSubjects = numel(taskFCArray); sampleMat = taskFCArray{1}; [numNodes,~,numTasks] = size(sampleMat); clear sampleMat;
elseif ~iscell(taskFCArray) 
    [numNodes,~,numTasks,numSubjects] = size(taskFCArray); 
end 

numNetworks = size(boundariesNew,1); 

%% Calculate GVC ("VC" per region)
% See: http://www.colelab.org/cole-etal-2013/#analysiscode
identityMat = eye(numNodes)~=0; identityMatFull = repmat(identityMat,[1 1 numTasks]); GVCmat = NaN(numNodes,numNodes,numSubjects);
for subjNum = 1:numSubjects
    if iscell(taskFCArray)
        thisSubject = taskFCArray{subjNum};
    elseif ~iscell(taskFCArray)
        thisSubject = taskFCArray(:,:,:,subjNum);
    end
    thisSubject(identityMatFull) = NaN; GVCmat(:,:,subjNum) = nanstd(thisSubject,0,3);
end

%% Cluster order & NaN-out within-network vals 

outNetMat = NaN(numNodes,numNodes,numSubjects);
for subjNum = 1:numSubjects
    thisSubjectReg = GVCmat(:,:,subjNum); thisSubjectReg(identityMat) = NaN; tempSubj = NaN(numNodes,numNodes);
    for netNum = 1:numNetworks
        endNode = boundariesNew(netNum,2); startNode = boundariesNew(netNum,1); 
        numNodesInNet = boundariesNew(netNum,3); tempMat = NaN(numNodesInNet,numNodes);
        for nodeNum = startNode:endNode
            thisNodeVec = thisSubjectReg(nodeNum,:); nilNodeStarter = startNode+1; thisNodeVec(nilNodeStarter:endNode) = NaN;
            startIdx = nodeNum-startNode+1; tempMat(startIdx,:) = thisNodeVec;
        end
        tempSubj(startNode:endNode,:) = tempMat;
    end
    outNetMat(:,:,subjNum) = tempSubj;
end

%% Between-network means (BVC)

bvcNodesSubjects = NaN(numNodes,numSubjects);
for subjNum = 1:numSubjects
    thisSubj = outNetMat(:,:,subjNum); regionMeans = nanmean(thisSubj,2); bvcNodesSubjects(:,subjNum) = regionMeans;
end
bvcNodes = nanmean(bvcNodesSubjects,2); 

%% Normalize

subjMeanVector = nanmean(bvcNodesSubjects); subjMeans = repmat(subjMeanVector,numNodes,1); 
subjSDVector = nanstd(bvcNodesSubjects); subjSDs = repmat(subjSDVector,numNodes,1);
bvcNodesSubjectsNormed = (bvcNodesSubjects - subjMeans)./subjSDs; 
bvcNodesNormed = nanmean(bvcNodesSubjectsNormed,2);

%% Cluster into networks: NORMALIZED

bvcNetsSubjectsNormed = NaN(numNetworks,numSubjects);
for netNum = 1:numNetworks
    startNode = boundariesNew(netNum,1); endNode = boundariesNew(netNum,2); thisCluster = bvcNodesSubjectsNormed(startNode:endNode,:);
    clusterMeans = nanmean(thisCluster); bvcNetsSubjectsNormed(netNum,:) = clusterMeans;
end

%% Average across subjects & get SEMS: NORMALIZED

bvcNetsNormed = NaN(numNetworks,1); netBGVCstdsNormed = NaN(numNetworks,1); netBVCsemsNormed = NaN(numNetworks,1);

for netNum = 1:numNetworks
    thisNet = bvcNetsSubjectsNormed(netNum,:); thisNetMean = nanmean(thisNet); thisNetStd = nanstd(thisNet);
    thisNetSEM = thisNetStd/(sqrt(numSubjects)); bvcNetsNormed(netNum,1) = thisNetMean;
    netBGVCstdsNormed(netNum,1) = thisNetStd; netBVCsemsNormed(netNum,1) = thisNetSEM;
end

%% Network bar graph: NORMALIZED VERSION
if visualizeResults==1
    netStruct; % Call up network colors/names 
    barFig = figure(1); get(gcf); barFig.Units = 'centimeters'; barFig.Position = [1,1,30,22];
    barFig.PaperUnits = 'normalized'; hold on; xVec = 1:numNetworks; yVecBGVC = bvcNetsNormed;
    minBVC = min(bvcNetsNormed); minLoc = find(bvcNetsNormed==min(bvcNetsNormed)); minVal = (minBVC - netBVCsemsNormed(minLoc)) - 0.1;
    maxBVC = max(bvcNetsNormed); maxLoc = find(bvcNetsNormed==max(bvcNetsNormed)); maxVal = (maxBVC + netBVCsemsNormed(maxLoc)) + 0.1;
    
    for netNum = 1:numNetworks
        h = bar(xVec(netNum),yVecBGVC(netNum));
        thisColor = [networkColors{netNum,2} networkColors{netNum,3} networkColors{netNum,4}];
        set(h,'FaceColor',thisColor,'facealpha',0.7); hold on;
    end
    
    ylim([minVal maxVal]); xlabel('Network','fontsize',16); ylabel('Mean BVC (z-scored)','fontsize',16);
    eBGVC = errorbar(1:numNetworks,yVecBGVC,netBVCsemsNormed,'.','linewidth',0.6); eBGVC.Color = 'k'; eBGVC.Marker = '.';
    set(gca,'box','off','color','none','xtick',xVec,'xticklabel',axesLabels,'xticklabelrotation',...
        90,'FontName','Avenir','fontsize',12);
    title('Mean BVC: Normalized (z-scored)','FontSize',18);
    
    figName = ['bvcBarGraph' conditionTag '_Normalized']; print(barFig, figName, '-dpng');
end

%% Cluster into networks: NON-NORMALIZED

bvcNetsSubjects = NaN(numNetworks,numSubjects);
for netNum = 1:numNetworks
    startNode = boundariesNew(netNum,1); endNode = boundariesNew(netNum,2); thisCluster = bvcNodesSubjects(startNode:endNode,:);
    clusterMeans = nanmean(thisCluster); bvcNetsSubjects(netNum,:) = clusterMeans;
end

%% Average across subjects & get SEMSs: NON-NORMALIZED

bvcNets = NaN(numNetworks,1); netBGVCstds = NaN(numNetworks,1); netBGVCsems = NaN(numNetworks,1);

for netNum = 1:numNetworks
    thisNet = bvcNetsSubjects(netNum,:); thisNetMean = nanmean(thisNet);
    thisNetStd = nanstd(thisNet); thisNetSEM = thisNetStd/(sqrt(numSubjects));
    bvcNets(netNum,1) = thisNetMean; netBGVCstds(netNum,1) = thisNetStd;
    netBGVCsems(netNum,1) = thisNetSEM;
end

%% Network bar graph: NON-NORMALIZED VERSION
    
if visualizeResults==1
    barFig2 = figure(2); get(gcf); barFig2.Units = 'centimeters'; barFig2.Position = [1,1,30,22];
    barFig2.PaperUnits = 'normalized'; hold on; xVec = 1:numNetworks; yVecBGVC = bvcNets;

    for netNum = 1:numNetworks
        h = bar(xVec(netNum),yVecBGVC(netNum));
        thisColor = [networkColors{netNum,2} networkColors{netNum,3} networkColors{netNum,4}];
        set(h,'FaceColor',thisColor,'facealpha',0.7); hold on;
    end
    
    maxBVCReg = max(bvcNets); maxLocReg = find(bvcNets==max(bvcNets)); maxValReg = (maxBVCReg + netBGVCsems(maxLocReg)) + 0.1;
    
    ylim([0 maxValReg]); 
    xlabel('Network','fontsize',16); ylabel('Mean BVC','fontsize',16);
    eBGVC = errorbar(1:numNetworks,yVecBGVC,netBGVCsems,'.','linewidth',0.6); eBGVC.Color = 'k'; eBGVC.Marker = '.';
    set(gca,'box','off','color','none','xtick',xVec,'xticklabel',axesLabels,'xticklabelrotation',...
        90,'FontName','Avenir','fontsize',12);
    title('Mean BVC: non-normalized','FontSize',18);
    
    figName = ['bvcBarGraph' conditionTag]; print(barFig2, figName, '-dpng');
end 

end 