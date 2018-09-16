%% A script to generate the variables used in the cartographic method, as well as the graphs 
% "Network Science Cartography of Cogntiive Control System Dynamics", C. Cocuzza et al., 2018 

% Sample data used here is n=2 for GitHub file size limits; original study used n=50 
% boundariesNew and nodeOrderNew were previously calculated for full dataset with restPartitionAdjuster.m 
% GVC variables were previously calculated and are loaded; see http://www.colelab.org/cole-etal-2013/#analysiscode
% NPD variables were previously calculated and are loaded; see npdAlgorithm.m
% Requires netStruct.m and related variables (see script) 

%% Load & pre-allocate
% CHANGE baseDir
baseDir = '/Users/carrisacocuzza/Desktop/qual_scripts/'; 

load([baseDir 'restFC_2Subjs.mat']); load([baseDir 'taskFC_2Subjs.mat']); 
load([baseDir 'boundariesNew.mat']); load([baseDir 'nodeOrderNew.mat']); 
load([baseDir 'gvcNetsNormed2Subjs.mat']); load([baseDir 'gvcNetsSubjectsNormed2Subjs.mat']); 
load([baseDir 'npdNetsZ_2Subjs.mat']); load([baseDir 'nodePrefIdxs_2Subjs.mat']); 
load([baseDir 'deviationVecNodes_2Subjs.mat']); load([baseDir 'restConsensus.mat']); 

[numNodes,~,numTasks,numSubjs] = size(ruleFC2); netStruct; 

%% Order regions in FC matrices 
restOrdered = NaN(numNodes,numNodes,numSubjs); taskOrdered = NaN(numNodes,numNodes,numTasks,numSubjs); 
for subjNum = 1:numSubjs
    subjRest = restFC2(:,:,subjNum); subjTask = ruleFC2(:,:,:,subjNum); 
    restCA = subjRest(nodeOrder,nodeOrder); restAdj = restCA(nodeOrderNew,nodeOrderNew); 
    for taskNum = 1:numTasks
        thisTask = subjTask(:,:,taskNum); taskCA = thisTask(nodeOrder,nodeOrder); taskAdj = taskCA(nodeOrderNew,nodeOrderNew);
        taskOrdered(:,:,taskNum,subjNum) = taskAdj; 
    end 
    restOrdered(:,:,subjNum) = restAdj; 
end 

%% Grand means of FC matrices 
taskMeansMid = tanh(nanmean(atanh(taskOrdered),4)); taskMeans = tanh(nanmean(atanh(taskMeansMid),3)); taskMeans(logical(eye(size(taskMeans)))) = NaN;
restMeans = tanh(nanmean(atanh(restOrdered),3)); restMeans(logical(eye(size(restMeans)))) = NaN; 

%% Global variability - between-network cartography of GVC-Task vs GVC-Region %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GVC-Region (same as GVC - which is done across tasks - but across regions)
tempMat = NaN(numNodes,numTasks); subjNodeGVC = NaN(numNodes,numSubjs); 
for subjNum = 1:numSubjs
    thisSubj = taskOrdered(:,:,:,subjNum); 
    for taskNum = 1:numTasks
        thisTask = thisSubj(:,:,taskNum); nodeStd = tanh(nanstd(atanh(thisTask),0,2)); tempMat(:,taskNum) = nodeStd; 
    end 
    subjNodeGVC(:,subjNum) = nanmean(tempMat,2); 
end 

% Z-score --> Cluster & SEM
subjMeanVector = nanmean(subjNodeGVC); subjMeans = repmat(subjMeanVector,numNodes,1); subjSDVector = nanstd(subjNodeGVC); subjSDs = repmat(subjSDVector,numNodes,1);
nodeGVCNormedSubjs = (subjNodeGVC - subjMeans)./subjSDs; 

nodeGVCnetMeans = NaN(numNetworks,1); nodeGVCnetSEMs = NaN(numNetworks,1); 
for netNum=1:numNetworks
    startNode = boundariesNew(netNum,1); endNode = boundariesNew(netNum,2); thisVal = nodeGVCNormedSubjs(startNode:endNode,:); 
    nodeGVCnetMeans(netNum,1) = nanmean(nanmean(thisVal));  nodeGVCnetSEMs(netNum,1) = (nanstd(nanmean(thisVal)))/(sqrt(numSubjs)); 
end 

%% Visualize - GVC-region vs GVC-Task Cartography
figH = figure(1); get(gcf); figH.Units = 'centimeters'; figH.Position = [1,1,30,22]; figH.PaperUnits = 'normalized';
%netNames = {'Vis1','Vis2','Smtr','CON','DAN','LAN','FPN','Aud','DMN','PMM','VMM','OAN'};

for netNum = 1:numNetworks
    s1=scatter(nodeGVCnetMeans(netNum),gvcNetsNormed(netNum),'filled'); thisCol = cell2mat(networkColors(netNum,2:end)); s1.MarkerFaceColor = thisCol;
    s1.Marker = 'd'; s1.SizeData = 250; s1.MarkerEdgeColor = 'k'; s1.MarkerEdgeAlpha = 0.1; 
    s1.MarkerFaceAlpha = 0.95;  hold on;
end 
xLab = xlabel({'Across-Region Variability (GVC-Region, z-scored)'}); yLab = ylabel({'Across-Task Variability (GVC-Task, z-scored)'});
smallTitle = title('Global Variability, Across 2 Dimensions'); xlim([-1.2 1.2]); ylim([-1.2 1.2]); hold on;
minX = min(nodeGVCnetMeans); maxX = max(nodeGVCnetMeans); minY = min(gvcNetsNormed); maxY = max(gvcNetsNormed);

lH = line([minX-1 maxX+1],[0 0]); lH.LineStyle = ':'; lH.LineWidth = 1.5; lH.Color = [0.2 0.2 0.2]; 
lH = line([0 0],[minY-1 maxY+1]); lH.LineStyle = ':'; lH.LineWidth = 1.5; lH.Color = [0.2 0.2 0.2]; 

set(gca,'color','none','box','off','fontname','avenir','fontsize',20); set(xLab,'fontsize',22); set(yLab,'fontsize',22); set(smallTitle,'fontsize',25); 

figName = 'global_cartography'; print(figH, figName, '-dpng'); 

%% Local, within-network cartography of partition integrity vs partition deviation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Integrity: jaccard index over task states  
jaccardRestTask = NaN(numNetworks,numTasks); 
for taskNum = 1:numTasks
    taskPref = nodePrefIdxs(:,taskNum); jaccardRestTask(:,taskNum) = jaccard(restConsensus,taskPref); 
end 
jaccardMean = nanmean(jaccardRestTask,2); jaccardNorm = zscore(jaccardMean); 

%% Visualize - NPD vs Jaccard Cartography - no labels
figH = figure(2); get(gcf); figH.Units = 'centimeters'; figH.Position = [1,1,30,22]; figH.PaperUnits = 'normalized';
netNames = {'Vis1','Vis2','Smtr','CON','DAN','LAN','FPN','Aud','DMN','PMM','VMM','OAN'};

for netNum = 1:numNetworks
    s1=scatter(npdNetsZ(netNum),jaccardNorm(netNum),'filled'); thisCol = cell2mat(networkColors(netNum,2:end)); s1.MarkerFaceColor = thisCol;
    s1.Marker = 'd'; s1.SizeData = 250; s1.MarkerEdgeColor = 'k'; s1.MarkerEdgeAlpha = 0.1; s1.MarkerFaceAlpha = 0.95; hold on;
end 
xLab = xlabel({'Network Partition Deviation (z-scored)'}); yLab = ylabel({'Network Partition Integrity (Jaccard, z-scored)'});
smallTitle = title('Local Integrity & Deviation'); xlim([-1.6 1.6]); ylim([-1.6 1.6]); hold on;

minX = min(npdNetsZ); maxX = max(npdNetsZ); minY = min(jaccardNorm); maxY = max(jaccardNorm);
lH = line([minY-1 maxY+1],[0 0]); lH.LineStyle = ':'; lH.LineWidth = 1.5; lH.Color = [0.2 0.2 0.2]; 
lH = line([0 0],[minY-1 maxY+1]); lH.LineStyle = ':'; lH.LineWidth = 1.5; lH.Color = [0.2 0.2 0.2]; 

set(gca,'color','none','box','off','fontname','avenir','fontsize',20); set(xLab,'fontsize',22); set(yLab,'fontsize',22); set(smallTitle,'fontsize',25); 

figName = 'local_cartography'; print(figH, figName, '-dpng'); 




