%% A script that uploads and organizes useful info regarding the CA partition (C. Cocuzza, 2018) 
% Spronk M., Ji JL., Kulkarni K., Repovs G., Anticevic A., Cole MW. (Preprint) "Mapping the human brain's cortical-subcortical functional network organization". bioRxiv. 

% *****IMPORTANT NOTES:******
%   - Change the baseDir variable for your local machine's path info 
%   - Assumes a 360 vector glasser parcellation (or 360 x 360 adjacency matrices; 180 regions per hemisphere; see Glasser, MF., Coalson, TS., Robinson, EC., Hacker, CD., 
%     Harwell, J., Yacoub, E., ? Van Essen, DC. (2016). A multi-modal parcellation of human cerebral cortex.  Nature, 536(7615), 171?178.)
%   - requires these files: parcel_order_partitionV4.mat; parcel_network_assignments.mat; parcel_network_assignments.txt; 
%
% Useful variables that are generated and uploaded into the workspace: 
%   - boundariesCA: matrix of region boundaries for CA partition scheme that is size: [12 x 3], where 12 = networks, and the 3 columns represent:
%       1) Beginning node number (index of regional start boundary) for each network
%       2) End node number (index of regional end boundary) for each network 
%       3) Network size for each network (e.g., the number of nodes in between the start and end boundaries, per network) 
%       * EXAMPLE (if partition had 178 regions in 5 networks): boundariesNew = [1 6 6; 7 60 54; 61 99 39; 100 155 56; 156 178 23];
%       * This is NOT the empirically-adjusted partition, this is the original CA partition (Spronk et al., 2017) see restPartitionAdjuster.m for more on adjusting 
%   - networkNames: an array of strings containing network names in order (1 through 12); useful for graphing (e.g., axes labels) 
%   - networkColors: a cell array of network colors that can be used in bar graphs, RGB triplet format 

%% Load & pre-allocate 
baseDir = '/Users/carrisacocuzza/Documents/MATLAB/ColeLabNetPartition_v4/'; % CHANGE ME 

load([baseDir 'parcel_order_partitionv4.mat']); nodeOrder = parcel_order_partitionv4; % CA partition node order 
load([baseDir 'parcel_network_assignments.mat']); % CA partition network numbers 
netAssign = csvread([baseDir 'parcel_network_assignments.txt']); % strings corresponding to the above 

netAssignOrdered = netAssign(nodeOrder); nodeIndices = modified_partition_v2; numNetworks = max(nodeIndices);

%% Sorted node/network vector 
endVec = [6 60 99 155 178 201 251 266 343 350 354 360]; sortedNetVec = NaN(360,1);
for netNum = 1:numNetworks
    thisEnder = endVec(netNum); thisNetIdxs = find(nodeIndices==netNum);
    numNodesInNet = numel(thisNetIdxs); thisStarter = thisEnder - (numNodesInNet-1);
    newClusterIdxs = repmat(netNum,[numNodesInNet 1]); sortedNetVec(thisStarter:thisEnder,1) = newClusterIdxs;
end 

%% Network info struct: names and colors 
network(1).name = {'Primary Visual'}; network(1).index = 1; network(1).colors = [0 0 255]; % DIVIDE BY 255 for RGB VALS (done below)
network(2).name = {'Secondary Visual'}; network(2).index = 2; network(2).colors = [100 0 255];
network(3).name = {'Somatomotor'}; network(3).index = 3; network(3).colors = [0 255 255];
network(4).name = {'Cingulo-Opercular'}; network(4).index = 4; network(4).colors = [153 0 153];
network(6).name = {'Language'}; network(6).index = 6; network(6).colors = [0 154 154];
network(9).name = {'Default'}; network(9).index = 9; network(9).colors = [255 0 0];
network(7).name = {'Frontoparietal'}; network(7).index = 7; network(7).colors = [255 255 0];
network(8).name = {'Auditory'}; network(8).index = 8; network(8).colors = [249 61 251];
network(10).name = {'Posterior Multimodal'}; network(10).index = 10; network(10).colors = [177 89 40];
network(5).name = {'Dorsal-attention'}; network(5).index = 5; network(5).colors = [0 255 0];
network(11).name = {'Ventral Multimodal (Unknown1)'}; network(11).index = 11; network(11).colors = [255 156 0];
network(12).name = {'Orbito-Affective (Unknown2)'}; network(12).index = 12; network(12).colors = [65 124 0];

%% Node (regional) boundaries: CA partition 
% Will show upper and lower node number for each, then the sample size per network (e.g., how many nodes are in each network) 
% Note that this is NOT the empirically adjusted resting-state preferences, this must be done separately (see restPartitionAdjuster.m)

boundariesCA = NaN(numNetworks,3);
for netNum = 1:numNetworks
    finder = find(sortedNetVec==netNum); vals = [min(finder) max(finder)]; valVec = min(finder):max(finder); sampSize = numel(valVec);
    boundariesCA(netNum,1:2) = vals; boundariesCA(netNum,3) = sampSize;
end 

%% Create Labels for network graph axes (or to cluster on node graphs; eg rFC matrices)
axesLabels = cell(1,numNetworks); 
for netNum = 1:numNetworks
    if isempty(network(netNum).name)
        axesLabels(:,netNum) = '';
    else
        axesLabels(1,netNum) = network(netNum).name;
    end
end
yAxisLabels = fliplr(axesLabels);

%% Network Names stored in cell array 
networkNames = {'Primary Visual';'Secondary Visual';'Somatomotor';'Cingulo-Opercular';'Dorsal Attention';'Language';'Frontoparietal';'Auditory';...
    'Default';'Posterior Multimodal';'Ventral Multimodal (Unknown 1)';'Oribito-Affective (Unknown 2)'};
%% Proper network colors for other uses
networkColors = cell(numNetworks,4);
for netNum = 1:numNetworks
    name = network(netNum).name; networkColors(netNum,1) = name; colorRGB = network(netNum).colors; trueColor = colorRGB./255;
    networkColors{netNum,2} = trueColor(1); networkColors{netNum,3} = trueColor(2); networkColors{netNum,4} = trueColor(3);
end 
%% Color map example fig to use for dynamic/manual selection of colors if needed
%colorFig = figure(100); get(gcf);
%colorFig.Units = 'centimeters'; colorFig.Position = [1 1 30 20];
%colorFig.PaperUnits = 'normalized';
%xColors = 1:numNetworks; yColors = ones(numNetworks,1);

%for netNum = 1:numNetworks
   % h = bar(xColors(netNum),yColors(netNum)); set(h,'FaceColor',(network(netNum).colors)./255,'facealpha',0.9); hold on;
%end 
%set(gca,'box','off','color','none','xtick',xColors,'xticklabel',axesLabels,'xticklabelrotation',45,'FontName','Avenir','yticklabel',[]); 
%title('Sample Color Map For Network Paritioner Version 4','fontsize',15);


