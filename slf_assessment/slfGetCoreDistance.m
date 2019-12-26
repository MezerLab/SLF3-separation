function dist = slfGetCoreDistance(fg1,fg2,method)
% Return the mean Euclidean distance betwen the tract cores of fg1 and fg2
% 'method' can be either 'closestnode' or 'nodebynode'.

nodesNum = 100;
core1 = dtiComputeSuperFiberRepresentation(fg1,[],nodesNum);
core2 = dtiComputeSuperFiberRepresentation(fg2,[],nodesNum);

core1 = fgFlip(core1,2,'min');
core2 = fgFlip(core2,2,'min');

if strcmp(method,'closestnode')
    % For each node, return its distance from the closest node to it in the
    % other core fg
    dist = pdist2(core1.fibers{1}',core2.fibers{1}','Euclidean','Smallest',1);
elseif strcmp(method,'nodebynode')
    % This will compute the pairwise distance between nodes, but it might be
    % problematic if one tracts starts a little far from the other, and then
    % all nodes are shifted and distances are increased
    dist = nan(1,nodesNum);
    for nI = 1:nodesNum
        dist(nI) = norm(core1.fibers{1}(:,nI)-core2.fibers{1}(:,nI));
    end
end

dist = mean(dist);