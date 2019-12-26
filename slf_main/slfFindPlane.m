function [yMost,coords] = slfFindPlane(fg)
% Find the approximate middle of the set
yStep = 0.5;
yMost = [];
if isempty(yMost)
    % (1) Find the yMost, y coordinate with most fibers (probably in the middle of the tracts)
    yMin = min(cellfun(@(x) min(x(2,:)), fg.fibers));
    yMax = max(cellfun(@(x) max(x(2,:)), fg.fibers));
    
    yCandidates = yMin:yStep:yMax; % A vector of candidate y values in steps of yStep mm
    numNodes = nan(size(yCandidates));
    for ii = 1:length(yCandidates )
        yTmp = yCandidates(ii);
        isNeighbor = cellfun(@(x) any(abs(x(2,:)-yTmp)<=1),fg.fibers);
        numNodes(ii) = sum(isNeighbor);
    end
    idx = find(numNodes == max(numNodes));
    if length(idx)>1
        idx = idx(ceil(length(idx)/2));
    end
    yMost = yCandidates(idx);
end

% Find the closest coordinate in each fiber to yMost
coords = nan(3,length(fg.fibers));

for fI = 1:length(fg.fibers)
    fiberCoords = fg.fibers{fI};
    yDist = abs(fiberCoords(2,:)-yMost);
    minFiber=round(length(fiberCoords)/5); % We don't want the first 5th of the tract (this is just done to ensure that we use coordinates from around the middle of the tract, and not start or end coordinates which may have travelled back to go through the plane a second time)
    maxFiber=round(length(fiberCoords)*4/5); % We don't want the first 5th of the tract (this is just done to ensure that we use coordinates from around the middle of the tract, and not start or end coordinates which may have travelled back to go through the plane a second time)
    yDist(1:minFiber)=inf;
    yDist(maxFiber:length(fiberCoords))=inf;
    idx = find(yDist < 2);
    if isempty(idx)
        continue % In case no close node was found, this will leave NaNs in the coordiantes of this streamline
    end
    zCoords = fiberCoords(3,idx); % Of all the close nodes, take the one with the highest z coordinate (this is pretty arbitrary)
    [~,zIdx] = max(zCoords);
    coords(:,fI) = fiberCoords(:,idx(zIdx));
end