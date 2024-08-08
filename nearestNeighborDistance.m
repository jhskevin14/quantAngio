function [meanDist, stdDist] = nearestNeighborDistance(centroids)
    % centroids: ?? ?? ?? (Nx2 ??)
    
    numPoints = size(centroids, 1);
    distances = pdist2(centroids, centroids);
    
    % ?? ???? ??? ??
    distances(1:numPoints+1:end) = Inf;
    
    nearestDist = min(distances, [], 2);
    meanDist = mean(nearestDist);
    stdDist = std(nearestDist);
end
