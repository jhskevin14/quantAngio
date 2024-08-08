function ripleysK = calculateRipleysK(centroids, maxDist, numBins)
    % centroids: ?? ?? ?? (Nx2 ??)
    % maxDist: ?? ??
    % numBins: ?? ?? ?
    
    numPoints = size(centroids, 1);
    distances = pdist2(centroids, centroids);
    
    ripleysK = zeros(numBins, 1);
    binEdges = linspace(0, maxDist, numBins+1);
    
    for i = 1:numBins
        binStart = binEdges(i);
        binEnd = binEdges(i+1);
        inBin = (distances > binStart) & (distances <= binEnd);
        ripleysK(i) = sum(inBin(:)) / numPoints;
    end
    
    ripleysK = cumsum(ripleysK);
end
