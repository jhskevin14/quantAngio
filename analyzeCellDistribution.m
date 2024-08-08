%% ?? ???
% imagePath = 'example_image.png'; % ??? ?? ??
% method = 'ripley'; % ?? ?? ??: 'ripley', 'nearest', 'frequency'
% analyzeCellDistribution(imagePath, method);

%%
function analyzeCellDistribution(image, method)
    % imagePath: ??? ?? ??
    % method: ?? ?? ('ripley', 'nearest', 'frequency')
    
    % ??? ????

    % image = imbinarize(image); % ???
    
    % ?? ?? ?? ??
    props = regionprops(image, 'Centroid');
    centroids = cat(1, props.Centroid);
    
    switch method
        case 'ripley'
            maxDist = 3000; % ?? ?? ??
            numBins = 200; % ?? ? ?
            ripleysK = calculateRipleysK(centroids, maxDist, numBins);
            figure;
            plot(linspace(0, maxDist, numBins), ripleysK);
            title('Ripley''s K Function');
            xlabel('Distance');
            ylabel('K');
        
        case 'nearest'
            [meanDist, stdDist] = nearestNeighborDistance(centroids);
            fprintf('Mean Nearest Neighbor Distance: %.2f\n', meanDist);
            fprintf('STD of Nearest Neighbor Distance: %.2f\n', stdDist);
        
        case 'frequency'
            gridSize = [10, 10]; % ??? ??
            [meanCells, stdCells] = frequencyAnalysis(image, gridSize);
            fprintf('Mean Cells per Grid: %.2f\n', meanCells);
            fprintf('STD of Cells per Grid: %.2f\n', stdCells);
        
        otherwise
            error('Unknown method: %s', method);
    end
end

