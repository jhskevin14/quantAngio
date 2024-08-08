function [meanCells, stdCells] = frequencyAnalysis(image, gridSize)
    % image: ???? ?? ???
    % gridSize: ??? ?? (?: [10, 10])
    
    [rows, cols] = size(image);
    rowStep = floor(rows / gridSize(1));
    colStep = floor(cols / gridSize(2));
    
    cellCounts = zeros(gridSize);
    
    for i = 1:gridSize(1)
        for j = 1:gridSize(2)
            rowStart = (i-1) * rowStep + 1;
            colStart = (j-1) * colStep + 1;
            rowEnd = min(i * rowStep, rows);
            colEnd = min(j * colStep, cols);
            
            cellCounts(i, j) = sum(sum(image(rowStart:rowEnd, colStart:colEnd)));
        end
    end
    
    meanCells = mean(cellCounts(:));
    stdCells = std(cellCounts(:));
end
