function [imFilt] = updateImage(image, lowV, highV, pos)
    
    [M, N] = size(image);
    imFilt = zeros([M N]);
    for i = 1 : M
        for j = 1 : N
            if image(i, j) >= lowV && image(i, j) <= highV
                imFilt(i, j) = 1;
            end
        end
    end
    subplot('Position', pos);
    imshow(imFilt, []); axis off;
end