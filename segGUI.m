function segGUI(image)

    image = mat2gray(image); % ??? ???

    % ?? ??? ??
    initialThreshold = 0.5;

    % Figure ??
    hFig = figure('Name', 'DAPI Segmentation with Threshold Adjustment', 'NumberTitle', 'off');
    
    % ??? ??
    hAx = axes('Parent', hFig);
    imshow(image, [], 'Parent', hAx);
    hold on
%    hold(hAx, 'on');
    
    % ???? ??
    uicontrol('Style', 'text', 'String', 'Threshold', 'Position', [20 20 60 20]);
    hSlider = uicontrol('Style', 'slider', 'Min', 0, 'Max', 1, 'Value', initialThreshold, ...
                        'Position', [100 20 300 20], 'Callback', @(src, event) updateSegmentation(src, image, hAx));
    
    % ?? ?????? ??
    updateSegmentation(hSlider, image, hAx);
end

function updateSegmentation(slider, image, hAx)
    % ???? ? ????
    threshold = get(slider, 'Value');
    
    % ???? ??? ??
    binaryImage = imbinarize(image, threshold);

    imshow(image, [], 'Parent', hAx);
    hold on
    % ?? ??? ?? ???? ??? ??
    boundaries = bwboundaries(binaryImage);
    for k = 1:length(boundaries)
        boundary = boundaries{k};
        plot(hAx, boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1);
    end
    hold off
    title(hAx, sprintf('Threshold: %.2f', threshold));
end
