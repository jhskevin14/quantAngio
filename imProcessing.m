function imProcessing(rfp_img)
    % Create the figure and UI controls
    hFig = figure('Name', 'Image Processing GUI', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 600]);

    % Axes for displaying images
    hAx1 = axes('Parent', hFig, 'Position', [0.05, 0.55, 0.4, 0.4]);
    hAx2 = axes('Parent', hFig, 'Position', [0.55, 0.55, 0.4, 0.4]);
    hAx3 = axes('Parent', hFig, 'Position', [0.3, 0.05, 0.4, 0.4]);

    % Gamma slider
    uicontrol('Style', 'text', 'Parent', hFig, 'Position', [50, 500, 100, 20], 'String', 'Gamma');
    hGamma = uicontrol('Style', 'slider', 'Parent', hFig, 'Position', [150, 500, 200, 20], ...
        'Min', 0.1, 'Max', 5, 'Value', 1, 'Callback', @updateImage);

    % FFT Bandpass sliders
    uicontrol('Style', 'text', 'Parent', hFig, 'Position', [50, 450, 100, 20], 'String', 'Low Cutoff');
    hLowCut = uicontrol('Style', 'slider', 'Parent', hFig, 'Position', [150, 450, 200, 20], ...
        'Min', 1, 'Max', 100, 'Value', 1, 'Callback', @updateImage);
    uicontrol('Style', 'text', 'Parent', hFig, 'Position', [50, 400, 100, 20], 'String', 'High Cutoff');
    hHighCut = uicontrol('Style', 'slider', 'Parent', hFig, 'Position', [150, 400, 200, 20], ...
        'Min', 1, 'Max', 100, 'Value', 100, 'Callback', @updateImage);

    % Initial display
    imshow(rfp_img, 'Parent', hAx1);
    title(hAx1, 'Original Image');
    updateImage();

    function updateImage(~, ~)
        % Get slider values
        gammaValue = get(hGamma, 'Value');
        lowCutValue = get(hLowCut, 'Value');
        highCutValue = get(hHighCut, 'Value');
        
        % Process the image to create Domain.tif
        domain_img = process_image(rfp_img, gammaValue, lowCutValue, highCutValue);
        imshow(domain_img, 'Parent', hAx2);
        title(hAx2, 'Domain Image');
        
        % Skeletonize the Domain.tif image to create Skel.tif
        skel_img = bwmorph(domain_img, 'skel', Inf);
        imshow(skel_img, 'Parent', hAx3);
        title(hAx3, 'Skeletonized Image');
    end
end

function processed_img = process_image(img, gammaValue, lowCutValue, highCutValue)
    % Apply Gamma correction
    gamma_corrected = imadjust(img, [], [], gammaValue);
    
    % Apply FFT Bandpass filter
    fft_img = fftshift(fft2(gamma_corrected));
    [rows, cols] = size(fft_img);
    center_row = floor(rows/2) + 1;
    center_col = floor(cols/2) + 1;

    for i = 1:rows
        for j = 1:cols
            distance = sqrt((i - center_row)^2 + (j - center_col)^2);
            if distance < lowCutValue || distance > highCutValue
                fft_img(i, j) = 0;
            end
        end
    end

    ifft_img = ifft2(ifftshift(fft_img));
    processed_img = mat2gray(abs(ifft_img));
    
    % Threshold the image
    processed_img = imbinarize(processed_img);
end
