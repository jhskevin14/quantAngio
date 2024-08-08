%% Initiation
clear;
close all;
clc;

code_dir = 'C:\Users\jhske\Dropbox\8_MatlabCode_HS\Angio_SID\quantAngio';%/Users/jang_hwanseok/Dropbox/MatlabCode_HS/Angio_SID/';
im_dir = pwd;

addpath(im_dir);
addpath(code_dir);

if exist('Results','dir') == 0
    mkdir('Results');
    mkdir('Results/Figs');
end
%
% Remove unnecessary folders
file_dir = dir([im_dir '/*.tif']);
i = 1;
for j = 1 : length(file_dir)
    if file_dir(j).name(1) == '.'
        omitdir(i) = j;
        i = i+1;
    end
end

if exist('omitdir') ~= 0
    file_dir(omitdir) = [];
end

if exist([im_dir '/Results/infoSave.mat']) == 0
    curInd = 1;
else
    load([im_dir '/Results/infoSave.mat']);
    curInd = fileInd;
end
%% Main
for fileInd = curInd : length(file_dir)

    imStr = file_dir(fileInd).name;
    im_temp = imread(imStr);
    shortStr = char(strcat(extractBefore(imStr, '_4X'), extractBetween(imStr, 'DAPI', 'f00')));

    % Domain segmentation

    clear btn_fn edt_1 hf sldr_1 txt_1
    
    filter_th = 0.02; % Threshold for Matlab Sobel filter
    dilate_rad = 20; % Radius of circular kernal used to dilate bright spots
    t = 1;

    [M, N] = size(im_temp);
    [~, xedge, yedge, domain] = findBoundary(im_temp, filter_th, dilate_rad);

    sldr_1.Value = filter_th;
    % sldr_2.Value = dilate_rad;

    hf = figure;
    set(hf, 'position', [500 10 1000 1000]);
    img_recall(1, im_temp, xedge, yedge);
    txt_1 = uicontrol(hf,'Style', 'text',...
            'String','Filter threshold',...
            'Position', [250 120 100 20]);
    edt_1 = uicontrol(hf,'Style', 'edit',...
            'String',num2str(filter_th),...
            'Position', [250 80 100 20]);
    sldr_1 = uicontrol(hf,'Style','slider',...
            'min',0.001,'max',1,'Value',filter_th, 'SliderStep',[0.0005 0.001],...
            'callback','filter_th = sldr_1.Value; edt_1.String=filter_th; [~, xedge, yedge, domain] = findBoundary(im_temp, filter_th, dilate_rad);img_recall(t, im_temp, xedge, yedge);',...
            'Position', [250 97 100 20]);
    btn_fn = uicontrol(hf,'Style','pushbutton',...
            'String', 'Apply', ...
            'Position', [800 80 50 20],...
            'Callback', 'close;');
    waitfor(hf);

    % Crop images

    im_temp = im_temp.*uint8(domain);
    s = regionprops(domain, 'Area', 'BoundingBox');

    [~,index] = max([s.Area]);

    bBox = s(index).BoundingBox;
    wth = 2500-1;
    imBox = [bBox(1)-round((wth-bBox(3))/2) bBox(2)-round((wth-bBox(4))/2) wth wth]; 

    im_crop = imcrop(im_temp, imBox);
    dom_crop = imcrop(domain, imBox);

    % Image filteration

    im_filt = imgaussfilt(im_crop,3);
    % imshow(im_filt)

    % Find Maxima

    clear btn_fn edt_1 hf sldr_1 txt_1
    

    maxima_th = 27;


    sldr_1.Value = maxima_th;

    hf = figure;
    set(hf, 'position', [500 10 1000 1000]);

    % imshow(im_crop, []);

    [centroids, maxima] = findMaxima(im_filt, sldr_1.Value);

    txt_1 = uicontrol(hf,'Style', 'text',...
            'String','Threshold',...
            'Position', [100 50 80 20]);
    edt_1 = uicontrol(hf,'Style', 'edit',...
            'String',num2str(maxima_th),...
            'Position', [200 50 100 20]);
    sldr_1 = uicontrol(hf,'Style','slider',...
            'min',1,'max',255,'Value',maxima_th, 'SliderStep',[0.01 0.1],...
            'callback','maxima_th = sldr_1.Value; edt_1.String=maxima_th; [centroids, maxima] = findMaxima(im_filt, maxima_th);',...
            'Position', [100 30 200 20]);
    btn_fn = uicontrol(hf,'Style','pushbutton',...
            'String', 'Apply', ...
            'Position', [700 40 50 20],...
            'Callback', 'close;');
    waitfor(hf);


    % Voronoi Segmentation

    [rows, cols] = size(im_crop);
    im_seg = zeros(rows, cols);

    for i = 1:rows
        for j = 1:cols
            distances = sqrt((centroids(:,1) - j).^2 + (centroids(:,2) - i).^2);
            [~, closestIdx] = min(distances);
            im_seg(i, j) = closestIdx;
        end
    end

    im_seg = im_seg.*dom_crop;
    im_seg(im_seg == 0) = nan;
    im_area = im_seg;

    for i = min(min(im_seg)) : max(max(im_seg))
        [row, col] = find(im_seg == i);
        for j = 1 :length(row)
            im_area(row(j), col(j)) = length(row);
        end
        areas(i) = length(row);
    end

    % Color map
    Cell_area = figure;
    set(Cell_area, 'color','w','position', [500 200 640 640]);
    set(gcf,'PaperPositionMode','auto');

    pcolor(im_area(:,:));
    set(gca,'Ydir','reverse');
    colormap(jet); shading interp; %flat;
    axis equal tight; axis off;
    c = colorbar('FontSize',15); c.Label.String = 'pixel'; c.Label.FontSize = 15;
    set(gca,'box','on','fontsize',11);
    set(gca,'color','k');
    set(gca,'CLim',[0 1500]);
    title('Voronoi Dimension','FontSize',15);

    print([im_dir '/Results/' shortStr '_VorDim'], '-dtiff');
    close;

        figure;
        violinplot(areas(areas ~= 0));
        title('Voronoi Dimension');
        xlabel(shortStr);
        ylabel('Pixel');

    print([im_dir '/Results/' shortStr '_VorDim_Violin'], '-dtiff');
    close;  

    % Quantifications
        
        meanVoro = mean(areas(areas ~= 0));
        stdVoro = std(areas(areas ~= 0));

        numPoints = size(centroids, 1);
        distances = pdist2(centroids, centroids);

    % RipleysK
        hDist = histogram(distances);
        maxDist = max(hDist.BinLimits);%3000;
        numBins = hDist.NumBins; close;
        ripleysK = zeros(numBins, 1);
        binEdges = linspace(0, maxDist, numBins+1);

        for i = 1:numBins
            binStart = binEdges(i);
            binEnd = binEdges(i+1);
            inBin = (distances > binStart) & (distances <= binEnd);
            ripleysK(i) = sum(inBin(:)) / numPoints;
        end

        ripleysK = cumsum(ripleysK);
        figure;
        plot(linspace(0, maxDist, numBins), ripleysK);
        title('Ripley''s K Function');
        xlabel('Distance');
        ylabel('K');

    print([im_dir '/Results/Figs/' shortStr '_RipleysK'], '-dtiff');
    close;  

    % Nearest Neighbor Distance
        distances(1:numPoints+1:end) = Inf;   
        nearestDist = min(distances, [], 2);   
        meanDist = mean(nearestDist);
        stdDist = std(nearestDist);
        figure;
        violinplot(nearestDist(nearestDist ~= 0));
        title('Nearest Neighbor Distance');
        xlabel(shortStr);
        ylabel('Pixel');

    print([im_dir '/Results/Figs/' shortStr '_NearND_Violin'], '-dtiff');
    close;

    % Cell Number Frequency in Grid
        gridSize = [30, 30]; % 30x30 square
        [rows, cols] = size(maxima);
        rowStep = floor(rows / gridSize(1));
        colStep = floor(cols / gridSize(2));

        cellCounts = zeros(gridSize);

        for i = 1:gridSize(1)
            for j = 1:gridSize(2)
                rowStart = (i-1) * rowStep + 1;
                colStart = (j-1) * colStep + 1;
                rowEnd = min(i * rowStep, rows);
                colEnd = min(j * colStep, cols);

                cellCounts(i, j) = sum(sum(maxima(rowStart:rowEnd, colStart:colEnd)));
            end
        end
        meanCells = mean(cellCounts(:));
        stdCells = std(cellCounts(:));
        figure;

        bar3(cellCounts)
        title('Cell Number Frequency');
        xlabel('Grid (au)');
        ylabel('Grid (au)');
        zlabel('Cell Number');

    print([im_dir '/Results/Figs/' shortStr '_CellFrequency'], '-dtiff');
    close;


    if exist([im_dir '\Results\Results.csv']) == 0
        fstrow = {'Sample', 'Voronoi Dimension (Average)', 'Voronoi Dimension (Std)', 'Nearest Neighbor Distance (Average)', 'Nearest Neighbor Distance (Std)', 'Cell number frequency (Average)', 'Cell number frequency (Std)'};
        writecell(fstrow, [im_dir '\Results\Results.csv']);
    end
    saveData = {shortStr, meanVoro, stdVoro, meanDist, stdDist, meanCells, stdCells};
    writecell(saveData, [im_dir '\Results\Results.csv'], 'WriteMode','append');


    save([im_dir '/Results/infoSave.mat'], 'fileInd');

end
