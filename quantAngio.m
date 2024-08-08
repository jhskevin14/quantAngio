%% Initiation
clear;
close all;

clc;

pix_size = 0.58;  % Ni2 20x  %0.49; % Olympus 10x %

code_dir = 'C:\Users\jhske\Dropbox\8_MatlabCode_HS\Angio_SID\';%/Users/jang_hwanseok/Dropbox/MatlabCode_HS/Angio_SID/';
im_dir = pwd;

addpath(im_dir);
addpath(code_dir);

if exist('Results','dir') == 0
    mkdir('Results');
    mkdir('Results/Figs');
end

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

% Call a sample
imStr = file_dir(fileInd).name;
im_temp = imread(imStr);
shortStr = char(strcat(extractBefore(imStr, '_4X'), extractBetween(imStr, 'RFP', 'f00')));

% Domain segmentation
clear btn_fn edt_1 hf sldr_1 txt_1

if fileInd == 1
filter_th = 0.0045; % Threshold for Matlab Sobel filter
end

dilate_rad = 20; % Radius of circular kernal used to dilate bright spots

t = 1;

[M, N] = size(im_temp);
[~, xedge, yedge, domain] = findBoundary(im_temp, filter_th, dilate_rad);

sldr_1.Value = filter_th;

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
    'min',0.001,'max',0.02,'Value',filter_th, 'SliderStep',[0.001 0.01],...
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
im_adj = im_crop .* uint8(dom_crop);
im_adj = imadjust(im_adj);
im_filt = imgaussfilt(im_adj,5);
im_filt2 = imadjust(im_filt);

sharpCoeff = [0 0 0;0 1 0;0 0 0]-fspecial('laplacian',0.2);
im_filt3 = imfilter(im_filt2, sharpCoeff, 'symmetric');

lowCutValue = 255;
highCutValue = 200;
toleranceValue = 100;

fft_img = fftshift(fft2(im_adj));
[rows, cols] = size(fft_img);
center_row = floor(rows/2) + 1;
center_col = floor(cols/2) + 1;

[X, Y] = meshgrid(1:cols, 1:rows);
distance = sqrt((X - center_col).^2 + (Y - center_row).^2);

low_cut_mask = 1 - exp(-((distance - lowCutValue).^2) / (2 * (toleranceValue^2)));
high_cut_mask = exp(-((distance - highCutValue).^2) / (2 * (toleranceValue^2)));
bandpass_mask = low_cut_mask .* high_cut_mask;

filtered_fft = fft_img .* bandpass_mask;

ifft_img = ifft2(ifftshift(filtered_fft));
im_fft = mat2gray(abs(ifft_img));
im_adfft = imadjust(im_fft);


% Brightness / Contrast
im_ft = im_adfft;
%    imshow(im_ft,[])

clear btn_fn edt_1 edt_2 hf sldr_1 txt_1 sldr_2 txt_2


if fileInd == 1
lowV = 0.28;
end

highV = 1;

sldr_1.Value = lowV;
sldr_2.Value = highV;
moWH = get(0, 'ScreenSize');

hf = figure;
set(hf, 'position', [0 0 moWH(3) moWH(4)]);
pos1 = [0 0.3 0.5 0.7];
fig1 = subplot('Position', pos1);
imshow(im_filt3, []); axis off;

pos2 = [0.1 0.18 0.3 0.1];
fig2 = subplot('Position', pos2);
imhist(im_adfft(im_adfft>0))
xlim([0 0.1])

pos3 = [0.5 0 0.5 1];
fig3 =  subplot('Position', pos3);

[im_th] = updateImage(im_ft, lowV, highV, pos3);

txt_1 = uicontrol(hf,'Style', 'text',...
    'String','Lower',...
    'Position', [200 100 100 20]);
txt_2 = uicontrol(hf,'Style', 'text',...
    'String','Upper',...
    'Position', [200 50 100 20]);

edt_1 = uicontrol(hf,'Style', 'edit',...
    'String',num2str(lowV),...
    'Position', [200 80 100 20]);
edt_2 = uicontrol(hf,'Style', 'edit',...
    'String',num2str(highV),...
    'Position', [200 30 100 20]);

sldr_1 = uicontrol(hf,'Style','slider',...
    'min',0,'max',1,'Value',lowV, 'SliderStep',[0.0001 0.001],...
    'callback','lowV = sldr_1.Value; edt_1.String=lowV; [im_th] = updateImage(im_ft, lowV, highV, pos3);',...
    'Position', [320 90 460 20]);

sldr_2 = uicontrol(hf,'Style','slider',...
    'min',0,'max',1,'Value',highV, 'SliderStep',[0.0001 0.001],...
    'callback','highV = sldr_2.Value; edt_2.String=highV; [im_th] = updateImage(im_ft, lowV, highV, pos3);',...
    'Position', [320 35 460 20]);

btn_fn = uicontrol(hf,'Style','pushbutton',...
    'String', 'Apply', ...
    'Position', [50 30 50 50],...
    'Callback', 'close;');
waitfor(hf);

imwrite(im_th,[im_dir '/Results/Figs/' shortStr '_Dom.tif'],'tif');


% Skeletonize

%    imshowpair(im_filt2, im_th, 'montage');
im_th = logical(im_th);
se = strel('disk',5);
im_dil = imdilate(im_th,se);

im_inv = ~im_dil;
im_cls = bwareaopen(im_inv, 100);
im_fin = ~im_cls;
im_fin = imerode(im_fin,se);

im_fin2 = bwareaopen(im_fin, 600);

im_skg = bwmorph(im_fin2,'thin',Inf);


imwrite(im_skg,[im_dir '/Results/Figs/' shortStr '_Skel.tif'],'tif');


% Quantification and Save

sSkg = regionprops(im_skg,'Area');

for k = 1 : length(sSkg)
    areas(k) = sSkg(k).Area;
end

maxSkg = max(areas);
numBrn = length(sSkg);

brPoints = bwmorph(im_skg,'branchpoints');
endPoints = bwmorph(im_skg,'endpoints');

numJun = length(brPoints);
numEnd = length(endPoints);
areaFin = length(find(im_fin2));
lengSkg = length(find(im_skg));
averTh = areaFin / lengSkg;

if exist([im_dir '\Results\Results.csv']) == 0
    fstrow = {'Sample', 'Total area', 'Total length', 'Average Thickness', ...
              'Branch Number', 'End points', ' Junction points', 'Longest Branch'};
    writecell(fstrow, [im_dir '\Results\Results.csv']);
end

saveData = {shortStr, areaFin, lengSkg, averTh, numBrn, numEnd, numJun, maxSkg};
writecell(saveData, [im_dir '\Results\Results.csv'], 'WriteMode','append');

save([im_dir '/Results/infoSave.mat'], 'fileInd', 'filter_th', 'lowV');

end
