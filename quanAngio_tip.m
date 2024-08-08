%% Initiation

clear;
close all;
clc;

% pix_size = 0.87;  % Ni2 20x  %0.49; % Olympus 10x %

code_dir = 'C:\Users\jhske\Dropbox\8_MatlabCode_HS\Angio_SID\';%/Users/jang_hwanseok/Dropbox/MatlabCode_HS/Angio_SID/';
im_dir = pwd;

addpath(im_dir);
addpath(code_dir);

if exist('Results','dir') == 0
    mkdir('Results');
    mkdir('Results/Figs');
end

% Remove unnecessary folders
file_dir = dir([im_dir '/*.jpg']);
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

%% Call a sample 
for fileInd = curInd : length(file_dir)

    fileInd = 10;
imStr = file_dir(fileInd).name;
im_temp = imread(imStr);
if size(im_temp,3) == 3
    im_temp = rgb2gray(im_temp);
end
shortStr = char(strcat(extractBetween(imStr, 'rfp_', '.jpg')));
%% Channel barrier detection

im_filt = imadjust(im_temp);
clear btn_fn edt_1 edt_2 hf sldr_1 txt_1 sldr_2 txt_2


if fileInd == 1
    filter_th = 155;
end

highV = 255;

sldr_1.Value = filter_th;
sldr_2.Value = highV;
moWH = get(0, 'ScreenSize');

hf = figure;
set(hf, 'position', [0 0 moWH(3) moWH(4)]);
pos1 = [0 0.3 0.5 0.7];
fig1 = subplot('Position', pos1);
imshow(im_filt, []); axis off;

pos2 = [0.1 0.18 0.3 0.1];
fig2 = subplot('Position', pos2);
imhist(im_filt(im_filt>0))
xlim([0 255])

pos3 = [0.5 0 0.5 1];
fig3 =  subplot('Position', pos3);

[im_th] = updateImage(im_filt, filter_th, highV, pos3);

txt_1 = uicontrol(hf,'Style', 'text',...
    'String','Lower',...
    'Position', [200 100 100 20]);
txt_2 = uicontrol(hf,'Style', 'text',...
    'String','Upper',...
    'Position', [200 50 100 20]);

edt_1 = uicontrol(hf,'Style', 'edit',...
    'String',num2str(filter_th),...
    'Position', [200 80 100 20]);
edt_2 = uicontrol(hf,'Style', 'edit',...
    'String',num2str(highV),...
    'Position', [200 30 100 20]);

sldr_1 = uicontrol(hf,'Style','slider',...
    'min',0,'max',255,'Value',filter_th, 'SliderStep',[0.01 0.1],...
    'callback','filter_th = sldr_1.Value; edt_1.String=filter_th; [im_th] = updateImage(im_filt, filter_th, highV, pos3);',...
    'Position', [320 90 460 20]);

sldr_2 = uicontrol(hf,'Style','slider',...
    'min',0,'max',255,'Value',highV, 'SliderStep',[0.01 0.1],...
    'callback','highV = sldr_2.Value; edt_2.String=highV; [im_th] = updateImage(im_filt, filter_th, highV, pos3);',...
    'Position', [320 35 460 20]);

btn_fn = uicontrol(hf,'Style','pushbutton',...
    'String', 'Apply', ...
    'Position', [50 30 50 50],...
    'Callback', 'close;');
waitfor(hf);
%%

[H,T,R] = hough(im_th);
P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(im_th,T,R,P,'FillGap',5,'MinLength',20);
pnts = 1;

for k = 1:length(lines)
    if abs(lines(k).theta) == 90
        sPnt(pnts,:) = lines(k).point1;
        ePnt(pnts,:) = lines(k).point2;
        pnts = pnts + 1;
    end
end

cutLine = mode([sPnt(:,2);ePnt(:,2)]);

%% Confirm the line
clear roi manLine btn_fn1 hf txt_1 btn_fn2;

close;

hf = figure;
set(hf, 'position', [500 10 1000 1000]);

imshow(im_temp,[]), hold on
plot([1 size(im_temp,2)], [cutLine cutLine],'--','LineWidth',2,'Color','green');

txt_1 = uicontrol(hf,'Style', 'text',...
    'String','Is the green line on the channel barrier?',...
    'Position', [250 120 250 20]);
btn_fn1 = uicontrol(hf,'Style','pushbutton',...
    'String', 'Yes', ...
    'Position', [250 80 50 20],...
    'Callback', 'manLine = [1 cutLine; size(im_temp,2) cutLine]; close;');

btn_fn2 = uicontrol(hf,'Style','pushbutton',...
    'String', 'No, I would draw a line manually.', ...
    'Position', [310 80 190 20],...
    'Callback', 'manLine = drawline; roi = manLine.Position;');

waitfor(hf);

if exist("roi") ~= 0
    manLine(:,2) = round(mean(roi(:,2))); 
end

im_cut = im_temp;
im_cut(1:manLine(1,2),:) = 0;


%% Image segmentation

im_adj = imadjust(im_cut);
im_filt = histeq(im_adj);

im_filt(im_filt < 160) = 160;
im_filt = imadjust(im_filt);


clear btn_fn edt_1 edt_2 hf sldr_1 txt_1 sldr_2 txt_2


if fileInd == 1
    lowV = 214;
end

highV = 255;

sldr_1.Value = lowV;
sldr_2.Value = highV;
moWH = get(0, 'ScreenSize');

hf = figure;
set(hf, 'position', [0 0 moWH(3) moWH(4)]);
pos1 = [0 0.3 0.5 0.7];
fig1 = subplot('Position', pos1);
imshow(im_filt, []); axis off;

pos2 = [0.1 0.18 0.3 0.1];
fig2 = subplot('Position', pos2);
imhist(im_filt(im_filt>0))
xlim([0 255])

pos3 = [0.5 0 0.5 1];
fig3 =  subplot('Position', pos3);

[im_th] = updateImage(im_filt, lowV, highV, pos3);

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
    'min',0,'max',255,'Value',lowV, 'SliderStep',[0.01 0.1],...
    'callback','lowV = sldr_1.Value; edt_1.String=lowV; [im_th] = updateImage(im_filt, lowV, highV, pos3);',...
    'Position', [320 90 460 20]);

sldr_2 = uicontrol(hf,'Style','slider',...
    'min',0,'max',255,'Value',highV, 'SliderStep',[0.01 0.1],...
    'callback','highV = sldr_2.Value; edt_2.String=highV; [im_th] = updateImage(im_filt, lowV, highV, pos3);',...
    'Position', [320 35 460 20]);

btn_fn = uicontrol(hf,'Style','pushbutton',...
    'String', 'Apply', ...
    'Position', [50 30 50 50],...
    'Callback', 'close;');
waitfor(hf);

im_cls = bwareaopen(im_th, 50);

se = strel('disk',6);
se2 = strel('disk',6);
im_dil = imdilate(im_cls,se);
im_fin = imerode(im_dil,se2);

domSave = im_fin;
domSave(1:manLine(1,2),:) = 1;

%% Domain image confirmation

clear recPos manBox btn_fn1 hf txt_1 btn_fn2;

close;

hf = figure;
set(hf, 'position', [0 0 moWH(3) moWH(4)]);

pos1 = [0.025 0.1 0.45 1];
fig1 = subplot('Position', pos1);
imshow(im_temp,[]), hold on

pos2 = [0.525 0.1 0.45 1];
fig2 = subplot('Position', pos2);
imshow(domSave,[]), hold on


txt_1 = uicontrol(hf,'Style', 'text',...
    'String','Dose the domain image look fine?',...
    'Position', [250 120 250 20]);
btn_fn1 = uicontrol(hf,'Style','pushbutton',...
    'String', 'Yes/Close', ...
    'Position', [225 80 75 20],...
    'Callback', 'close;');

btn_fn2 = uicontrol(hf,'Style','pushbutton',...
    'String', 'No, I want to remove particles.', ...
    'Position', [310 80 190 20],...
    'Callback', 'manBox = drawrectangle; recPos = round(manBox.Position);');

waitfor(hf);

if exist("recPos") ~= 0
    domSave(recPos(2):recPos(2)+recPos(4),recPos(1):recPos(1)+recPos(3)) = 0; 
    im_fin(recPos(2):recPos(2)+recPos(4),recPos(1):recPos(1)+recPos(3)) = 0; 
end
imwrite(domSave,[im_dir '/Results/Figs/' shortStr '_Dom.tif'],'tif');
%%



im_skg = bwmorph(im_fin,'thin',Inf);

im_skg_bk = im_skg;
im_fin_bk = im_fin;

[M, N] = size(im_temp);

for j = 1 : N
    [temRow,~] = find(im_skg(:,j));
    im_skg(1:2*min(temRow)- manLine(1,2),j) = 0;
    im_fin(1:2*min(temRow)- manLine(1,2),j) = 0;
end


% skeleton > 10, branch > 200
im_skopen = bwareaopen(im_skg, 10);
im_fnopen = bwareaopen(im_fin, 200);


imwrite(im_skopen,[im_dir '/Results/Figs/' shortStr '_Skel.tif'],'tif');

%% Quantification and Save
sSkg = regionprops(im_skopen,'Area');
sFn = regionprops(im_fnopen,'Area');

brPoints = bwmorph(im_skopen,'branchpoints');
endPoints = bwmorph(im_skopen,'endpoints');

[epRow,~] = find(endPoints);
maxProt = max(epRow) - min(epRow);

for k = 1 : length(sSkg)
    areas(k) = sSkg(k).Area;
end

maxSkg = max(areas);
minSkg = min(areas);
stdSkg = std(areas);

numJun = length(find(brPoints));
numEnd = length(find(endPoints));
areaFin = length(find(im_fnopen));
lengSkg = length(find(im_skopen));
aveThk = areaFin / lengSkg;
numBrn = length(sFn);

if exist([im_dir '\Results\Results.csv']) == 0
    fstrow = {'Sample', 'Total area', 'Total length', 'Average Thickness',...
              'Branch Number', 'End points', ' Junction points',...
              'Longest Branch', 'Shortest Branch', 'Std length'};
    writecell(fstrow, [im_dir '\Results\Results.csv']);
end

saveData = {shortStr, areaFin, lengSkg, aveThk, numBrn, numEnd, numJun, maxSkg, minSkg, stdSkg};
writecell(saveData, [im_dir '\Results\Results.csv'], 'WriteMode','append');

save([im_dir '/Results/infoSave.mat'], 'fileInd', 'filter_th', 'lowV');

end