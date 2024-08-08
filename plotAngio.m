%% Initiation

clear;
close all;
clc;

% Read data
[results, rsDir] = uigetfile('*.csv');
rdData = readcell([rsDir results]);

% Remove duplicate samples
[~,ia,~] = unique(rdData(:,1));
ia = sort(ia);
rdData = rdData(ia,:);

% Extract group names
groups = extractBefore([rdData(2:end,1)], '_'); %4); %
grCat = unique(groups);

% Plot and save each metric
grayCol = colormap(gray(length(grCat)+4)); close;

for i = 2 : size(rdData,2)
    for j = 1 : length(grCat)
        IndexC = strfind(groups,cell2mat(grCat(j)));
        grInd = find(not(cellfun('isempty',IndexC)));
        valTemp = cell2mat(rdData(grInd+1,i));
        bar(j,mean(valTemp),'FaceColor', grayCol(j+2,:),'LineWidth', 1.5); hold on
        errorbar(j,mean(valTemp),std(valTemp), 'Color', 'k', 'LineWidth', 1.5, 'CapSize', 15);
    end
        xticks(1:length(grCat));
        xticklabels(grCat);
        ylabel(cell2mat(rdData(1,i)));
        set(gca,'FontSize',15);

        print([rsDir cell2mat(rdData(1,i))], '-dtiff');
        close;  
end

[mM, fF] = mode(cell2mat(groups));

for i = 2 : size(rdData,2)
    allTemp = NaN([length(grCat), fF]);
    for j = 1 : length(grCat)
        IndexC = strfind(groups,cell2mat(grCat(j)));
        grInd = find(not(cellfun('isempty',IndexC)));
        valTemp = cell2mat(rdData(grInd+1,i));
        allTemp(j,1:length(valTemp)) = valTemp';
    end
        violinplot(allTemp', grCat, 'ViolinColor', grayCol);
        
        xticks(1:length(grCat));
        xticklabels(grCat);
        ylabel(cell2mat(rdData(1,i)));
        set(gca,'FontSize',15);

        print([rsDir cell2mat(rdData(1,i)) '_vin'], '-dtiff');
        close;  
end

