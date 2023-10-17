% combine the big gene expression data and calculate the correlation
 
 % Parameters
roinum = 5000;  % Number of ROIs
vertex = 535475;

% Load labels
labelf = load('label_5k.mat');
label = labelf.idx;

% Read data
surfdata = [];

for num = 0:14
    % Data files 1-14
    dataname = sprintf('sm15_all_gene_%02d000.func.gii', num);
    gg = gifti(dataname);
    surfdata1 = gg.cdata(:, 1:1000);
    surfdata = [surfdata, surfdata1];
    clear gg;
end

% Data file 15
dataname15 = 'sm15_all_gene_15000.func.gii';
gg15 = gifti(dataname15);
surfdata15 = gg15.cdata(:, :);
all_surfdata = [surfdata0, surfdata, surfdata15];

% Calculate mean of the surf data
timelength = size(all_surfdata, 2);
mean_all_surfdata = zeros(roinum, timelength);

for r = 1:roinum   
    mean_all_surfdata(r, :) = mean(all_surfdata(label == r, :));
end

% Save mean_all_surfdata
save('macaque_all_gene.mat', 'mean_all_surfdata', '-v7.3');

% Calculate correlations and mean them
cor_ceb = corrcoef(mean_all_surfdata');
save('macaque_ceb_allgene_cor.mat', 'cor_ceb', '-v7.3');

