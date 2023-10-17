 
% input the gene name and generate the gene expression gifti file 
% Add the path to the gifti-main toolbox
addpath ./gifti-main

% Define the species (choose one)
species = 'macaque';
% species = 'marmoset';
% species = 'mouse1';
% species = 'mouse2';

% Load data for macaque
if strcmp(species, 'macaque')
    % Load label data
    label = load('label_5k.mat');
    label = label.idx;
    
    % Load gene expression data
    all_gene_expression = load('ceb_gene_expression_5k.mat');
    all_gene_expression = all_gene_expression.mean_all_surfdata;
    
    % Load example gene names
    gene_name = readmatrix('ceb_gene_name.csv', 'OutputType', 'string');
    
    % Specify gene names to choose
    gene_name_choose = ["NTF3", "CNBD1", "MEIS2", "KCNV1"];
    
    % Extract the chosen gene expressions
    all_gene_expression_choose = all_gene_expression(:, ismember(gene_name(:, 2), gene_name_choose));
    
    % Extract the names of the chosen genes
    all_gene_expression_choose_name = gene_name(ismember(gene_name(:, 2), gene_name_choose), 2);
    
    % Load the gradient data
    gg2 = gifti('gradient1_5k.func.gii');
else
    disp('Species is not macaque');
end

% Save gene expression data as gifti files for each selected gene
roinum = 5000;

for gene_num = 1:size(all_gene_expression_choose, 2)
    % Copy the label data
    label2 = label;
    
    % Replace labels with gene expression values
    for rr = 1:roinum
        label2(label2 == rr) = all_gene_expression_choose(rr, gene_num);
    end
    
    % Define the filename for saving
    save_name = strcat(species, '_fmri_corrected_gene_', char(all_gene_expression_choose_name(gene_num)), '.func.gii');
    
    % Set the gifti data
    gg2.cdata = label2;
    
    % Save the gifti file
    saveg(gg2, save_name);
end

