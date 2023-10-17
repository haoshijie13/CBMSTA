 
%calculate the correlation of fmri gradient  and gene expression

% Add paths for required packages
addpath ./gifti-main

% Define species
species = 'macaque';
% Options:
% species = 'marmoset';
% species = 'mouse1';
% species = 'mouse2';

% Choose the top and last 1000 genes
max_cor_num = 1000;

% Load data

% Import fMRI gradient data
gradient_fmri = gifti('ceb_gradient1_to_tran.func.gii');
gradient_fmri = gradient_fmri.cdata;

% Import label data
label = load('label_5k.mat');
label = label.idx;

% Import all gene expression data
all_gene_expression = load('ceb_gene_expression_5k.mat');
gene_name = readmatrix('ceb_gene_name.csv', 'OutputType', 'string');

% Calculate mean gradient by label and compute correlation with gene expression

roinum = 5000;
mean_gradient_fmri = zeros(roinum, 1);

for r = 1:roinum   
    mean_gradient_fmri(r, :) = mean(gradient_fmri(find(label==r), :));
end

% Select appropriate gene expression data based on species
if strcmp(species, 'macaque')
    all_gene_expression = all_gene_expression.mean_all_surfdata;
else 
    all_gene_expression = all_gene_expression.mean_surfdata;
end

% Compute correlation between fMRI gradient and gene expression
fmri_cor_geneexpress = corr(mean_gradient_fmri, all_gene_expression);

% Sort the correlations and record the indices
[fmri_cor_genes, fmri_cor_genes_idx] = sort(fmri_cor_geneexpress); 

% Remove NaN values
NANvalue = isnan(fmri_cor_genes);
fmri_cor_genes(NANvalue) = [];
fmri_cor_genes_idx(NANvalue) = [];

% Choose the top and last highly correlated genes
max_fmri_cor_genes_idx = fmri_cor_genes_idx(end-max_cor_num+1:end)'; % Top 1000, from low to high
min_fmri_cor_genes_idx = fmri_cor_genes_idx(1:max_cor_num)'; % Last 1000, from high to low
max_fmri_cor_genes = fmri_cor_genes(end-max_cor_num+1:end)'; % Top 1000, from low to high
min_fmri_cor_genes = fmri_cor_genes(1:max_cor_num)'; % Last 1000, from high to low

% Select the gene names for the top and last correlated genes
max_fmri_cor_genes_name = gene_name(max_fmri_cor_genes_idx, 2);
min_fmri_cor_genes_name = gene_name(min_fmri_cor_genes_idx, 2);

% Save results to CSV files
final_max_name_file = [max_fmri_cor_genes_name, max_fmri_cor_genes, max_fmri_cor_genes_idx];
final_min_name_file = [min_fmri_cor_genes_name, min_fmri_cor_genes, min_fmri_cor_genes_idx];

writematrix(final_max_name_file, strcat('macaque_fmri_top', num2str(max_cor_num), '_genename.csv'));
writematrix(final_min_name_file, strcat('macaque_fmri_last', num2str(max_cor_num), '_genename.csv'));
i_cor_genes);

 fmri_cor_genes(NANvalue)=[];
 
 fmri_cor_genes_idx(NANvalue)=[];
 

 % choose the top and last related gene

max_fmri_cor_genes_idx = fmri_cor_genes_idx(1,end-max_cor_num+1:end)'; % Top list ,from low to high

min_fmri_cor_genes_idx = fmri_cor_genes_idx(1,1:max_cor_num)';%from high to low

max_fmri_cor_genes = fmri_cor_genes(1,end-max_cor_num+1:end)'; % Top list ,from low to high

min_fmri_cor_genes = fmri_cor_genes(1,1:max_cor_num)';%from high to low



 % choose the top  and last gene names

max_fmri_cor_genes_name = gene_name(max_fmri_cor_genes_idx,2);

min_fmri_cor_genes_name = gene_name(min_fmri_cor_genes_idx,2);


final_max_name_file =  [max_fmri_cor_genes_name,max_fmri_cor_genes,max_fmri_cor_genes_idx]

final_min_name_file =  [min_fmri_cor_genes_name,min_fmri_cor_genes,min_fmri_cor_genes_idx]


writematrix(final_max_name_file, strcat('macaque_fmri_top',num2str(max_cor_num),'_genename.csv'));
writematrix(final_min_name_file, strcat('macaque_fmri_last',num2str(max_cor_num),'_genename.csv'));




