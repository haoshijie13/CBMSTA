 % diffusion embedding  for the transcriptomic data 

 % Add necessary paths
addpath('./gifti-main');
addpath('/BrainSpace-0.1.2');

% Load data
gg2 = gifti('macaque_totalRNA_scale100_ceb.shape.func.gii');

roinum = 5000;  % Number of ROIs
vertex = 535475;

labelf = load('label_5k.mat');   
label = labelf.idx;

gene_f = load('macaque_gene_exprassion_5k.mat');
gene = gene_f.mean_all_surfdata;

cor_ceb = corrcoef(gene');

% Perform diffusion embedding
gm = GradientMaps();
gm = gm.fit(cor_ceb);

% Plot the components
scree_plot(gm.lambda{1});
lambda = gm.lambda{1};
save('macaque_gene_lambda_5k.mat', 'lambda');

% Save gradients
for gra = 1:8
    label2 = label;

    for rr = 1:roinum
        label2(label2 == rr) = gm.gradients{1}(rr, gra);
    end

    ROI_name = sprintf('macaque_gene_gradient%d_5k.func.gii', gra);
    gg2.cdata = label2;
    saveg(gg2, ROI_name);
end

