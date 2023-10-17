% radar the gradient in transcriptomic space 
 % Add necessary paths
addpath ./gifti-main

% Define the number of ROIs
ROI_num = 15;

% Load data for macaque
macaque_gradient_gene = gifti('macaque_gene_gradient1_5k.func.gii');
macaque_gradient_fmri = gifti('macaque_fmri_gradient1_5k.func.gii');
macaque_gradient_fmri_ct = gifti('macaque_fmri_ct_gradient1_5k');
macaque_roi = gifti('macaque_ceb_ROI.func.gii');
macaque_roi = macaque_roi.cdata;

% Combine macaque and marmoset ROIs to match mouse
macaque_roi(find(macaque_roi == 2)) = 1;
macaque_roi(find(macaque_roi == 5)) = 4;

% Replace non-zero values with sorted values
macaque_roi = replaceNonZeroWithSorted(macaque_roi);

% Load data for marmoset
marmoset_gradient_gene = gifti('marmoset_gene_gradient1_5k.func.gii');
marmoset_gradient_fmri = gifti('marmoset_fmri_gradient1_5k.func.gii');
marmoset_gradient_fmri_ct = gifti('marmoset_fmri_ct_gradient1_5k.func.gii');
marmoset_roi = gifti('marmoset_ceb_ROI.func.gii');
marmoset_roi = marmoset_roi.cdata;

% Combine ROIs for matching mouse
marmoset_roi(find(marmoset_roi == 2)) = 1;
marmoset_roi(find(marmoset_roi == 5)) = 4;

% Replace non-zero values with sorted values
marmoset_roi = replaceNonZeroWithSorted(marmoset_roi);

% Load data for mouse1
mouse1_gradient_gene = gifti('mouse1_gene_gradient1_5k.func.gii');
mouse1_gradient_fmri = gifti('mouse1_fmri_gradient1_5k.func.gii');
mouse1_gradient_fmri_ct = gifti('mouse1_fmri_ct_gradient1_5k.func.gii');
mouse1_roi = gifti('mouse1_ROI.shape.gii');
mouse1_roi = mouse1_roi.cdata;

% Replace non-zero values with sorted values
mouse1_roi = replaceNonZeroWithSorted(mouse1_roi);

% Load data for mouse2
mouse2_gradient_gene = gifti('mouse2_gene_gradient1_5k.func.gii');
mouse2_gradient_fmri = gifti('mouse2_fmri_gradient1_5k.func.gii');
mouse2_gradient_fmri_ct = gifti('mouse2_fmri_ct_gradient1_5k.func.gii');
mouse2_roi = gifti('mouse2_ceb_ROI.func.gii');
mouse2_roi = mouse2_roi.cdata;

% Replace non-zero values with sorted values
mouse2_roi = replaceNonZeroWithSorted(mouse2_roi);

% Initialize arrays to store gradient data
macaque_gradient_3type = zeros(size(macaque_gradient_gene.cdata, 1), 3);
marmoset_gradient_3type = zeros(size(marmoset_gradient_gene.cdata, 1), 3);
mouse1_gradient_3type = zeros(size(mouse1_gradient_gene.cdata, 1), 3);
mouse2_gradient_3type = zeros(size(mouse2_gradient_gene.cdata, 1), 3);

% Populate gradient arrays
macaque_gradient_3type(:, 1) = macaque_gradient_gene.cdata;
macaque_gradient_3type(:, 2) = macaque_gradient_fmri.cdata;
macaque_gradient_3type(:, 3) = macaque_gradient_fmri_ct.cdata;

marmoset_gradient_3type(:, 1) = marmoset_gradient_gene.cdata;
marmoset_gradient_3type(:, 2) = marmoset_gradient_fmri.cdata;
marmoset_gradient_3type(:, 3) = marmoset_gradient_fmri_ct.cdata;

mouse1_gradient_3type(:, 1) = mouse1_gradient_gene.cdata;
mouse1_gradient_3type(:, 2) = mouse1_gradient_fmri.cdata(:, 1);
mouse1_gradient_3type(:, 3) = mouse1_gradient_fmri_ct.cdata(:, 1);

mouse2_gradient_3type(:, 1) = mouse2_gradient_gene.cdata;
mouse2_gradient_3type(:, 2) = mouse2_gradient_fmri.cdata(:, 1);
mouse2_gradient_3type(:, 3) = mouse2_gradient_fmri_ct.cdata(:, 1);


% Initialize arrays to store mean gradient data by ROI
macaque_gradient_3type_mean = zeros(ROI_num, 3);
marmoset_gradient_3type_mean = zeros(ROI_num, 3);
mouse1_gradient_3type_mean = zeros(ROI_num, 3);
mouse2_gradient_3type_mean = zeros(ROI_num, 3);

% Calculate mean gradient by ROI
for i = 1:ROI_num
    for j = 1:3
        idx = find(macaque_roi == i);
        macaque_gradient_3type_mean(i, j) = mean(macaque_gradient_3type(idx, j));
        
        idx = find(marmoset_roi == i);
        marmoset_gradient_3type_mean(i, j) = mean(marmoset_gradient_3type(idx, j));
        
        idx = find(mouse1_roi == i);
        mouse1_gradient_3type_mean(i, j) = mean(mouse1_gradient_3type(idx, j));
        
        idx = find(mouse2_roi == i);
        mouse2_gradient_3type_mean(i, j) = mean(mouse2_gradient_3type(idx, j));
    end
end

% Rescale the mean gradient data
rescale_min = 0.2;
rescale_max = 1;

for i = 1:3
    macaque_gradient_3type_mean(:, i) = rescale(macaque_gradient_3type_mean(:, i), rescale_min, rescale_max);
    marmoset_gradient_3type_mean(:, i) = rescale(marmoset_gradient_3type_mean(:, i), rescale_min, rescale_max);
    mouse1_gradient_3type_mean(:, i) = rescale(mouse1_gradient_3type_mean(:, i), rescale_min, rescale_max);
    mouse2_gradient_3type_mean(:, i) = rescale(mouse2_gradient_3type_mean(:, i), rescale_min, rescale_max);
end

% Define lobe labels for radar plots
%lobe_label = {'I-II', 'III', 'IV-V', 'VI', 'VII', 'VIII', 'IX', 'X', 'Fl', 'Cop', 'Par', 'Sim', 'Crus I', 'Crus II', 'PFl'};             
lobe_label = {'I-II', 'III', 'IV-V', 'VI', 'VII', 'VIII', 'IX', 'Cop', 'Par', 'Crus I', 'Crus II', 'Sim', 'X', 'Fl', 'PFl'};

% Define mapping for specific labels (ap)
mapping = [1, 2, 3, 4, 5, 6, 7, 10, 11, 13, 14, 12, 8, 9, 15];

% Rearrange the data for AP representation
all_gradient_ap = all_gradient(:, mapping);
all_gradient = all_gradient_ap;
all_gradient_num = size(all_gradient, 1);

% Create radar plots
for i = 1:3
    j = i + 3;
    y = i + 6;
    z = i + 9;
    figure
    ijyz = [i j y z];
    X = all_gradient(ijyz, :);
    RC = radarChart(X);
    RC.PropName = lobe_label;
    RC.ClassName = class_name(ijyz);
    RC = RC.draw();
    RC.legend();

    % Save radar plots as PDFs
    filename = strcat('radar_ap_scale_percent_', num2str(i), '.pdf');
    print('-dpdf', filename);
end



