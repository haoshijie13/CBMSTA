% Diffusion embedding for fMRI data at the subject level

% Add necessary paths
addpath('./gifti-main');
addpath('/BrainSpace-0.1.2');

% Set the working directory
cd('data_dir/');

% Load the template surface and fMRI data
gg2 = gifti('Tmeansub-IONm01_ses-01_LR.func.gii');
roinum = 5000;  % Number of ROIs

% Load label information
labelf = load('label_5k.mat');
label = labelf.idx;
vertex = size(label, 1);
sub_num = 39;

% Process data for each subject
for sub = 1:sub_num
    s1 = 'sub-*m';
    s2 = num2str(sub);
    s3 = '*';
    s4 = 'func.gii';
    session_name = strcat(s1, s2, s3, s4);
    session_ls = dir(session_name);  % Find the smooth* files
    session_num = length(session_ls);
    
    cor_ceb_allses = zeros(roinum, roinum);
    
    for ses = 1:session_num
        dataname = session_ls(ses).name;
        gg1 = gifti(dataname);
        surfdata = gg1.cdata;

        % Compute mean surface data for each ROI
        timelength = size(surfdata, 2);
        mean_surfdata = zeros(roinum, timelength);
        
        for r = 1:roinum
            mean_surfdata(r, :) = mean(surfdata(find(label == r), :));
        end

        % Compute correlations and accumulate them
        mean_surfdata = mean_surfdata';
        cor_ceb = corrcoef(mean_surfdata);
        cor_ceb_allses = cor_ceb_allses + cor_ceb;
    end

    % Calculate the mean correlation for the subject
    cor_ceb_mean = cor_ceb_allses ./ session_num;

    % Save the mean correlation
    s5 = 'sub-m';
    s6 = num2str(sub);
    s7 = '_mean_cor_5k.mat';
    savename = strcat(s5, s6, s7);
    save(savename, 'cor_ceb_mean', '-v7.3');

    % Perform diffusion embedding
    gm = GradientMaps();
    gm = gm.fit(cor_ceb_mean);

    % Save gradient maps
    for gra = 1:8
        label2 = label;

        for rr = 1:roinum
            label2(label2 == rr) = gm.gradients{1}(rr, gra);
        end

        s5 = 'sub-m';
        s6 = num2str(sub);
        s7 = '_gradient';
        s8 = num2str(gra);
        s9 = '_5k.func.gii';
        save_name = strcat(s5, s6, s7, s8, s9);
        gg2.cdata = label2;
        saveg(gg2, save_name);
    end
end

% Calculate the group correlation and perform group-level diffusion embedding
cor_allsub = [];

for sub = 1:sub_num
    s1 = 'sub-m';
    s2 = num2str(sub);
    s3 = '_mean_cor_5k.mat';
    sub_name = strcat(s1, s2, s3);
    cor = load(sub_name);
    cor_allsub = cor_allsub + cor.cor_ceb_mean;
end

cor_allsub_mean = cor_allsub ./ sub_num;

% Save the group-level correlation
save('sub-m1-39_mean_cor_5k.mat', 'cor_allsub_mean', '-v7.3');

% Perform diffusion embedding at the group level
gm = GradientMaps();
gm = gm.fit(cor_allsub_mean);

% Plot the components (you can uncomment these lines)
% scree_plot(gm.lambda{1});
% gradient_in_euclidean(gm.gradients{1}(:, 1:3));

% Save gradient maps for group-level
for gra = 1:8
    label2 = label;

    for rr = 1:roinum
        label2(label2 == rr) = gm.gradients{1}(rr, gra);
    end

    s5 = 'sub_m1-39_mean_5k_cor';
    s6 = '_gradient';
    s7 = num2str(gra);
    s8 = '.func.gii';
    save_name = strcat(s5, s6, s7, s8);
    gg2.cdata = label2;
    saveg(gg2, save_name);
end

