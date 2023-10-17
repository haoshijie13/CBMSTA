% find the same gene between species

% Read gene lists from CSV files for macaque and marmoset
list1 = readmatrix('macaque_fmri_corrected_gene_top1000.csv', 'OutputType', 'string');
list1 = list1(:, 1);

list2 = readmatrix('marmoset_fmri_corrected_gene_top1000.csv', 'OutputType', 'string');
list2 = list2(:, 1);

% Initialize cell arrays to store matching gene names
matches1 = {};
matches2 = {};

% Loop through both lists to find matching gene names
for i = 1:numel(list1)
    str1 = list1{i};
    
    for j = 1:numel(list2)
        str2 = list2{j};
        
        % Compare two strings, case insensitively
        if strcmpi(str1, str2)
            matches1{end+1} = str1;
            matches2{end+1} = str2;
        end
    end
end

% Display matching gene names
if ~isempty(matches1)
    disp(matches1);
    disp(matches2);
else
    disp('Nan');
end

% Save overlapping gene names to a CSV file
list1 = readmatrix('gene_list1.csv', 'OutputType', 'string');
macaque_marmoset_ls_1 = list1(ismember(list1, matches1), :);
writematrix(macaque_marmoset_ls_1, 'overlap_gene_name.csv');

