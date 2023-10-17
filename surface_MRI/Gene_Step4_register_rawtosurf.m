% ==register the 2D slice of all gene expression to 2D slice from 3D volume nifti

movfolder = 'total-RNA gene expression_path/';   
movfiles = dir(fullfile(movfolder, '*.nii'));

geneexpfolder = 'all gene expression_nifti_raw_path';
geneexpfiles = dir(fullfile(geneexpfolder, '*.nii'));

fixfolder = '2D_slice_from_3D_volume_nifti_path/';
fixfiles = dir(fullfile(fixfolder, '*.nii.gz'));

savepath = 'result_path/';

for i = 1: numel(movfiles)

    curmovID = movfiles(i).name;
    
    curID = curmovID(1:4);
    
    disp(curID)
    
    curfixID = fixfiles(i).name;

    curgeneID = strcat(curID,'.nii');

    % read the file
    
    movfile = fullfile(movfolder, curmovID);   % reg totalRNA file in 3D
    
    fixfile = fullfile(fixfolder, curfixID);  % raw totalRNA file

    geneexpfile = fullfile(geneexpfolder, curgeneID);  
      

    fiximg = niftiread(fixfile);
    fixsize = size(fiximg);
    fiximg = squeeze(fiximg);

    movimg = squeeze(niftiread(movfile));

    mov_first5000 = squeeze(niftiread(geneexpfile));

    % Rotate if needed (Example: flipping horizontally)
 
    movimg = flipud(movimg); % x-axis flip
    mov_first5000 = flipud(mov_first5000);% x-axis flip


    % Remove high values from the data if needed
    % Example: Limiting the intensity values to 2000
    % if i == 2
    %     fiximg(fiximg >= 2000) = 2000;
    % end


    % Register the images
    
    movresult = registerImages(movimg, fiximg);

    % Register the images again if needed
   
    % movresult2 = registerImages(movresult.RegisteredImage, fiximg);

    % check the result by show the img
    obj = imshowpair(fiximg, fiximg);
      
    %  figure
    %  obj = imshowpair(movimg, movimg);

      
    % ====apply the Transformation of the registration
        
    regresult_first5000 = imwarp(mov_first5000, movresult.Transformation, 'OutputView', imref2d(size(fiximg)));
           
    % if the dim is not true
    % regresult_first5000 = permute(regresult_first5000, [1, 3, 2]);

    info = niftiinfo(fixfile);
    cursavepath = [savepath, filesep, 'reg_gene_cluster/', curID];

    if exist(cursavepath, 'dir') == 0
        mkdir(cursavepath)
    end

    for j = 1: size(regresult_first5000, 3)
            
        % niftiwrite(single(regresult_first5000(:, j, :)), [cursavepath, '/', curID, '-', num2str(j, '%05d')], info, 'Compressed', true)
        niftiwrite(single(regresult_first5000(:, :, j)), [cursavepath, '/', curID, '-', num2str(j, '%05d')],'Compressed', true)

    end
end

