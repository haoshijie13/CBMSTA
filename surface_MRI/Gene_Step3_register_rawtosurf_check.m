% Register the 2D slice of total-RNA gene expression to 2D slice from 3D volume nifti

% Define data paths
movfolder = 'total-RNA gene expression_path/'; % Path to moving (registered) images
movfiles = dir(fullfile(movfolder, '*.nii'));

fixfolder = '2D_slice_from_3D_volume_nifti_path/'; % Path to fixed (reference) images
fixfiles = dir(fullfile(fixfolder, '*.nii.gz'));

save_dir = 'check_result_path/'; % Path to save registration result images

% Loop through the image files
for i = 1:numel(movfiles)
    curmovID = movfiles(i).name;
    curID = curmovID(1:4);
    
    disp(curID)
    
    curfixID = fixfiles(i).name;

    % Read the image files
    movfile = fullfile(movfolder, curmovID); % Registered totalRNA file in 3D
    fixfile = fullfile(fixfolder, curfixID); % Raw totalRNA file in 2D
    
    fiximg = niftiread(fixfile);
    fixsize = size(fiximg);
    fiximg = squeeze(fiximg);

    movimg = squeeze(niftiread(movfile));

    % Rotate if needed (Example: flipping horizontally)
    movimg = flipud(movimg); % x-axis flip

    % Remove high values from the data if needed
    % Example: Limiting the intensity values to 2000
    % if i == 2
    %     fiximg(fiximg >= 2000) = 2000;
    % end

    % Register the images
    movresult = registerImages(movimg, fiximg);

    % Register the images again if needed
   
    % movresult2 = registerImages(movresult.RegisteredImage, fiximg);

    % Check the registration result by showing the images
    obj = imshowpair(fiximg, fiximg);
    saveas(obj, [save_dir, curID, 'fix.png'])
    
    obj = imshowpair(fiximg, movresult.RegisteredImage);
    saveas(obj, [save_dir, curID, 'fixresult.png'])
    
    obj = imshowpair(fiximg, movresult.RegisteredImage, 'montage');
    saveas(obj, [save_dir, curID, 'both.png'])
end

% Check the registration result, 'fix.png' and 'fixresult.png' should be in
% the same space, and 'both.png' shows complete overlap between two images.

