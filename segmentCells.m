function segmentCells(imagePath, cellPath)
% segment cells in IF images under imagePath

if ~exist(cellPath,'dir')
    mkdir(cellPath);
end

readList = dir([imagePath '/*_green.jpg']); % read all protein channel
for i=1:length(readList)
    fprintf([num2str(i) '/' num2str(length(readList)) ' segmenting IF images...\n']);
    
    writePath = strrep([cellPath '/' readList(i,1).name], 'jpg', 'png'); % save cells in writePath
    nucPath = strrep([imagePath '/' readList(i,1).name], 'green', 'blue'); % nucleus channel
    erPath = strrep([imagePath '/' readList(i,1).name], 'green', 'yellow'); % ER channel
    
    if exist(writePath,'file')
        continue;
    end

    nucim = imread(nucPath); % nucleus
    cellim = imread(erPath); % cell
    
    IMAGEPIXELSIZE = 0.05; % um/px
    MINNUCLEUSDIAMETER = 6; % um
    MAXNUCLEUSDIAMETER = 14; % um
    
    regions = segmentation(nucim, cellim, MINNUCLEUSDIAMETER, MAXNUCLEUSDIAMETER, IMAGEPIXELSIZE);

    % saving results
    imwrite(regions, writePath);
end
fprintf('Done.\n');

function regions = segmentation( nucim, cellim, MINNUCLEUSDIAMETER, MAXNUCLEUSDIAMETER, IMAGEPIXELSIZE)

minarea = (MINNUCLEUSDIAMETER/IMAGEPIXELSIZE/2)^2*pi;
maxarea = (MAXNUCLEUSDIAMETER/IMAGEPIXELSIZE/2)^2*pi;

% convert nucleus and ER images from rgb to gray level 
nucim = rgb2gray(nucim); 
cellim = rgb2gray(cellim);

% determine foreground seeds
fgs = fgseeds(nucim, MINNUCLEUSDIAMETER, MAXNUCLEUSDIAMETER, IMAGEPIXELSIZE);

% process cell image
cellim_proc = getcellimg(cellim);

% determine background seeds
bgs = bgseeds( cellim_proc, MINNUCLEUSDIAMETER, IMAGEPIXELSIZE);

% combine seeds
seeds = fgs + bgs;

% perform seeded watershed
MA = max(cellim_proc(:));
regions = seededwatershed( MA-cellim_proc,seeds,4);