%% These scripts accompany the manuscript:
%  Gulbinaite et al. (2024) "Spatiotemporal resonance in mouse visual
%  cortex" Curr Biol
%
%   MATLAB version used 2019b
%% Code for recording ALIGNEMENT across days and Allen CCF

clear all
close all
homedir = '...\Data\1_Glu\'; % add path to the data
cd(homedir)
filelist = dir('*.tif');
filelist = {filelist.name};

Bregref = [10 30]; % reference Bregma point on Allen CCF

% Find the first recording 
% Note: "0_Focused_HD" is not the first recording
% but HD anatomical photos
fileidx = regexp(filelist,'1_i*.*Focused_video');
fileidx = find(~cellfun(@isempty,fileidx));

img = imread([homedir filelist{fileidx}],'tiff');
img = double(img);
img = (2^16-1)*img/max(max(img));
img = uint16(img);
img = imresize(img, 0.125);

%% adjust image contrast to see bregma/lambda better
h1 = figure; hold on
imshow(img);
imcontrast(h1)
%% Select bregma and lambda points on the brain image
title('Click on Bregma on your image, then press enter','fontsize',12);
[x,y] = getpts(h1);
Bregpoint=round([x y]);
hold on
plot(Bregpoint(:,1),Bregpoint(:,2),'xr','linewidth',2,'markersize',5);

title('Click on Lambda on your image, then press enter','fontsize',12);
[x,y] = getpts(h1);
Lambpoint=round([x y]);
plot(Lambpoint(1),Lambpoint(2),'xm','linewidth',2,'markersize',5);

% Calculate Translation relative to Bregref 
% (this way bregma is at the same spot across animals and recordings)
dx = Bregref(1)-Bregpoint(1);
dy = Bregref(2)-Bregpoint(2);

% Calculate Rotation (such that bregma and lambda is on the straight line)
Lambpoint2 = [0 Lambpoint(2)];     % another point on the lambda line
pt = [Bregpoint 0];
v1 = [Lambpoint 0];
v2 = [Lambpoint2 0];
distance = norm(cross(v1-v2,pt-v2))/norm(v1-v2);

angle2rotate = acosd(distance/sqrt( (Bregpoint(1)-Lambpoint(1))^2 + (Bregpoint(2)-Lambpoint(2))^2) );

%% Translate and rotate the image

% rotate
rotated = rotateAround(img, Bregpoint(2), Bregpoint(1),  -angle2rotate,'bicubic');
figure
imshow(rotated)

% translate
img_transf = imtranslate(rotated,[dx, dy]);

% plot the final result
h3 = figure;
imshow(img_transf)
hold on
plot(Bregref(1),Bregref(2),'xm','linewidth',2,'markersize',5);
% imcontrast(h3)

% Save the image and all the reference points
imwrite(img_transf, [homedir filelist{fileidx}(1:end-17) 'rotated_translated2.tif'])
save([homedir filelist{fileidx}(1:end-17) 'rotated_translated2.mat'],'img_transf','angle2rotate','Bregref','Bregpoint','Lambpoint')

%% Co-register
clear all 
close all

% Files and directories
homedir = '...\Data\1_Glu\';
writedir = '...\Data\1_Glu\';

cd(homedir)
datafiles   = dir('*video.tif');
datafiles = {datafiles.name};

reference_img_file = dir('*translated2.tif');
img = imread([homedir reference_img_file.name],'tiff');
disp(reference_img_file.name)

% Plot the image
% h1 = figure; hold on
% imshow(img);
% imcontrast(h1)

%% LOOP across recordings 

for filei = 2:length(datafiles)

    
    %% importing TO-BE-TRANSFORMED image / unregistered
    disp(datafiles{filei})
    img_unreg = imread([homedir datafiles{filei}],'tiff');
    
    img_unreg = double(img_unreg);
    img_unreg = (2^16-1)*img_unreg/max(max(img_unreg));
    img_unreg = uint16(img_unreg);
    img_unreg = imresize(img_unreg, 0.125);
    
    %% Selecting alignment points (movingPoints,fixedPoints) and exporting it to the workspace 
    % TIP: alignement between recordings is easiest done by selecting 
    % blood vessel bifurcation points
    cpselect(img_unreg,img)
    
%     cpselect(img_unreg,img,movingPoints,fixedPoints)

    %% Finding the transformation
    mytform = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity'); % 
    
    % Applying the transformation to the brain image
    img_reg = imwarp(img_unreg,mytform,'outputview',imref2d(size(img)));
    
    %% Plot to verify the quality of alignement 
    f = imfuse(img,img_unreg);
    h1 = figure; imshow(f,'InitialMagnification','fit');
    title('before-registration');
    figure_name = [writedir datafiles{filei}(1:end-4) '_before2' ];
    saveas(h1, figure_name, 'jpeg') %Save figure
    
    f = imfuse(img,img_reg);
    h2 = figure; imshow(f,'InitialMagnification','fit');
    title('after-registration');

    figure_name = [writedir datafiles{filei}(1:end-4) '_adjusted2' ];
    saveas(h2, figure_name, 'jpeg') %Save figure
    
    %% Save the image and the transformation
    
    imwrite(img_reg,[writedir datafiles{filei}(1:end-4) '_adjusted2.tif'],'tiff')
    
    outfilename = [writedir datafiles{filei}(1:end-4) '_adjusted2.mat' ];
    save(outfilename, 'mytform','movingPoints','fixedPoints')
    
end