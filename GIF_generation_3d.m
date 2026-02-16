%% Script to generate gifs for 3D Simulations
clc; clear;

%% Inputs
NameFile = 'LevelSet_Standard_0_3D_v2.*.png';   
NameOutput = 'GIF_LevelSet_Standard_0_3D_v3.gif';
output= fullfile("3D_files/3D_Gifs/",NameOutput);

fps = 28;               
delay = 1/fps;
scale= 0.75;

%% File lecture
files = dir(fullfile('3D_files/PNG_Images/',NameFile));
files = sort({files.name});   

%% SELECT REGION TO CAPTURE
img = imread(files{1});
imshow(img)
rect = getrect;

%% GIF GENERATION
for k = 1:length(files)
    img = imread(files{k});
    img_crop = imcrop(img, rect);
    img_small = imresize(img_crop, scale);

    [A,map] = rgb2ind(img_small,256);

    if k == 1
        imwrite(A,map,output,'gif', ...
            'LoopCount',inf, ...
            'DelayTime',delay);
    else
        imwrite(A,map,output,'gif', ...
            'WriteMode','append', ...
            'DelayTime',delay);
    end
end

disp('GIF created')
