%% Script to generate gifs for 3D Simulations

%% Inputs
NameFile = 'LevelSet_45_v2_3D.*.png';   
NameOutput = 'GIF_LevelSet_45_3D_v2.gif';
output= fullfile("3D_files/3D_Gifs/",NameOutput);

fps = 28;               
delay = 1/fps;

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
    [A,map] = rgb2ind(img_crop,256);

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
