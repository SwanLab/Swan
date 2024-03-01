clear; clc; close all;

load("Test.mat");

targ = targets(1,:);
outp = outputs(1,:);

% Reshape del vector a una matriz de 28x28
img_target = reshape(targ, 28, 28);
img_output = reshape(outp, 28, 28);

% Mostrar la imagen
figure(1)
imshow(img_target);
figure(2)
imshow(img_output);