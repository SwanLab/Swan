clc;
clear;
close all;


filePathRecognition = fullfile('AlbertTFG files/PatternRecognition.csv');
filePathNormalized = fullfile('AlbertTFG files/PatternNormalized.csv');

patternRecognition = readmatrix(filePathRecognition);
patternNormalized = readmatrix(filePathNormalized);


recog = [];
norm  =[];

for i = 1:9
    recog = cat(3, recog, patternRecognition(8*(i-1)+1:8*(i), 2:end));
    norm = cat(3, norm, patternNormalized(8*(i-1)+1:8*(i), 2:end));
end

diffs = {};
diff = [];

for i = 1:9
    aux = {};
    au = [];
    for j = 1:9
        aux = cat(2, aux, {norm(:,:,i)-norm(:,:,j)});
        au = cat(2, au, norm(:,:,i)-norm(:,:,j));
    end
    diffs = cat(1, diffs, aux);
    diff = cat(1, diff, au);
end

figure
b=bar3(abs(diff));
colorbar
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end

title("Abs")
xlabel("Position i in array")
ylabel("Position j in array")