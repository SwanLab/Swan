clc; clear all; close all;

nFeatures = 784;
nTest = 3000;
load("Test.mat");

%Ha de buscar quin dels 784 és el que te un 1 i saber la posició. Llavors
%comprovar si quadra amb el test la casella. Si al test hi ha un 0 què fem?
%si al test hi ha un 1 doncs casella vàlida posem result = 1;

result = zeros(nFeatures,nFeatures);

    for i = 1:nTest 

        j = find(outputs(i, :) >= 0.5);

        %L'hem trobat; Ara volem saber quin era el target per aquesta foto

        k = find(targets(i, :) == 1);

        result(k,j) = result(k,j) + 1;

    end

    quo = 0;
    
    for u = 1:nFeatures
        quo = quo + result(u,u);
    end

    perc = (quo/nTest)*100;


surf(result);
colormap(jet);
disp(['Output value vs expected percentatge is ', num2str(perc), '%']);