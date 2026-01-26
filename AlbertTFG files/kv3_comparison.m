clc;
clear;
close all;

NNname = "AlbertTFG files/comparisons/k case 3/100x100/deg1 test50 lambda0 rate 0_01 layers256x1_16 err0.0058.mat";
plotName = 'pol = 1; test = 50; learning = 0.01; hidden = 256x(1,16); error = 0.0058';


filePathShort = fullfile('AlbertTFG files/comparisons/k case 3/20x20/K_data case 3.csv');
filePathLong = fullfile('AlbertTFG files/comparisons/k case 3/100x100/K_data case 3 100x100.csv');

opt = load(NNname);
opt = opt.opt;

k3_short = readmatrix(filePathShort);
k3_long  = readmatrix(filePathLong);


%% Short

error = [];

for i = 1:size(k3_short,1)
    predicted = opt.computeOutputValues(k3_short(i,1));
    error = cat(1, error, k3_short(i,2:end)-predicted );
end

error_short = error;

nx = 3;
ny = 2;


figure
for i = 1:size(k3_short,1)

    aux = [];
    for j = 1:8
        aux = cat(1, aux, error_short(   i,8*(j-1)+1 : 8*(j)   ));
    end
    maxDiff = max(max(abs(aux)));
    subplot(nx, ny, i);
        
    h = bar3c(abs(aux));
    zlabel("abs(diff)")
    title("r = "+ num2str(k3_short(i,1)))
    colormap('turbo');
    clim([min(min(abs(aux))), max(max(abs(aux)))])
    colorbar()
        if maxDiff >10^-5
        subtitle("max diff = " + num2str(maxDiff),'Color', 'red')
    else
        subtitle("max diff = " + num2str(maxDiff),'Color', 'green')
    end


    

    
end
sgtitle(sprintf(append('Error real vs predicted\n', plotName)))
pos = get(gcf, 'Position');
set(gcf, 'Position',pos+[-100 -400 500 500])

%% Long


error = [];

for i = 1:size(k3_long,1)
    predicted = opt.computeOutputValues(k3_long(i,1));
    error = cat(1, error, k3_long(i,2:end)-predicted );
end
error_long = error;


r_long = k3_long(:,1)';
arrayPos = 1:1:(size(k3_short,2)-1);

[A, R] = meshgrid(arrayPos, r_long);

figure
s = surf(R,A,abs(error_long));
s.EdgeColor = 'none';
zscale log
colorbar
colormap('turbo');
clim([min(min(abs(error_long))), max(max(abs(error_long)))])
view(0,90)

xlabel("Radius of inclusion")
ylabel("Position in array")
zlabel("Difference between predicted and real")

title("Abs of real vs predicted")
subtitle(plotName)



