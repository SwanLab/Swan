%% visualizaionSVD.m

load('SVD_Results.mat');

% PLot the V columns grouped in 10
step=10;
Nwindow=ceil(size(V,2)/step);
idx=1;

for j=1:Nwindow
    figure('Position',[75 100 1400 600]);
    tiledlayout(2,5,'TileSpacing','compact','Padding','compact');
    for i=1:step
        ax=nexttile;
        plot(radii,V(:,idx), 'LineWidth', 1.5);
        xlabel('r');
        ylabel("V(:,"+idx);
        title("V-"+ idx+" ; r="+radii(1,idx));
        grid on
        idx=idx+1;
    end
end

