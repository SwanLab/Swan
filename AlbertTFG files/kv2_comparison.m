clc;
clear;
close all;

NNname = "AlbertTFG files/comparisons/k case 2/100x100/k case2 100x100.mat";
plotName = 'pol = 1; test = 50; learning = 0.01; hidden = 256x(1,16)';


filePath = fullfile('AlbertTFG files/comparisons/k case 2/100x100/K_data case 2 100x100.csv');


opt = load(NNname);
opt = opt.opt;

k2 = readmatrix(filePath);



%% Short

% error = [];
% 
% for i = 1:size(k2,1)
%     predicted = opt.computeOutputValues(k2(i,1:9));
%     error = cat(1, error, k2(i,10:end)-predicted );
% end
% 
% nx = 2;
% ny = 3;
% 
% 
% figure
% for i = 1:size(k2,1)/8
% 
%     aux = error( 8*(i-1)+1:8*i ,:  );
%     maxDiff = max(max(abs(aux)));
%     subplot(nx, ny, i);
% 
%     h = bar3c(abs(aux));
%     zlabel("abs(diff)")
%     title("r = "+ num2str(k2(i,1)))
%     colormap('turbo');
%     clim([min(min(abs(aux))), max(max(abs(aux)))])
%     colorbar()
%         if maxDiff >10^-5
%         subtitle("max diff = " + num2str(maxDiff),'Color', 'red')
%     else
%         subtitle("max diff = " + num2str(maxDiff),'Color', 'green')
%     end
% 
% 
% 
% 
% 
% end
% sgtitle(sprintf(append('Error real vs predicted\n', plotName)))
% pos = get(gcf, 'Position');
% set(gcf, 'Position',pos+[-100 -400 500 500])



%% Long


error = [];

for i = 1:size(k2,1)
    predicted = opt.computeOutputValues(k2(i,1));
    error = cat(1, error, k2(i,10:end)-predicted );
end
error_long = error;

k2_order = []

error_formated = [];
for i=1:size(k2,1)/8
    err = [];
    k2_o = [];
    for j = 1:8
        err = cat(2, err, error_long(8*(i-1)+j,:) ) ;
        k2_o = cat(2, k2_o, k2(8*(i-1)+j,10:end) );
    end
    error_formated = cat(1, error_formated, err);
    k2_order = cat(1, k2_order, k2_o);
end



r_long = k2(:,1)';
r_long = unique(r_long);
arrayPos = 1:1:64;

[A, R] = meshgrid(arrayPos, r_long);

figure
s = surf(R,A,abs(error_formated));
s.EdgeColor = 'none';
zscale log
colorbar
colormap('turbo');
clim([min(min(abs(error_formated))), max(max(abs(error_formated)))])
view(0,90)

xlabel("Radius of inclusion")
ylabel("Position in array")
zlabel("Difference between predicted and real")

title("Abs of real vs predicted")
subtitle(plotName)



