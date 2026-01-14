clc
clear all
close all


load ('XSNAPdataRVEe_3.mat','XSNAP')
X_data_show = XSNAP ;
% 
% ALEXP =1 ; 
% if ALEXP == 1
%     for i = 1:size(X_data_show,2)
%         n(i) = bin2dec(num2str(X_data_show(:,i)'));
%     end
%     plot(n)
% end


% r = 5 ;
% ang = 0:0.4:2*pi ;
% x  =r*cos(ang); y = r*sin(ang);
% X_data_show = [x;y];
D = L2_distance(X_data_show, X_data_show, 1);

%To run Isomap with K = 7 neighbors/point, and produce embeddings
%in dimensions 1, 2, ..., 10, type these commands:
options = struct('dims',1:15,'overlay',1,'comp',1,'display',1,'verbose',1);
%  options.dims = 1:10;
n_fcn ='k'; 'epsilon' ;
n_size =2;

[Y, R, E] = Isomap(D, n_fcn, n_size, options);

% nfigure = 100 ;
%  
% figure(nfigure)
% hold on
% fn = ['PICTURES/pict_',num2str(9),'.jpg'];
% dispim(fn, nfigure, [0 0])

% 
%  if ~isempty(twod)
%         figure(100);
%         hold on;
%         plot(Y.coords{twod}(1,:), Y.coords{twod}(2,:), 'ro');
%        for iii = 1:length(Y.coords{twod}(1,:))
%            text(Y.coords{twod}(1,iii), Y.coords{twod}(2,iii), num2str(iii));
%        end
%         if (overlay == 1)
%             gplot(E(Y.index, Y.index), [Y.coords{twod}(1,:); Y.coords{twod}(2,:)]');
%             title('Two-dimensional Isomap embedding (with neighborhood graph).');
%         else
%             title('Two-dimensional Isomap.');
%         end
%         hold off;
%     end




% ISOMAP   Computes Isomap embedding using the algorithm of
%             Tenenbaum, de Silva, and Langford (2000).
%
% [Y, R, E] = isomap(D, n_fcn, n_size, options);
%
% Input:
%    D = N x N matrix of distances (where N is the number of data points)
%    n_fcn = neighborhood function ('epsilon' or 'k')
%    n_size = neighborhood size (value for epsilon or k)
%
%    options.dims = (row) vector of embedding dimensionalities to use
%                        (1:10 = default)
%    options.comp = which connected component to embed, if more than one.
%                        (1 = largest (default), 2 = second largest, ...)
%    options.display = plot residual variance and 2-D embedding?
%                        (1 = yes (default), 0 = no)
%    options.overlay = overlay graph on 2-D embedding?
%                        (1 = yes (default), 0 = no)
%    options.verbose = display progress reports?
%                        (1 = yes (default), 0 = no)
%
% Output:
%    Y = Y.coords is a cell array, with coordinates for d-dimensional embeddings
%         in Y.coords{d}.  Y.index contains the indices of the points embedded.
%    R = residual variances for embeddings in Y
%    E = edge matrix for neighborhood graph