% Author: Chloros2 <chloros2@gmx.de>
% Created: 2019-02-11
%
% Copyright (C) 2018 Chloros2 <chloros2@gmx.de>
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see
% <https://www.gnu.org/licenses/>.
%

clear all;
addpath('ffmatlib');

[p,b,t,nv,nbe,nt,labels]=ffreadmesh('capacitor_2d.msh');

figure;
ffpdeplot(p,b,t,'Boundary','on', ...
                'Title','All boundaries / edges');
figure;
ffpdeplot(p,b,t,'Boundary','on','BDLabels',labels,'BDShowText','on', ...
                'Title','With Text - using the labels variable');

figure;
ffpdeplot(p,b,t,'Boundary','on','BDLabels',labels,'BDShowText','on', ...
                'BDTextSize',16,'BDTextWeight','bold', ...
                'Title','Text size and weight');

figure;
ffpdeplot(p,b,t,'Boundary','on','BDLabels',[3,4,5], ...
                'BDColors',[[1,0,0];[0,1,0];[0,0,1]],'BDShowText','on', ...
                'Title','Three column matrix of RGB triplets (3x3)');

plotcols =  [[1,0,1];[153,101,21]./255];
figure;
ffpdeplot(p,b,t,'Boundary','on','BDLabels',[3,4], ...
                'BDColors',plotcols,'BDShowText','on', ...
                'Title','Three column matrix of RGB triplets (3x2)');

plotcols = rand(numel(labels),3);
figure;
ffpdeplot(p,b,t,'Boundary','on','BDLabels',labels, ...
                'BDColors',plotcols,'BDShowText','on', ...
                'Title','All boundaries are colored by random');

figure;
ffpdeplot(p,b,t,'Boundary','on','BDColors','b', ...
                'Title','A single character like b specifies color for all boundaries');

figure;
ffpdeplot(p,b,t,'Boundary','on','BDColors','green', ...
                'Title','A name like green specifies the color for all boundaries');
figure;
ffpdeplot(p,b,t,'Boundary','on','BDLabels',[3,4],'BDColors','g','BDShowText','on', ...
                'Title','A single character like r works also for a subset of labels');
