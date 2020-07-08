classdef ray
%+========================================================================+
%|                                                                        |
%|           OPENRAY - LIBRARY FOR TRI-DIMENSIONAL RAY TRACING            |
%|           openRay is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : ray.m                                         |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2018                                    |
%| ( === ) |   SYNOPSIS   : Ray tracing class definition                  |
%|  `---'  |                                                              |
%+========================================================================+

properties    
    pos = {};            % SPATIAL POSITION AT EACH STEP 
    dir = [];            % DIRECTION AT FINAL STEP
    iel = [];            % COLLISION ELEMENT INDICE AT EACH STEP     
    msh = [];            % RAY TRACING SPACE
end

methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function ray = ray(mesh,Xsrc,dir)
        if (numel(dir) == 1)
            Nray       = dir;
            ray.pos{1} = ones(Nray,1) * Xsrc;
            ray.dir    = raySource(Nray);
        else
            Nray       = size(dir,1);
            ray.pos{1} = ones(Nray,1) * Xsrc;            
            ray.dir    = dir ./ (sqrt(sum(dir.^2,2)) * [1 1 1]);
        end
        ray.iel = zeros(Nray,1);
        ray.msh = mesh;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot(varargin)
        ray = varargin{1};
        if (nargin == 1)
            ord = length(ray.pos)-1;
        else
            ord = min(varargin{2},length(ray.pos)-1);
        end
        rayPlot(ray,ord)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LENGTH
    function s = length(ray)
        s = size(ray.dir,1);
    end
    
    % SIZE
    function s = size(ray)
        s = [size(ray.dir,1) length(ray.pos)];
    end
    
    % DISTANCE
    function dst = dst(ray)
       dst = zeros(length(ray),1);
       for i = 2:length(ray.pos)
           tmp = ray.pos{i} - ray.pos{i-1};
           dst = dst + sqrt(sum(tmp.^2,2));
       end
    end
    
    % IMAGES
    function [img,nrg] = image(ray,mic,rad,rMax)
        [img,nrg] = rayImage(ray,mic,rad,rMax);
    end

    % FIR 8
    function fir = bank(ray,N,fs)
        load('odeon.mat')
        f = 0.5.*(odeon.frq(1:end-1)+odeon.frq(2:end));
        f = [2/fs.*f 1];
        fir(:,1) = firls(N,[0 f(1) f(1) 1],[1 1 0 0]);
        for i = 1:length(f)-1
            fir(:,i+1) = firls(N,[0 f(i) f(i) f(i+1) f(i+1) 1],[0 0 1 1 0 0]);
        end
        fir(:,end) = firls(N,[0 f(end-1) f(end-1) 1],[0 0 1 1]);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODIFIER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RAY-TRACER
    function ray = tracer(varargin)
        ray = varargin{1};
        if (nargin == 1)
            ord  = length(ray.pos);
            rMax = 1e8;
        elseif (nargin == 2)
            ord  = varargin{2};
            rMax = 1e8;
        elseif (nargin == 3)
            ord  = varargin{2};
            rMax = varargin{3};
        end
        ray = rayTracer(ray,ord,rMax);
    end
    
    % SUB-RAY   
    function ray = sub(ray,I)       
       for i = 1:length(ray.pos)
           ray.pos{i} = ray.pos{i}(I,:);
       end
       ray.dir = ray.dir(I,:);
       ray.iel = ray.iel(I,:);       
    end   
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INDICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MEASURE
    function [I,src] = measure(varargin)
        ray  = varargin{1};
        Xmes = varargin{2};
        rad  = varargin{3};
        if (nargin < 4)
            rMax = 1e8;
        else
            rMax = varargin{4};
        end
        [I,src] = rayMeasure(ray,Xmes,rad,rMax);
    end  
    
    % SPHERE INTERSECTION
    function I = inSphere(ray,Xmes,rad)
        I = rayInSphere(ray,Xmes,rad);
    end      
    
    % DYNAMIC
    function I = dynamic(ray)
        I = find(sum(ray.dir.^2,2)~=0);
    end  
end
end
