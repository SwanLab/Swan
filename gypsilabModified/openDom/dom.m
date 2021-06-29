classdef dom
%+========================================================================+
%|                                                                        |
%|              OPENDOM - LIBRARY FOR NUMERICAL INTEGRATION               |
%|           openDom is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal & Francois Alouges (c) 2017-2018.          |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%|             francois.alouges@polytechnique.edu                         |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : dom.m                                         |
%|    #    |   VERSION    : 0.61                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Domain class definition                       |
%|  `---'  |                                                              |
%+========================================================================+

properties 
    gss = 1;       % NUMBER OF INTEGRATION POINTS
    msh = [];      % MESH OF THE DOMAIN
end

methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function domain = dom(varargin)
        if nargin > 0
            % Mesh
            domain     = dom;
            domain.msh = varargin{1};
            
            % Gaussian formula
            if (nargin == 2)
                domain.gss = varargin{2};
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot(varargin)
        domain = varargin{1};
        spc    = '.b';
        if (nargin == 2)
            spc = varargin{2};
        end
        X = domain.qud;
        plot3(X(:,1),X(:,2),X(:,3),spc)
    end
    
    function plotNrm(varargin)
        domain = varargin{1};
        spc    = 'b';
        if (nargin == 2)
            spc = varargin{2};
        end
        Xqud = domain.qud;
        Vnrm = domain.qudNrm;
        quiver3(Xqud(:,1),Xqud(:,2),Xqud(:,3),Vnrm(:,1),Vnrm(:,2),Vnrm(:,3),spc);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LENGTH
    function l = length(domain)
        l = size(domain.qud,1);
    end
    
    % SIZE
    function s = size(varargin)
        s = size(varargin{1}.qud);
        if (nargin == 2)
           s = s(varargin{2}); 
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QUADRATURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GAUSSIAN QUADRATURE
    function [X,W,elt2qud] = qud(domain)
        [X,W,elt2qud] = domQuadrature(domain);
    end
        
    % GAUSSIAN NORMALES
    function N = qudNrm(domain)
        x = domReference(domain);
        N = zeros(size(x,1)*size(domain.msh.elt,1),size(domain.msh.vtx,2));
        for j = 1:size(x,1)
            idx      = (j:size(x,1):size(x,1) * size(domain.msh.elt,1))';
            N(idx,:) = domain.msh.nrm;
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTEGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INTEGRAL FCT
    function I = integral(varargin)
        switch nargin
            case 2
                I = domIntegral2(varargin);                
            case 3
                I = domIntegral3(varargin);                
            case 4
                I = domIntegral4(varargin);
            case 5
                I = domIntegral5(varargin);
            case 6
                I = domIntegral6(varargin);
            case 7
                I = domIntegral7(varargin);
            otherwise
                error('dom.m : unavailable case')
        end
    end
    
    % SINGULAR REGULARIZAION
    function S = regularize(varargin)
        mesh = varargin{end}.msh;
        if (size(mesh,2) == 3)
            S = domRegularize3D(varargin);
        elseif (size(mesh,2) == 2) && is2d(mesh)
            S = domRegularize2D(varargin);
        end
    end
    
    % INTERPOLATION
    function I = interpolate(varargin)
        domain = varargin{1};
        u      = varargin{2};
        v      = varargin{3};
        M      = integral(domain,u,u);
        Fv     = integral(domain,u,v);
        if (nargin == 4)
            f = varargin{4};
        else
            f = 1;
        end
        if iscell(Fv)
            I{1} = M \ (Fv{1} * f);
            I{2} = M \ (Fv{2} * f);
            I{3} = M \ (Fv{3} * f);
        else
            I = M \ (Fv * f);
        end
    end
        
    % L2 AND H1 ERRORS
    function err = diff(domain, fe, sol, ref, type )
        err = domDifference(domain, fe, sol, ref, type );
    end
end
end
