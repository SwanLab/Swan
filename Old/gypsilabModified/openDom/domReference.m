function [x,w] = domReference(domain)
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
%|    #    |   FILE       : domReference.m                                |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Discrete quadrature for edges, triangle       |
%|  `---'  |                and tetrahedre of reference                   |
%+========================================================================+

% Particles mesh
if (size(domain.msh.elt,2) == 1)
    error('domReference.m : unavailable case')
    
% Edge mesh
elseif (size(domain.msh.elt,2) == 2)
    if (domain.gss == 1)
        x = 0.5;
        w = 1;

    elseif (domain.gss == 2)
        x = [0.5*(1-1/sqrt(3)) ; 0.5*(1+1/sqrt(3))];
        w = [0.5 ; 0.5];

    elseif (domain.gss == 3)
        x = [0.5*(1-sqrt(3/5)) ; 0.5 ; 0.5*(1+sqrt(3/5))];
        w = [5/18 ; 4/9 ; 5/18];

    elseif (domain.gss == 4)
        a = sqrt(3/7 + 2/7*sqrt(6/5));
        b = sqrt(3/7 - 2/7*sqrt(6/5));
        w1 = (18-sqrt(30))/72;
        w2 = (18+sqrt(30))/72;
        x = [0.5*(1-a) ; 0.5*(1-b) ; 0.5*(1+b) ; 0.5*(1+a)] ;
        w = [w1 ; w2 ; w2 ; w1];

    elseif (domain.gss == 5)
        a = 1/3*sqrt(5 + 2*sqrt(10/7));
        b = 1/3*sqrt(5 - 2*sqrt(10/7));
        w1 = (322-13*sqrt(70))/1800;
        w2 = (322+13*sqrt(70))/1800;
        x = [0.5*(1-a) ; 0.5*(1-b) ; 0.5 ; 0.5*(1+b) ; 0.5*(1+a)] ;
        w = [w1 ; w2 ;64/225 ; w2 ; w1];

    else
        error('domReference.m : unavailable case')
    end
    
% Triangular mesh
elseif (size(domain.msh.elt,2) == 3)
    if (domain.gss == 1)
        x = [1/3  1/3];
        w = 1;
        
    elseif (domain.gss == 3)
        x = [1/6 1/6 
             2/3 1/6
             1/6 2/3];
        w = [1/3 ; 1/3 ; 1/3];
        
    elseif (domain.gss == 7)
        a = (155-sqrt(15))/1200;
        b = (155+sqrt(15))/1200;
        x = [1/3                1/3
            (6-sqrt(15))/21     (6-sqrt(15))/21
            (6-sqrt(15))/21     (9+2*sqrt(15))/21 
            (9+2*sqrt(15))/21   (6-sqrt(15))/21
            (6+sqrt(15))/21     (6+sqrt(15))/21
            (6+sqrt(15))/21     (9-2*sqrt(15))/21
            (9-2*sqrt(15))/21   (6+sqrt(15))/21];
        w = [9/40 ; a ; a ; a ; b ; b ; b];
        
    elseif domain.gss == 12
        A  = 0.063089014491502;
        B  = 0.249286745170910;
        C  = 0.310352451033785;
        D  = 0.053145049844816;
        P1 = 0.025422453185103;
        P2 = 0.058393137863189;
        P3 = 0.041425537809187;
        x(1:domain.gss,1) = [A 1-2*A A     B 1-2*B B     C D 1-C-D 1-C-D C     D];
        x(1:domain.gss,2) = [A A     1-2*A B B     1-2*B D C C     D     1-C-D 1-C-D];
        w = [P1 P1 P1 P2 P2 P2 P3 P3 P3 P3 P3 P3]'*2;
    
    else
        error('domReference.m : unavailable case')
    end
    
% Tetrahedron mesh
elseif (size(domain.msh.elt,2) == 4)
    if (domain.gss == 1)
        x = [1/4 1/4 1/4];
        w = 1;
        
    elseif (domain.gss == 4)
        a = (5-sqrt(5))/20;
        b = (5+3*sqrt(5))/20;
        x = [a a a;
            a a b;
            a b a;
            b a a];
        w = [1/4 1/4 1/4 1/4]';
            
    elseif (domain.gss == 5)
        a = 1/4;
        b = 1/6;
        c = 1/2;
        x = [a a a;
            b b b;
            b b c;
            b c b;
            c b b];
        w = [-4/5 9/20 9/20 9/20 9/20]';
            
    elseif (domain.gss == 15)
        a  = 1/4;
        b1 = (7+sqrt(15))/34; 
        b2 = (7-sqrt(15))/34;
        c1 = (13-3*sqrt(15))/34;
        c2 = (13+3*sqrt(15))/34;
        d  = (5-sqrt(15))/20;
        e = (5+sqrt(15))/20;
        x = [a a a;
            b1 b1 b1;
            b1 b1 c1;
            b1 c1 b1;
            c1 b1 b1;
            b2 b2 b2;
            b2 b2 c2;
            b2 c2 b2;
            c2 b2 b2;
            d d e;
            d e d;
            e d d;
            d e e;
            e d e;
            e e d];
        w1 = 6*(2665-14*sqrt(15))/226800;
        w2 = 6*(2665+14*sqrt(15))/226800;
        w = [48/405 w1 w1 w1 w1 w2 w2 w2 w2 30/567 30/567 30/567 30/567 30/567 30/567]';
        
    else
        error('domReference.m : unavailable case')
    end
    
% Unknown type
else
    error('domReference.m : unavailable case')
end
end
