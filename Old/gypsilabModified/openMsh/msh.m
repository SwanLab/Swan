classdef msh
%+========================================================================+
%|                                                                        |
%|                 OPENMSH - LIBRARY FOR MESH MANAGEMENT                  |
%|           openMsh is part of the GYPSILAB toolbox for Matlab           |
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
%|    #    |   FILE       : msh.m                                         |
%|    #    |   VERSION    : 0.53                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 21.06.2019                                    |
%| ( === ) |   SYNOPSIS   : Mesh class definition                         |
%|  `---'  |                                                              |
%+========================================================================+

properties 
    vtx = [];      % VERTEX COORDINATES (3 dimensions)
    elt = [];      % ELEMENTS LIST (particles, edges, triangles or tetrahedron)
    col = [];      % ELEMENT COLOR (group number)
end

methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function mesh = msh(varargin)
        % Read file
        if (length(varargin) == 1) && ischar(varargin{1})
            mesh = msh;
            file = varargin{1}; 
            ext  = file(end-2:end);
            if strcmp(ext,'msh')
                [mesh.vtx,mesh.elt,mesh.col] = mshReadMsh(file);
            elseif strcmp(ext,'ply')
                [mesh.vtx,mesh.elt] = mshReadPly(file);
            elseif strcmp(ext,'stl')
                [mesh.vtx,mesh.elt] = mshReadStl(file);
            elseif strcmp(ext,'vtk')
                [mesh.vtx,mesh.elt] = mshReadVtk(file);
            elseif strcmp(file(end-3:end),'mesh')
                [mesh.vtx,mesh.elt,mesh.col] = mshReadMesh(file);
            else
                error('msh.m : unavailable case')
            end
            if isempty(mesh.col)
                mesh.col = zeros(size(mesh.elt,1),1);
            end
            
        % Inout only vertex
        elseif (length(varargin) == 1) && isnumeric(varargin{1})
            mesh.vtx = varargin{1};
            mesh.elt = (1:size(mesh.vtx,1))';
            mesh.col = zeros(size(mesh.elt,1),1);
        
        % Input vertex and elements   
        elseif (length(varargin) == 2)
            mesh.vtx = varargin{1};
            mesh.elt = varargin{2};
            mesh.col = zeros(size(mesh.elt,1),1);
        
        % Input vertex, elements and colours     
        elseif (length(varargin) == 3)
            mesh.vtx = varargin{1};
            mesh.elt = varargin{2};
            mesh.col = varargin{3};
        end
        
        % Clean mesh
        mesh = clean(mesh);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLEAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function mesh = clean(varargin)
        mesh = varargin{1};
        if (nargin==1)
            mesh = mshClean(mesh,[]);
        else
            mesh = mshClean(mesh,varargin{2});            
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot(varargin)
        mesh = varargin{1};
        if (nargin == 1)
            mshPlot(mesh,[])
        else
            V = varargin{2};
            if (numel(V) == size(mesh.elt,1)) && (~ischar(V))
                mesh.col = V;
                mshPlot(mesh,[])
            else
                mshPlot(mesh,V)
            end
        end
    end
    
    function plotNrm(varargin)
        mesh = varargin{1};
        spc  = 'b';
        if (nargin == 2)
            spc = varargin{2};
        end            
        Xctr = mesh.ctr;
        Vnrm = mesh.nrm;
        quiver3(Xctr(:,1),Xctr(:,2),Xctr(:,3),Vnrm(:,1),Vnrm(:,2),Vnrm(:,3),spc);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LENGTH
    function l = length(mesh)
        l = size(mesh.elt,1);
    end
    
    % SIZE
    function s = size(varargin)
        s = size(varargin{1}.elt);
        if (nargin == 2)
           s = s(varargin{2}); 
        end
    end
    
    % STEP
    function l = stp(mesh)
        mesh = mesh.edg;
        l    = mesh.vtx(mesh.elt(:,2),:) - mesh.vtx(mesh.elt(:,1),:);
        l    = sqrt(sum(l.^2,2));
        l    = [min(l) max(l) mean(l) std(l)];
    end 
    
    % 2 DIMENSIONS
    function b = is2d(mesh)
        b = (max(abs(mesh.vtx(:,3))) < 1e-12);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ELEMENT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CENTER
    function X = ctr(mesh)
        X = zeros(size(mesh.elt,1),size(mesh.vtx,2));
        for i = 1:size(mesh.elt,2)
            X = X + (1/size(mesh.elt,2)) .* mesh.vtx(mesh.elt(:,i),:);
        end
    end
    
    % ND-VOLUME
    function V = ndv(mesh)
        V = mshNdvolume(mesh);
    end
    
    % NORMALS
    function N = nrm(mesh)
        T = mesh.tgt;
        if (size(mesh,2) == 3)
            N = cross(T{1},T{2},2);
            N = N ./ (sqrt(sum(N.^2,2)) * [1 1 1]);
        elseif (size(mesh,2) == 2) && is2d(mesh)
            N = T * [0 -1 0 ; 1 0 0 ; 0 0 1]';
        else
            error('msh.m : unavailable case')
        end
    end
    
    % EDGES NORMALS
    function Nu = nrmEdg(mesh)
        if (size(mesh,2) == 3)
            Nu = cell(1,3);
            for i = 1:3
                Nu{i} = cross(mesh.tgt{i},mesh.nrm,2);
            end
        elseif (size(mesh,2) == 2)
            T     = mesh.tgt; 
            Nu{1} = T;
            Nu{2} = -T;
        else
            error('msh.m : unavailable case')
        end
    end
    
    % TANGENTS
    function T = tgt(mesh)
        if (size(mesh,2) == 3)
            T = cell(1,3);
            for i = 1:3
                ip1  = mod(i,3)+1;
                ip2  = mod(ip1,3)+1;
                A    = mesh.vtx(mesh.elt(:,ip1),:);
                B    = mesh.vtx(mesh.elt(:,ip2),:);
                T{i} = (B-A)./(sqrt(sum((B-A).^2,2))*[1 1 1]);
            end
        elseif (size(mesh,2) == 2)
            A = mesh.vtx(mesh.elt(:,1),:);
            B = mesh.vtx(mesh.elt(:,2),:);
            T = (B-A)./(sqrt(sum((B-A).^2,2))*[1 1 1]);
        else
            error('msh.m : unavailable case')
        end
    end
    
    % SWAP
    function mesh = swap(varargin)
        mesh = varargin{1};
        Ielt = 1:size(mesh.elt,1);
        if (nargin == 2)
            Ielt = varargin{2};
        end
        if (size(mesh,2) == 2)
            mesh.elt(Ielt,:) = mesh.elt(Ielt,[2 1]);
        elseif (size(mesh,2) == 3)
            mesh.elt(Ielt,:) = mesh.elt(Ielt,[2 1 3]);
        else
            error('msh.m : unavailable case')
        end
    end
    
    % COLOURS
    function mesh = color(mesh,c)
        mesh.col(:) = c;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBMESH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SUBMESHING
    function mesh = sub(mesh,Ielt)
        mesh.elt = mesh.elt(Ielt,:);
        mesh.col = mesh.col(Ielt);
        mesh     = clean(mesh);
    end  
    
    % FACES
    function [mesh,elt2fce] = fce(mesh)
        [mesh,elt2fce] = mshFace(mesh);
    end
    
    % EDGES
    function [mesh,elt2edg] = edg(mesh)
        [mesh,elt2edg] = mshEdge(mesh);
    end
    
    % PARTICLES
    function [mesh,elt2prt] = prt(mesh)
        [mesh,elt2prt] = mshParticle(mesh);
    end
    
    % BOUNDARY
    function mesh = bnd(mesh)
        mesh = mshBoundary(mesh);
    end
        
    % MIDPOINT
    function [mesh,Ir] = midpoint(varargin)
        mesh = varargin{1};
        if (nargin == 1)
            I = (1:length(mesh))';
        else
            I = varargin{2};
        end
        [mesh,Ir] = mshMidpoint(mesh,I);
    end
    
    % REFINE WITH MIDPOINT
    function mesh = refine(varargin)
        mesh = varargin{1};
        ref  = sum(mesh.ndv);
        if (nargin == 1)
            mesh = midpoint(mesh);
        else
            mesh = mshRefine(mesh,varargin{2});
        end
        sol = sum(mesh.ndv);
        if (norm(ref-sol)/norm(ref) > 1e-12)
            error('msh.m : unavailable case')
        end
    end
    
    % HIERARCHICAL TREE
    function tree = tree(varargin)
        mesh = varargin{1};
        if (nargin == 1)
            typ = 'octree';
            Nlf = 1;
            fig = 0;
        elseif (nargin == 2)
            typ = varargin{2};
            Nlf = 1;
            fig = 0;
        elseif (nargin == 3)
            typ = varargin{2};
            Nlf = varargin{3};
            fig = 0;
        elseif (nargin == 4)
            typ = varargin{2};
            Nlf = varargin{3};
            fig = varargin{4};
        end
        tree = mshTree(mesh,typ,Nlf,fig);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALGEBRA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIGNATURE
    function M = sgn(mesh)
        M = cell(1,size(mesh.elt,2));
        for i = 1:length(M)
            M{i} = mesh.vtx(mesh.elt(:,i),:);
        end
        M = sort(cell2mat(M),2);
        M = [M,mesh.ctr];
        M = single(M);
    end
    
    % UNICITY
    function [meshC,IA,IC] = unique(mesh)
        M         = sgn(mesh);
        [~,IA,IC] = unique(M,'rows','stable');
        meshC     = mesh.sub(IA);
    end
    
    % INTERSECTION
    function [meshC,IA,IB] = intersect(meshA,meshB)
        A         = sgn(meshA);
        B         = sgn(meshB);
        [~,IA,IB] = intersect(A,B,'rows','stable');
        meshC     = meshA.sub(IA);
    end
    
    % DIFFERENCE
    function [meshC,IA] = setdiff(meshA,meshB)
        A      = sgn(meshA);
        B      = sgn(meshB);
        [~,IA] = setdiff(A,B,'rows','stable');
        meshC  = meshA.sub(IA);
    end
    
    % UNION
    function [meshC,IA,IB] = union(meshA,meshB)
        A         = sgn(meshA);
        B         = sgn(meshB);
        [~,IA,IB] = union(A,B,'rows','stable');
        meshA     = meshA.sub(IA);
        meshB     = meshB.sub(IB);
        vtxC      = [meshA.vtx ; meshB.vtx];
        eltC      = [meshA.elt ; meshB.elt + size(meshA.vtx,1)];
        colC      = [meshA.col ; meshB.col];
        meshC     = msh(vtxC,eltC,colC);
    end    
    
    % ISMEMBER
    function [IA,IB] = ismember(meshA,meshB)
        A       = sgn(meshA);
        B       = sgn(meshB);
        [IA,IB] = ismember(A,B,'rows');
    end
    
    % FUNCTION
    function mesh = fct(mesh,fct)
        mesh.vtx = fct(mesh.vtx);
    end
    
    % SHUFFLE
    function mesh = shuffle(mesh,varargin)
        Nvtx = size(mesh.vtx,1);        
        if (nargin == 2)
            RPV = (1:Nvtx)';
        else
            RPV = randperm(Nvtx);
        end
        mesh.vtx(RPV,:) = mesh.vtx;
        mesh.elt        = RPV(mesh.elt);         
        RPE      = randperm(length(mesh));
        mesh.elt = mesh.elt(RPE,:);
        mesh.col = mesh.col(RPE,:);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSFORMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TRANSLATION
    function mesh = translate(mesh,U)
        mesh.vtx = mesh.vtx + ones(size(mesh.vtx,1),1)*U;
    end
    
    % ROTATION
    function mesh = rotate(mesh,U,phi)
        N        = U./norm(U);
        mesh.vtx = cos(phi) * mesh.vtx + ...
            (1-cos(phi)) .* ((mesh.vtx*N')*N) + ...
            sin(phi) .* cross(ones(size(mesh.vtx,1),1)*N,mesh.vtx,2);
    end
    
    % SPLIT
    function [mesh1,mesh2] = split(mesh,X0,U)
        [mesh1,mesh2] = mshSplit(mesh,X0,U);
    end
end
end
