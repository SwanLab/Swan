classdef Tutorial14_parametricMeshJoiner2D < handle

    properties (Access = private)
        nSubdomains
        r
        xmin
        xmax
        ymin
        ymax
        cx
        cy
        Nr
        Ntheta
        tolSameNode
        shape
    end

    methods (Access = public)

        function obj = Tutorial14_parametricMeshJoiner2D()
            close all
            obj.init()
            [mSbd,refMesh]      = obj.createSubDomainMeshes(); 
            [mD,~,~,~,discMesh] = obj.mergeSubdomainMeshes(mSbd,refMesh);
            mD.plot

            % refMesh is the mesh from 1st domain, just need info of number of nodes,
            % as we assume all domains have the same number of nodes so far.

            % I leave this comment here. If the merge for other meshes 
            % coming from gid doesn't work, it could be a uniquetol missing
            % (instead of unique)in DomainMeshComputer updateCoord
        end

    end

    methods (Access = private)

        function init(obj)
            obj.nSubdomains  = [15 5]; %nx ny
            obj.tolSameNode  = 1e-10;
            rng(0);                 % seed value for random distribution
            maxr             = 0.1; % max value of radius
            minr             = 0.9; % min value of radius
            obj.r            = (maxr - minr) * rand(obj.nSubdomains(2),obj.nSubdomains(1)) + minr;
            obj.xmin         = -1;
            obj.xmax         = 1;
            obj.ymin         = -1;
            obj.ymax         = 1;
            obj.cx           = 0;  % center of the inclusion
            obj.cy           = 0;  % center of the inclusion
            obj.Nr           = 7;  % number of dividions along the radial direction
            obj.Ntheta       = 14; % dividions in the angular direction
            obj.shape        = 'Circle';
        end

        function  [mSbd,refMesh] = createSubDomainMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);
            Lx = obj.xmax-obj.xmin;
            Ly = obj.ymax-obj.ymin;
            for jDom = 1:nY
                for iDom = 1:nX
                    meshOrigin = createRectangularMeshWithHole(...
                                   obj.r(jDom,iDom),obj.xmax,obj.xmin,obj.ymax,obj.ymin,...
                                   obj.Nr,obj.Ntheta,obj.cx,obj.cy,obj.shape);
                    mIJ        = obj.shiftMeshCoord(meshOrigin,Lx,Ly,iDom,jDom);
                    mSbd{jDom,iDom} = mIJ;
                end
            end
            refMesh = mSbd{1,1};
        end

        function m = shiftMeshCoord(obj,meshOrigin,Lx,Ly,iDom,jDom)
            coord0       = meshOrigin.coord;
            s.coord(:,1) = coord0(:,1)+Lx*(iDom-1);
            s.coord(:,2) = coord0(:,2)+Ly*(jDom-1);
            s.connec     = meshOrigin.connec;
            m            = Mesh.create(s);
        end

        function [mD,iC,lG,iCR,discMesh] = mergeSubdomainMeshes(obj,mSbd,refMesh)
            s.nsubdomains   = obj.nSubdomains;
            s.meshReference = refMesh;
            s.tolSameNode   = obj.tolSameNode;
            s.meshSbd       = mSbd;
            m = MeshJoiner(s);
            [mD,iC,~,lG,iCR,discMesh] = m.create();
        end

    end

end
