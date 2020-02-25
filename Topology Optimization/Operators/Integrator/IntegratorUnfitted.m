classdef IntegratorUnfitted < Integrator
    
    properties (GetAccess = public, SetAccess = protected)
        meshUnfitted
        meshBackground
    end
    
    properties (Access = private)
        backgroundInterp
        unfittedInterp
        quadrature
        unfittedQuad
        integratorCut
    end
    
    methods (Access = public, Abstract)
        integrate(obj)
    end
    
    methods (Access = public)
        
        function obj = IntegratorUnfitted(cParams)
            obj.init(cParams);
            obj.meshUnfitted   = obj.mesh;
            obj.meshBackground = obj.mesh.meshBackground;
            cParams.meshBackground = obj.meshBackground;
            obj.integratorCut = IntegratorCutMesh(cParams);
        end
        
        function A = integrateUnfittedMesh(obj,F,meshUnfitted)
            if exist('meshUnfitted','var')
                %obj.updateMeshes(meshUnfitted);
                %   obj.updateBackgroundMesh();
                %   obj.updateUnfittedMesh();
            end
            A = obj.integrate(F);
            A2 = obj.integratorCut.integrate(F);
            norm(A(:)-A2(:))
        end
        
        
        
    end
    
    methods (Static, Access = public)
        
        
    end
    
    methods (Access = protected)
        
        function shapeValues = evaluateCutShapes(obj,F1)
            obj.createBackgroundInterpolation();
            obj.createUnfittedInterpolation();
            obj.computeThisQuadrature();
            obj.computeUnfittedGaussPoints();
           
            shapes = obj.computeShapes();
            dvolum = obj.computeJacobians();
            Fgauss = obj.interpolateFunctionInGaussPoints(F1,shapes);
            
             
            
            dV = dvolum;
            fdV(1,:,:) = Fgauss.*dV;   
            shapes = permute(shapes,[1 3 2]);
            shapeValues = bsxfun(@times,shapes,fdV);
            shapeValues = sum(shapeValues,3);
            shapeValues = shapeValues';
            
        end
        
        
        function M2 = rearrangeOutputRHS(obj,shapes)
            npnod = obj.meshBackground.npnod;
            nnode = obj.meshBackground.nnode;
            
            M2 = zeros(npnod,1);
            for inode = 1:nnode
                M2 = M2 + accumarray(obj.meshBackground.connec(:,inode),shapes(:,inode),[npnod,1],@sum,0);
            end
        end
        
        function itIs = isLeveSetCuttingMesh(obj)
            itIs = ~isempty(obj.meshUnfitted.backgroundCutCells);
        end
        
        
    end
    
    methods (Static, Access = private)
        
        function xGaus = computePosGP(coord,interp,quad)
            interp.computeShapeDeriv(quad.posgp);
            xGaus = zeros(size(coord,3),quad.ngaus,size(coord,1));
            for igaus = 1:quad.ngaus
                for idime = 1:size(coord,3)
                    xGaus(idime,igaus,:) = coord(:,:,idime)*interp.shape(:,igaus);
                end
            end
        end
        
    end
    
    
    methods (Access = private)
        
        function computeUnfittedGaussPoints(obj)
            coord = obj.meshUnfitted.subcellIsoCoords;
            inter = obj.unfittedInterp;
            quad  = obj.quadrature;
            quadU = obj.computePosGP(coord,inter,quad);
            obj.unfittedQuad = quadU;
        end           
        
        function shapes = computeShapes(obj)
            xGauss = obj.unfittedQuad(:,:,:);
            obj.backgroundInterp.computeShapeDeriv(xGauss);
            shapes = obj.backgroundInterp.shape;                        
        end
        
        function computeThisQuadrature(obj)
            type = obj.meshUnfitted.geometryType;
            obj.quadrature = obj.computeQuadrature(type);
        end
        
        function createBackgroundInterpolation(obj)
            mesh = obj.meshBackground;
            int = Interpolation.create(mesh,'LINEAR');
            obj.backgroundInterp = int;
        end
        
        function createUnfittedInterpolation(obj)
            mesh = obj.meshUnfitted;
            int = Interpolation.create(mesh,'LINEAR');
            obj.unfittedInterp = int;
        end
        
        function dvolu = computeJacobians(obj)
           s.mesh = obj.meshUnfitted;
           g = Geometry.create(s);
           g.computeGeometry(obj.quadrature,obj.unfittedInterp)
           dvolu = g.dvolu;
        end
        
        function Fgaus = interpolateFunctionInGaussPoints(obj,F,shapes)
            connec = obj.meshUnfitted.globalConnec;
            nCell = obj.meshUnfitted.nelem;
            nnode = obj.meshBackground.nnode;         
            ngaus = size(shapes,2);
            Fgaus = zeros(nCell,ngaus);
          %  shapes2 = permute(shapes,[2 3 1]);
             for inode = 1:nnode
                nodes  = connec(:,inode);
                Fnodes = F(nodes,1);
                for igaus = 1:ngaus
                   shape = squeeze(shapes(inode,igaus,:)); 
                   Fgaus(:,igaus) = Fgaus(:,igaus) + (shape.*Fnodes);
                end
               % shape2 = shapes2(:,:,inode);
              %  Fgaus2 = bsxfun(@times,shape2,Fnodes);
              %  Fgaus2 = sum(Fgaus2,2);
             end
        end
        
        
       
    end
    
end

