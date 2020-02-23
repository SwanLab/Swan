classdef IntegratorCutMesh < Integrator
    
    properties (GetAccess = public, SetAccess = protected)
        cutMesh
    end
    
    properties (Access = private)
        Fnodal
        Fprojection
        unfittedInterp
        quadrature
        shapes
        RHScells
        RHScellsCut
        backgroundMesh
    end
    
    methods (Access = public)
        
        function obj = IntegratorCutMesh(cParams)
            obj.init(cParams);
            obj.backgroundMesh = cParams.meshBackground;
            obj.cutMesh        = obj.mesh;
        end
        
        function A = integrate(obj,Fnodal,quad)
           obj.Fnodal = Fnodal;
           if nargin == 3
               obj.quadrature = quad;
           else
               type = obj.mesh.geometryType;
               obj.quadrature = obj.computeQuadrature(type);
           end
            obj.computeElementalRHS();
            obj.assembleSubcellsInCells();
            A = obj.assembleIntegrand();
        end
        
    end
    
    methods (Access = private)
        
        
        function computeElementalRHS(obj)
            obj.createUnfittedInterpolation();
            obj.computeShapesInUnfittedGaussPoints();
            int = obj.integrateFwithShapeFunction();
            obj.RHScellsCut = int;            
        end
        
        function assembleSubcellsInCells(obj)
            nnode = obj.backgroundMesh.nnode;
            nelem = obj.backgroundMesh.nelem;
            cellNum = obj.cutMesh.cellContainingSubcell;
            totalInt = zeros(nelem,nnode);
            for iNode = 1:nnode
                int = obj.RHScellsCut(:,iNode);
                intGlobal  = accumarray(cellNum,int,[nelem,1],@sum,0);
                totalInt(:,iNode) = totalInt(:,iNode) + intGlobal;
            end
            obj.RHScells = totalInt;
        end
        
        function f = assembleIntegrand(obj)
            integrand = obj.RHScells;
            npnod  = obj.backgroundMesh.npnod;
            nnode  = obj.backgroundMesh.nnode;
            connec = obj.backgroundMesh.connec;
            f = zeros(npnod,1);
            for inode = 1:nnode
                int = integrand(:,inode);
                con = connec(:,inode);
                f = f + accumarray(con,int,[npnod,1],@sum,0);
            end
        end
        
        function itIs = isLeveSetCuttingMesh(obj)
            itIs = ~isempty(obj.cutMesh.backgroundCutCells);
        end
        
    end
    
    methods (Access = private)
        
        function computeShapesInUnfittedGaussPoints(obj)
            xG  = obj.computeUnfittedGaussPoints();
            int = obj.createBackgroundInterpolation();
            int.computeShapeDeriv(xG);
            obj.shapes = int.shape;
        end
        
        function int = createBackgroundInterpolation(obj)
            m = obj.backgroundMesh;
            int = Interpolation.create(m,'LINEAR');
        end
        
        function xGauss = computeUnfittedGaussPoints(obj)
            quad  = obj.quadrature;
            coord = obj.cutMesh.subcellIsoCoords;
            obj.unfittedInterp.computeShapeDeriv(quad.posgp);
            xGauss = zeros(size(coord,3),quad.ngaus,size(coord,1));
            shape = obj.unfittedInterp.shape;
            for igaus = 1:quad.ngaus
                for idime = 1:size(coord,3)
                    xGauss(idime,igaus,:) = coord(:,:,idime)*shape(:,igaus);
                end
            end
        end
        
        function createUnfittedInterpolation(obj)
            int = Interpolation.create(obj.mesh,'LINEAR');
            obj.unfittedInterp = int;
        end
        

        function Fproj = integrateFwithShapeFunction(obj)
            Fgauss = obj.interpolateFunctionInGaussPoints();
            dvolume = obj.computeDvolume();            
            fdV = (Fgauss.*dvolume);
            nelem = obj.mesh.nelem;
            nnode = obj.backgroundMesh.nnode;            
            Fproj = zeros(nnode,nelem);
            for igaus = 1:obj.quadrature.ngaus
               fdVig = fdV(igaus,:);
               shapeIg = squeeze(obj.shapes(:,igaus,:));            
               Fproj = Fproj + bsxfun(@times,shapeIg,fdVig);
            end
            Fproj = Fproj';
        end
                
        function Fgaus = interpolateFunctionInGaussPoints(obj)
            connec = obj.cutMesh.globalConnec;
            nCell  = obj.cutMesh.nelem;
            nnode = obj.backgroundMesh.nnode;
            ngaus = obj.quadrature.ngaus;
            Fgaus = zeros(ngaus,nCell);
             for inode = 1:nnode
                nodes  = connec(:,inode);
                Fnodes = obj.Fnodal(nodes,1);
                for igaus = 1:ngaus
                   shape = squeeze(obj.shapes(inode,igaus,:)); 
                   Fgaus(igaus,:) = Fgaus(igaus,:) + (shape.*Fnodes)';
                end
             end
        end

        function dvolume = computeDvolume(obj)
            s.mesh = obj.cutMesh;
            g = Geometry.create(s);
            g.computeGeometry(obj.quadrature,obj.unfittedInterp)
            dvolume = g.dvolu;
            dvolume = dvolume';
        end
        
    end
    
end

