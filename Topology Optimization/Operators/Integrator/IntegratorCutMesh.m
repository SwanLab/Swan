classdef IntegratorCutMesh < Integrator
    
    properties (Access = private)
        Fnodal
        quadrature
        RHScells
        RHScellsCut
        backgroundMesh
        xGauss
    end
    
    methods (Access = public)
        
        function obj = IntegratorCutMesh(cParams)
            obj.init(cParams);
            obj.backgroundMesh = cParams.meshBackground;
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
            obj.computeUnfittedGaussPoints();            
            obj.computeShapeFunctions();
            int = obj.integrateFwithShapeFunction();
            obj.RHScellsCut = int;
        end
        
        function computeUnfittedGaussPoints(obj)
            q = obj.quadrature;
            m = obj.mesh.cutMeshOfSubCellLocal;
            obj.xGauss = m.computeXgauss(q.posgp);
        end
        
        function int = integrateFwithShapeFunction(obj)
            Fgauss  = obj.computeFinGaussPoints();
            dV      = obj.computeDvolume;
            fdV     = (Fgauss.*dV);
            shapes  = obj.computeShapeFunctions();
            int = obj.initIntegrand();
            for igaus = 1:obj.quadrature.ngaus
                fdv   = fdV(igaus,:);
                shape = shapes(:,:,igaus);
                int = int + bsxfun(@times,shape,fdv);
            end
            int = transpose(int);
        end        
        
        function fG = computeFinGaussPoints(obj)
            f = createFeFunction(obj);
            fG = f.interpolateFunction(obj.xGauss);
            fG = permute(fG,[2 3 1]);            
        end
        
        function f = createFeFunction(obj)            
            m = obj.mesh.cutMeshOfSubCellGlobal();
            s.fNodes = obj.Fnodal;
            s.mesh = m;
            f = FeFunction(s);
        end    
        
        function dV = computeDvolume(obj)
            q  = obj.quadrature;
            dV = obj.mesh.computeDvolume(q);            
        end
        
        function shapes = computeShapeFunctions(obj)
            m = obj.backgroundMesh;
            int = Interpolation.create(m,'LINEAR');
            int.computeShapeDeriv(obj.xGauss);
            shapes = permute(int.shape,[1 3 2]);            
        end       
        
        function int = initIntegrand(obj)
            nelem = obj.mesh.nelem;
            nnode = obj.backgroundMesh.nnode;
            int = zeros(nnode,nelem);            
        end
        
        function assembleSubcellsInCells(obj)
            nnode = obj.backgroundMesh.nnode;
            nelem = obj.backgroundMesh.nelem;
            cellNum = obj.mesh.cellContainingSubcell;
            totalInt = zeros(nelem,nnode);
            for iNode = 1:nnode
                int = obj.RHScellsCut(:,iNode);
                intGlobal = accumarray(cellNum,int,[nelem,1],@sum,0);
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
        
    end

end

