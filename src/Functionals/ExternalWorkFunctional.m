classdef ExternalWorkFunctional < handle
    
    properties (Access = private)
        mesh
    end

    properties (Access = private)
        bMesh
        testU
        bFunfExt
        bFunU
    end
    
    methods (Access = public)
        
        function obj = ExternalWorkFunctional(cParams)
            obj.init(cParams)
        end
        
        function F = computeCost(obj,u,fExt,quadOrder)
            int = Integrator.create('Function',obj.bMesh.mesh,quadOrder);
            obj.computeFunsInBoundary(u,fExt);
            F = int.compute(DP(obj.bFunU,obj.bFunfExt));
        end
        
        function Ju = computeGradient(obj,u,fExt,quadOrder)
            obj.computeFunsInBoundary(u,fExt);
            m = obj.bMesh.mesh;
            Ju = IntegrateRHS(@(v) DP(v,obj.bFunfExt),obj.testU,m,'Domain',quadOrder);
            Ju = obj.reducedToFull(Ju);
        end
        
    end
    
    methods (Access = private)

        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.bMesh = obj.mesh.createBoundaryMesh{4}; %% Accounting for loads of Top boundary {4} %%

            u = cParams.testSpace.u;
            obj.testU    = LagrangianFunction.create(obj.bMesh.mesh,u.ndimf,u.order);      
            obj.bFunfExt = LagrangianFunction.create(obj.bMesh.mesh,u.ndimf,'P1');
            obj.bFunU    = LagrangianFunction.create(obj.bMesh.mesh,u.ndimf,u.order);  
        end

        function computeFunsInBoundary(obj,u,fExt) 
            nodes = unique(obj.bMesh.globalConnec);
            if isempty(fExt)
                obj.bFunU.setFValues(u.fValues(nodes,:));
            else
                obj.bFunfExt.setFValues(fExt.fValues(nodes,:));
                obj.bFunU.setFValues(u.fValues(nodes,:));
            end
        end

        function JuFull = reducedToFull(obj,Ju)
            nNodes = obj.mesh.nnodes;
            nDim = obj.mesh.ndim;
            nDofs = nNodes*obj.mesh.ndim;
            JuFull = zeros(nDofs,1);
            nodes = unique(obj.bMesh.globalConnec);
            nNodesB = length(nodes);
            ForceDofs = zeros(nNodesB*nDim,1);
            for iDim = 1:nDim
                ForceDofs(iDim:nDim:(end-nDim+iDim)) = nDim*(nodes-1)+iDim;
            end
            JuFull(ForceDofs) = Ju;
        end

    end
    
end