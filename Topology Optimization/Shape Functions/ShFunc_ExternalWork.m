classdef ShFunc_ExternalWork < handle
    
    properties (Access = private)
        mesh
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = ShFunc_ExternalWork(cParams)
            obj.init(cParams)
        end
        
        function F = computeFunction(obj,u,fExt,quadOrder)
            bMeshUp = obj.mesh.createBoundaryMesh{4};
            int = Integrator.create('Function',bMeshUp.mesh,quadOrder);
            
            [u,fExt] = obj.adaptFuns(u,fExt);
            F = int.compute(u.*fExt);
        end
        
        function Ju = computeGradient(obj,u,fExt,quadOrder)
            bMesh = obj.mesh.createBoundaryMesh{4};
            s.mesh = bMesh.mesh;
            s.quadType = quadOrder;
            s.type = 'ShapeFunction';
            RHS = RHSintegrator.create(s);

            [u,fExt] = adaptFuns(obj,u,fExt);
            test = LagrangianFunction.create(bMesh.mesh,u.ndimf,u.order);
            Ju = RHS.compute(fExt,test);
            Ju = obj.reducedToFull(Ju,bMesh);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
        end
        
        function [uFun,fExtFun] = adaptFuns(obj,u,fExt) %% ADAPTING TO TOP BOUNDARY %%
            bMesh = obj.mesh.createBoundaryMesh{4};
            nodes = unique(bMesh.globalConnec);
            if isempty(fExt)
                uFun = LagrangianFunction.create(bMesh.mesh,u.ndimf,'P1');
                uFun.setFValues(u.fValues(nodes,:));
                fExtFun = LagrangianFunction.create(bMesh.mesh,u.ndimf,'P1');
            else
                fExtFun = LagrangianFunction.create(bMesh.mesh,u.ndimf,'P1');
                fExtFun.fValues = fExt.fValues(nodes,:);
                uFun = LagrangianFunction.create(bMesh.mesh,u.ndimf,'P1');
                uFun.setFValues(u.fValues(nodes,:));
            end
        end

        function JuFull = reducedToFull(obj,Ju,bMesh)
            nNodes = obj.mesh.nnodes;
            nDim = obj.mesh.ndim;
            nDofs = nNodes*obj.mesh.ndim;
            JuFull = zeros(nDofs,1);
            nodes = unique(bMesh.globalConnec);
            nNodesB = length(nodes);
            ForceDofs = zeros(nNodesB*nDim,1);
            for iDim = 1:nDim
                ForceDofs(iDim:nDim:(end-nDim+iDim)) = nDim*(nodes-1)+iDim;
            end
            JuFull(ForceDofs) = Ju;
        end

    end
    
end