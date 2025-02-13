classdef ShFunc_ExternalWork2 < handle

   properties (Access = private)
       mesh
       quadOrder
   end

   properties (Access = private)
       bMesh
       test
       RHS
   end
     
   methods (Access = public)
      
       function obj = ShFunc_ExternalWork2(cParams)
           obj.init(cParams);
           obj.createRHSintegrator();
       end

       function setTestFunction(obj,u)           
            obj.test = LagrangianFunction.create(obj.bMesh.mesh, u.ndimf, u.order);           
        end 
      
        function F = computeFunction(obj,u,fExt)
          bMesh = obj.mesh.createBoundaryMesh{4}; %CARA SUPERIOR
          int = Integrator.create('Function',bMesh.mesh,obj.quadOrder);
          [u, fExt] = obj.adaptFuns(u,fExt);
          F = int.compute(u.*fExt);
       end
      
       function Ju = computeResidual(obj,u,fExt)
            [~,fExt] = adaptFuns(obj,u,fExt);
            Ju = obj.RHS.compute(fExt,obj.test);
            Ju = obj.reducedToFull(Ju,obj.bMesh);
       end
      
   end
  
   methods (Access = private)
      
       function init(obj,cParams)
         obj.mesh = cParams.mesh;
         obj.quadOrder = cParams.quadOrder;
         obj.bMesh = obj.mesh.createBoundaryMesh{4}; %CARA SUPERIOR
       end

       function createRHSintegrator(obj)
            s.mesh = obj.bMesh.mesh;
            s.quadType = obj.quadOrder;
            s.type = 'ShapeFunction';
            obj.RHS = RHSIntegrator.create(s);
       end

        function [uFun,fExtFun] = adaptFuns(obj,u,fExt) %% ADAPTING TO TOP BOUNDARY %%
            nodes = unique(obj.bMesh.globalConnec);
            if isempty(fExt)
                uFun = LagrangianFunction.create(obj.bMesh.mesh,u.ndimf,'P1');
                uFun.setFValues (u.fValues(nodes,:));
                fExtFun = LagrangianFunction.create(obj.bMesh.mesh,u.ndimf,'P1');
            else
                fExtFun = LagrangianFunction.create(obj.bMesh.mesh,u.ndimf,'P1');
                fExtFun.fValues = fExt.fValues(nodes,:);
                uFun = LagrangianFunction.create(obj.bMesh.mesh,u.ndimf,'P1');
                uFun.fValues = u.fValues(nodes,:);
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

