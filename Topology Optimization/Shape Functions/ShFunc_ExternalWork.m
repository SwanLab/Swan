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
            int = Integrator.create('ScalarProduct',bMeshUp.mesh,quadOrder);
            
            [u,fExt] = obj.adaptFuns(u,fExt);
            F = int.compute(u,fExt);
        end
        
        function Ju = computeGradient(obj,u,fExt,quadOrder)
            bMesh = obj.mesh.createBoundaryMesh{4}.mesh;
            s.mesh = bMesh;
            s.quadType = quadOrder;
            s.type = 'ShapeFunction';
            RHS = RHSintegrator.create(s);

            [u,fExt] = adaptFuns(obj,u,fExt);
            test = LagrangianFunction.create(bMesh,u.ndimf,u.order);
            Ju = RHS.compute(fExt,test);
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
                uFun.fValues = u.fValues(nodes,:);
                fExtFun = LagrangianFunction.create(bMesh.mesh,u.ndimf,'P1');
            else
                fExtFun = LagrangianFunction.create(bMesh.mesh,u.ndimf,'P1');
                fExtFun.fValues = fExt.fValues(nodes,:);
                uFun = LagrangianFunction.create(bMesh.mesh,u.ndimf,'P1');
                uFun.fValues = u.fValues(nodes,:);
            end
        end
    end
    
end