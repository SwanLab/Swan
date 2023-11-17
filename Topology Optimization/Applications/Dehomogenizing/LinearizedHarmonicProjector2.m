classdef LinearizedHarmonicProjector2 < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        eta
        epsilon
        internalDOFs
        massMatrix
        stiffnessMatrix
    end
    
    properties (Access = private)
        mesh
        boundaryMesh
    end
    
    methods (Access = public)
        
        function obj = LinearizedHarmonicProjector2(cParams)
            obj.init(cParams)
            obj.createInternalDOFs();
            obj.computeMassMatrix();
            obj.stiffnessMatrix();
        end

        function [lRes,hRes] = evaluateAllResiduals(obj,bBar,b)
            lRes = obj.evaluateLossResidual(bBar,b);
        end

        function lRes = evaluateLossResidual(obj,bBar,b)
          difB      = abs(b.fValues - bBar.fValues);
          s.fValues = difB;
          s.mesh    = obj.mesh;
          lRes = P1Function(s);
        end

        function evaluateHarmonicResidual(obj,b)
            


        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.mesh             = cParams.mesh;
           obj.boundaryMesh     = cParams.boundaryMesh;
           obj.eta     = (obj.mesh.computeMeanCellSize)^2;                            
           obj.epsilon = 100;
        end
        
        function computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.test  = P1Function.create(obj.mesh, 1);
            s.trial = P1Function.create(obj.mesh, 1);
            s.quadratureOrder = 'QUADRATICMASS';
            lhs = LHSintegrator.create(s);
            M = lhs.compute();           
            obj.massMatrix = M;
        end        

        function computeStiffnessMatrix(obj)        
            s.test  = P1Function.create(obj.mesh, 1);
            s.trial = P1Function.create(obj.mesh, 1);
            s.mesh         = obj.mesh;
            s.type         = 'StiffnessMatrix';
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
            obj.stiffnessMatrix = K;
        end   
        
        function createInternalDOFs(obj)
           b     = obj.boundaryMesh;
           iDOFs = setdiff(1:obj.mesh.nnodes,b); 
           obj.internalDOFs = iDOFs;
        end   

        function createStiffNessWithFunction(obj,f)
            s.test   = P1Function.create(obj.mesh, 1);
            s.trial = P1Function.create(obj.mesh, 1);
            s.mesh  = obj.mesh;
            s.type  = 'StiffnessMatrix';
            lhs = LHSintegrator.create(s);
            K = lhs.compute();
            obj.stiffnessMatrix = K;
        end
        
        
    end
    
end