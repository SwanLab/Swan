classdef VolumeMesh < Mesh
    
    properties (Access = public)
        geometryType = 'Volume';
    end
    
    properties (Access = private)
        cParams;
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = VolumeMesh(cParams)
            obj = obj@Mesh(cParams);
            obj.initVol(cParams)
        end

        function detJ = computeJacobianDeterminant(obj,xV)
            J = obj.computeJacobian(xV);
            detJ = MatrixVectorizedInverter.computeDeterminant(J);
        end

        function invJ = computeInverseJacobian(obj,xV)
            J = obj.computeJacobian(xV);
            invJ = MatrixVectorizedInverter.computeInverse(J);
        end


        function plot(obj)
            gPar.type         = 'Full';
            g                 = GeometricalFunction(gPar);
            phiFun            = g.computeLevelSetFunction(obj);
            lsCircle          = phiFun.fValues;
            
            sUm.backgroundMesh = obj;
            sUm.boundaryMesh   = obj.createBoundaryMesh;
            uMesh              = UnfittedMesh(sUm);
            uMesh.compute(lsCircle);
            uMesh.plot
        end
        
    end
    
    methods (Access = private)
        
        function initVol(obj,cParams)
        end
        
    end
    
end