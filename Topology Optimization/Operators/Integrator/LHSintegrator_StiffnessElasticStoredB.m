classdef LHSintegrator_StiffnessElasticStoredB < LHSintegrator

    properties (Access = private)
        Btot
        material
        geometry
    end

    methods (Access = public)
        
        function obj = LHSintegrator_StiffnessElasticStoredB(cParams)
            obj.init(cParams);
            obj.material = cParams.material;
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
            obj.computeB();
        end

        function LHS = compute(obj)
            CmatTot = obj.assemblyCmat();
            LHS = obj.computeStiffness(CmatTot);
        end
        
        function setMaterialC(obj, Cmat)
            obj.material.C = Cmat;
        end
    end
    
   methods (Access = protected)
        
        function lhs = computeElementalLHS(obj)
        end

        function lhs = assembleMatrix(obj)
        end
        
   end
    
   methods (Access = private)

       function createGeometry(obj)
           s.mesh = obj.mesh;
           g      = Geometry.create(s);
           quad   = obj.quadrature;
           interp = obj.interpolation;
           g.computeGeometry(quad,interp);
           obj.geometry = g;
       end

       function computeB(obj)
           s.dim          = obj.dim;
           s.geometry     = obj.geometry;
           s.globalConnec = obj.globalConnec;
           BMC  = BMatrixComputer(s);
           obj.Btot = BMC.compute();
       end

       function CmatTot = assemblyCmat(obj)
           Cmat = obj.material.C;
           dvol = obj.geometry.dvolu;
           s.dim = obj.dim;
           s.globalConnec = [];
           s.nnodeEl = [];
           assembler = Assembler(s);
           CmatTot = assembler.assembleC(Cmat, dvol);
       end

       function K = computeStiffness(obj,CmatTot)
           B  = obj.Btot;
           CB = CmatTot*B;
           K  = B'*CB;
       end

   end
    
end