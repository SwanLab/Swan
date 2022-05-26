classdef LHSintegrator < handle

    properties (Access = protected)
        mesh
        quadrature
        interpolation
        dim
        globalConnec
    end

    methods (Access = public, Static)
        
        function obj = create(s)
            f = LHSintegratorFactory();
            obj = f.create(s);
        end

    end

    methods (Access = public)
        
        function q = getQuadrature(obj)
            q = obj.quadrature;
        end
 
    end

    methods (Access = protected)
     
        function init(obj,cParams)
            obj.dim          = cParams.dim;
            obj.mesh         = cParams.mesh;
            obj.globalConnec = cParams.globalConnec;
        end
        
       function createQuadrature(obj)
           quad = Quadrature.set(obj.mesh.type);
           quad.computeQuadrature('LINEAR'); % QUADRATIC LINEAR
           obj.quadrature = quad;
       end

        function createInterpolation(obj)
            int = obj.mesh.interpolation;
            int.computeShapeDeriv(obj.quadrature.posgp);
            obj.interpolation = int;
        end

        function LHS = assembleMatrix(obj,LHSelem)
            s.dim          = obj.dim;
            s.globalConnec = obj.globalConnec;
            s.nnodeEl      = obj.interpolation.nnode;
            assembler = Assembler(s);
            LHS = assembler.assemble(LHSelem);
        end

    end
    
end