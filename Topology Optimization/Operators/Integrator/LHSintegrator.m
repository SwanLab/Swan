classdef LHSintegrator < handle

    properties (Access = protected)
        mesh
        quadrature
        interpolation
        dim
        globalConnec
        dofsInElem
        material
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
            obj.dofsInElem   = cParams.dofsInElem;
        end
        
       function createQuadrature(obj)
           quad = Quadrature.set(obj.mesh.type);
           quad.computeQuadrature('LINEAR');
           obj.quadrature = quad;
       end

        function createInterpolation(obj)
            int = Interpolation.create(obj.mesh,'LINEAR');
            int.computeShapeDeriv(obj.quadrature.posgp);
            obj.interpolation = int;
        end

        function LHS = assembleMatrix(obj,LHSelem)
            s.dim          = obj.dim;
            s.globalConnec = obj.globalConnec;
            s.dofsInElem   = obj.dofsInElem;
            assembler = Assembler(s);
            LHS = assembler.assemble(LHSelem);
        end

    end
    
end