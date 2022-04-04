classdef NewFilter < handle
    
    properties (GetAccess = public, SetAccess = protected)
        diffReacProb
        ngaus
        nelem
    end
    
    properties (Access = protected)
        x
        x_reg
    end
    
    properties (GetAccess = protected, SetAccess = protected)
        Poper  
        geometry
        quadrature   
        mesh
        nnode
        npnod
        shape  
        quadratureOrder
    end
    
    properties (Access = protected)
        interp
        M
        Kernel
    end

    methods (Access = public, Static)

        function obj = create(cParams)
            f = FilterFactory();
            obj = f.create(cParams);
        end

    end

    methods(Access = public)

        function varInitialization(obj,cParams)
            obj.init(cParams);
            obj.createMassMatrix(cParams);
            obj.preProcess();
            obj.createPoperator(cParams);
        end

        function preProcess(obj)
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
            obj.storeParams();
        end
        
    end

    methods (Access = public, Static)

        function pB = createDiffReacProblem(cParams)
            s = cParams.femSettings;
            s.mesh = cParams.mesh;
            switch s.scale
                case 'MACRO'
                    pB = NewDiffReactProblem(s);
                case 'MICRO'
                    pB = NewDiffReactProblemMicro(s);
            end
        end

    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.createDiffReacProblem(cParams);
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end
        
        function itHas = xHasChanged(obj,x)
            itHas = ~isequal(x,obj.x);
        end
        
        function updateStoredValues(obj,x,x0)
            obj.x = x;
            obj.x_reg = x0;
        end
        
    end
    
    methods (Access = private)
        
        function createInterpolation(obj)
            obj.interp = Interpolation.create(obj.mesh,'LINEAR');
        end
        
        function createGeometry(obj)
            s.mesh = obj.mesh;
            obj.geometry = Geometry.create(s);
            obj.geometry.computeGeometry(obj.quadrature,obj.interp);
        end
        
        function createQuadrature(obj)
            obj.quadrature = Quadrature.set(obj.mesh.type);
            obj.quadrature.computeQuadrature(obj.quadratureOrder);
        end
        
        function storeParams(obj)
            obj.nelem = obj.mesh.nelem;
            obj.nnode = obj.mesh.nnode;
            obj.npnod = obj.mesh.npnod;
            obj.ngaus = obj.quadrature.ngaus;
            obj.shape = obj.interp.shape;
        end

        function createMassMatrix(obj,cParams)
            obj.diffReacProb = obj.createDiffReacProblem(cParams);
            obj.M = obj.diffReacProb.getM();
        end

        function createPoperator(obj,cPar)
            cParams.nelem  = obj.mesh.nelem;
            cParams.nnode  = obj.mesh.nnode;
            cParams.npnod  = obj.mesh.npnod;
            cParams.connec = obj.mesh.connec;
            cParams.diffReactEq = cPar.femSettings;
            obj.Poper = Poperator(cParams);
        end

    end

end