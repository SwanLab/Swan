classdef EdgeFunctionInterpolator < handle
    
    properties (Access = private)
        edgeMesh
        fInEdge
    end
    
    properties (Access = private)
        mesh
        fNodes
    end
    
    methods (Access = public)
        
        function obj = EdgeFunctionInterpolator(cParams)
            obj.init(cParams)            
        end
        
        function fE = compute(obj)
            obj.createP1FunctionInEdges();
            fE = obj.interpolateInMiddleEdge();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.edgeMesh = cParams.edgeMesh;
            obj.fNodes   = cParams.fNodes;
        end
        
        function f = createP1FunctionInEdges(obj)
            s.mesh    = obj.edgeMesh;
            s.fValues = obj.fNodes;
            s.functionType = 'P1';
            s.order   = 'P1';
            f = LagrangianFunction(s);
            obj.fInEdge = f;
        end
        
        function fE = interpolateInMiddleEdge(obj)
            m = obj.edgeMesh;
            q = Quadrature.create(m,1);
            xV = q.posgp;
            f  = obj.fInEdge;
            fE = squeeze(f.evaluate(xV));            
        end
        
    end
    
end