classdef TriangleLinear_New < handle
    
    properties (GetAccess = public, SetAccess = protected)
        order
        ndime
        nnode
        pos_nodes
    end
    
    properties (Access = private)
        shapeFun
        shapeDer
        mesh
    end
    
    methods (Access = public)
        
        function obj = TriangleLinear_New(cParams)
            obj.init();
            obj.createShapeFunctions();
            obj.createShapeDerivatives();
        end
        
        function shape = computeShapeFunctions(obj,xV)
            shape = obj.shapeFun.evaluate(xV);
        end
        
        function deriv = computeShapeDerivatives(obj,xV)
            deriv = obj.shapeFun.evaluate(xV);
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.ndime = 2;
            obj.nnode = 3;
            obj.pos_nodes = [0 0; 1 0; 0 1];
            
            s.coord = obj.pos_nodes;
            s.connec = [1 2 3];
            obj.mesh = ReferenceMesh(s); % loop...
        end
        
        function createShapeFunctions(obj)
            N1 = @(x) 1- x(1,:,:) - x(2,:,:);
            N2 = @(x) x(1,:,:);
            N3 = @(x) x(2,:,:);
            obj.shapeFun = AnalyticalFunction.create(@(x) [N1(x),N2(x),N3(x)],...
                3, obj.mesh);
        end
        
        function createShapeDerivatives(obj)
            % sth like analyticaltensor?
            dN1dXi = @(x) -1;
            dN2dXi = @(x) +1;
            dN3dXi = @(x)  0;
            dN1dEta = @(x) -1;
            dN2dEta = @(x)  0;
            dN3dEta = @(x) +1;
            obj.shapeDer = AnalyticalFunction.create(@(x) [dN1dXi(x),dN2dXi(x),dN3dXi(x);
                                                        dN1dEta(x),dN2dEta(x),dN3dEta(x);],...
                                                        [2 3], obj.mesh);
        end
        
    end
    
end
