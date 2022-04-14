classdef BoundaryCoordinatesCalculator < handle
    
    properties (Access = private)
        c
        theta 
        nodes
        vertCoord
        div
        MStransition
    end
    
    properties (Access = public)
        boundCoord
    end
    
    methods (Access = public)
        
        function obj = BoundaryCoordinatesCalculator(cParams)
            obj.init(cParams);
            obj.computeBoundary();
        end
        
        function computeBoundary(obj)
            obj.initVertex();
            obj.obtainMasterBound();
            obj.obtainSlaveBound();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.c = cParams.c;
            obj.theta = cParams.theta;
            obj.nodes = cParams.nodes;
            obj.vertCoord = cParams.vertCoord;
            obj.div = cParams.div;
            obj.boundCoord = zeros(cParams.nodes.bound,2);
        end
        
        function initVertex(obj)
            nsides = obj.nodes.vert;
            obj.boundCoord(1:nsides,:) = obj.boundCoord(1:nsides,:)+obj.vertCoord;
        end
        
        function obtainMasterBound(obj)
            nsides = obj.nodes.vert;
            cont = obj.nodes.vert+1;
            for iMaster = 1:nsides/2
                c0 = obj.vertCoord(iMaster,:);
                for iDiv = 1:obj.div(iMaster)-1
                    leng = obj.c(iMaster)/obj.div(iMaster);
                    angle = obj.theta(iMaster);
                    pos = NodeCoordinatesComputer.computeThePosition(c0,leng,angle);
                    obj.boundCoord(cont,:) = obj.boundCoord(cont,:)+pos;
                    cont = cont+1;
                    c0 = pos;
                end
            end
            obj.MStransition = cont;
        end
        
        function obtainSlaveBound(obj)
            nsides = obj.nodes.vert;
            cont = obj.MStransition;
            for iSlave = 1:nsides/2
                c0 = obj.vertCoord(nsides/2+iSlave,:);
                for iDiv = 1:obj.div(iSlave)-1
                    leng = obj.c(iSlave)/obj.div(iSlave);
                    angle = obj.theta(iSlave)+180;
                    pos = NodeCoordinatesComputer.computeThePosition(c0,leng,angle);
                    obj.boundCoord(cont,:) = obj.boundCoord(cont,:)+pos;
                    cont = cont+1;
                    c0 = pos;
                end
            end
        end
        
    end

end