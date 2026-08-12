classdef BoundaryCoordinatesCalculator < handle
    
    properties (Access = private)
        latticeVectors
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
            
            nVectors = size(obj.latticeVectors, 1);
            
            
            if nVectors == 2
                obj.computeBoundaryParallelogram();
            else
                obj.obtainMasterBound();
                obj.obtainSlaveBound();
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.latticeVectors = cParams.latticeVectors;
            obj.nodes = cParams.nodes;
            obj.vertCoord = cParams.vertCoord;
            obj.div = cParams.div;
            obj.boundCoord = zeros(cParams.nodes.bound,2);
        end
        
        function initVertex(obj)
            nsides = obj.nodes.vert;
            obj.boundCoord(1:nsides,:) = obj.boundCoord(1:nsides,:) + obj.vertCoord;
        end
        
        function computeBoundaryParallelogram(obj)
            v1 = obj.latticeVectors(1,:);
            v2 = obj.latticeVectors(2,:);
            P1 = obj.vertCoord(1,:);  
            P2 = obj.vertCoord(2,:);  
            P3 = obj.vertCoord(3,:);  
            P4 = obj.vertCoord(4,:); 
            cont = obj.nodes.vert + 1;
            nDiv1 = obj.div(1);
            for iDiv = 1:nDiv1-1
                t = iDiv / nDiv1;
                obj.boundCoord(cont,:) = P1 + t * v1;
                cont = cont + 1;
            end
            nDiv2 = obj.div(2);
            for iDiv = 1:nDiv2-1
                t = iDiv / nDiv2;
                obj.boundCoord(cont,:) = P2 + t * v2;
                cont = cont + 1;
            end
            
            obj.MStransition = cont;
            
           
            for iDiv = 1:nDiv1-1
                t = iDiv / nDiv1;
                obj.boundCoord(cont,:) = P3 - t * v1;
                cont = cont + 1;
            end
            
            
            for iDiv = 1:nDiv2-1
                t = iDiv / nDiv2;
                obj.boundCoord(cont,:) = P4 - t * v2;
                cont = cont + 1;
            end
        end
        
        function obtainMasterBound(obj)
            nsides = obj.nodes.vert;
            cont = obj.nodes.vert + 1;
            
            for iMaster = 1:nsides/2
                v = obj.latticeVectors(iMaster, :);
                c0 = obj.vertCoord(iMaster, :);
                
                for iDiv = 1:obj.div(iMaster)-1
                    t = iDiv / obj.div(iMaster);
                    pos = c0 + t * v;
                    obj.boundCoord(cont,:) = obj.boundCoord(cont,:) + pos;
                    cont = cont + 1;
                end
            end
            obj.MStransition = cont;
        end
        
        function obtainSlaveBound(obj)
            nsides = obj.nodes.vert;
            cont = obj.MStransition;
            
            for iSlave = 1:nsides/2
                v = -obj.latticeVectors(iSlave, :);
                c0 = obj.vertCoord(nsides/2 + iSlave, :);
                
                for iDiv = 1:obj.div(iSlave)-1
                    t = iDiv / obj.div(iSlave);
                    pos = c0 + t * v;
                    obj.boundCoord(cont,:) = obj.boundCoord(cont,:) + pos;
                    cont = cont + 1;
                end
            end
        end
        
    end

end