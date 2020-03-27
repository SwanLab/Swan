classdef BestSubCellCaseSelector < handle
    
    properties (Access = private)
        coord
        nodesSubCases
        nSubCases
        nSubCellNodes
        nElemInCase
        nSubCellsByQuad
    end
    
    methods (Access = public)
        
        function obj = BestSubCellCaseSelector(cParams)
            obj.init(cParams)
        end
        
        function nodesQ = compute(obj)
            nodesQ = obj.initNodeQ();
            imax = obj.computeBetterSubCaseOption();
            for iCase = 1:obj.nSubCases
                isActive = imax == iCase;
                nodesS = obj.nodesSubCases(:,:,iCase,isActive);
                nodesQ(:,:,isActive) = nodesS;
            end
        end
        
    end
        
    methods (Access = private)
        
        function init(obj,cParams)
            obj.coord         = cParams.coord;
            obj.nodesSubCases = cParams.nodesSubCases;
            obj.nSubCellsByQuad = size(obj.nodesSubCases,1);
            obj.nSubCellNodes   = size(obj.nodesSubCases,2);
            obj.nSubCases       = size(obj.nodesSubCases,3);
            obj.nElemInCase     = size(obj.nodesSubCases,4);
        end
        
        function nodes = initNodeQ(obj)
            nodes = zeros(obj.nSubCellsByQuad,obj.nSubCellNodes,obj.nElemInCase);
        end
        
        function imax = computeBetterSubCaseOption(obj)
            qT = zeros(obj.nSubCases,obj.nElemInCase);
            for isubCase = 1:obj.nSubCases
                nodesS = obj.nodesSubCases(:,:,isubCase,:);
                q = obj.computeCaseQuality(nodesS);
                qT(isubCase,:) = sum(q,1);
            end
            [~,imax] = max(qT);
        end
        
        function q = computeCaseQuality(obj,nodeSubCases)
            q = zeros(obj.nSubCases,obj.nElemInCase);
            nodeSubCases = permute(nodeSubCases,[4 2 1 3]);
            for iCase = 1:obj.nSubCases
                nodes = nodeSubCases(:,:,iCase);
                q(iCase,:) = obj.computeQuality(nodes);
            end
        end
        
        function q = computeQuality(obj,nodes)
            s.connec = nodes;
            s.coord = obj.coord;
            m = Mesh().create(s);            
            q = m.computeElementQuality();
        end
      
        
    end
    
end