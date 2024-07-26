classdef BestSubCellCaseSelector < handle
    
   properties (Access = private)
        coord
        nodesSubCases
        nSubCases
        nElemInCase
    end
    
    methods (Access = public)
        
        function obj = BestSubCellCaseSelector(cParams)
            obj.init(cParams)
        end
        
        function imax = compute(obj)
            qT = zeros(obj.nSubCases,obj.nElemInCase);
            for isubCase = 1:obj.nSubCases
                nodesS = obj.nodesSubCases(:,:,isubCase,:);
                q = obj.computeCaseQuality(nodesS);
                qT(isubCase,:) = obj.computeMeritFunction(q);
            end
            [~,imax] = max(qT);
        end
        
    end
        
    methods (Access = private)
        
        function init(obj,cParams)
            obj.coord         = cParams.coord;
            obj.nodesSubCases = cParams.nodesSubCases;
            obj.nSubCases       = size(obj.nodesSubCases,3);
            obj.nElemInCase     = size(obj.nodesSubCases,4);
        end
   
        function q = computeCaseQuality(obj,nodeSubCases)
            q = zeros(obj.nSubCases,obj.nElemInCase);
            nodeSubCases = permute(nodeSubCases,[4 2 1 3]);
            for iCase = 1:obj.nSubCases
                nodes = nodeSubCases(:,:,iCase);
                q(iCase,:) = obj.computeQuality(nodes);
            end
        end
        
        function m = computeMeritFunction(obj,q)
            quality     = sum(q,1);
            meanQuality = quality/obj.nSubCases;
            desviation = zeros(size(quality));
            for iCase = 1:obj.nSubCases
                desv = (q(iCase,:) - meanQuality).^2;
                desviation = desviation + desv;
            end                        
            homogeonity = -sqrt(desviation);
            m = 0.9*quality + 0.1*homogeonity;
        end
        
        function q = computeQuality(obj,nodes)
            s.connec = nodes;
            s.coord = obj.coord;
            if size(obj.coord,2) == 3
               s.kFace = -1;
            end
            m = Mesh.create(s);
            q = m.computeElementQuality();
        end
      
        
    end
    
end