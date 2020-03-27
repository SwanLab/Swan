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
            for isubCase = 1:obj.nSubCases
                isSubCaseActive = imax == isubCase;
                nodesSubCase = obj.nodesSubCases(:,:,isubCase,isSubCaseActive);
                nodesQ(:,:,isSubCaseActive) = nodesSubCase;
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
                nodesT = obj.nodesSubCases(:,:,isubCase,:);
                q = obj.computeQuality(nodesT,obj.nElemInCase);
                qT(isubCase,:) = sum(q,1);
            end
            [~,imax] = max(qT);
       end        
        
        function q = computeQuality(obj,allConnec,nElemInCase)
            coordT = obj.coord;
            A = zeros(2,nElemInCase);
            q = zeros(2,nElemInCase);
            for isubElem = 1:2
                
                nodes = squeeze(allConnec(isubElem,:,:))';
                
                xA = coordT(nodes(:,1),1);
                yA = coordT(nodes(:,1),2);
                
                xB = coordT(nodes(:,2),1);
                yB = coordT(nodes(:,2),2);
                
                xC = coordT(nodes(:,3),1);
                yC = coordT(nodes(:,3),2);
                
                A(isubElem,:) =1/2*abs((xA -xC).*(yB-yA)-(xA-xB).*(yC-yA));
                Lab = (xA - xB).^2 + (yA - yB).^2;
                Lcb = (xC - xB).^2 + (yC - yB).^2;
                Lac = (xA - xC).^2 + (yA - yC).^2;
                L = Lab + Lcb + Lac;
                q(isubElem,:) = 4*sqrt(3)*A(isubElem,:)./L';
            end
        end       
       
   end
    
   
    
    
end