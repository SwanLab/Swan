classdef VademecumDifferenceMakerInTxiRhoPlane < handle
    
    properties (Access = public)
       vademecumDiff         
    end
    
    properties (Access = private)
        vademecums
        communNodes
        index
        iIndex
        jIndex
        txi1
        txi2
        rho1
        rho2
        tensor1
        tensor2
        tensor
        tensorDiff
    end
    
    methods (Access = public)
        
        function obj = VademecumDifferenceMakerInTxiRhoPlane(cParams)
            obj.init(cParams)
        end
        
        function compute(obj)
            obj.computeTxiVariable();
            obj.computeCommunNodes();
            obj.computeTxiRhoInCommunNodes();
            obj.computeVademecumDiff();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.index = [1 1; 2 2;3 3; 1 2;2 3;1 3];
            obj.vademecums{1} = cParams.firstVademecum;
            obj.vademecums{2} = cParams.secondVademecum;
        end
        
        function computeCommunNodes(obj)
            nodes1 = obj.vademecums{1}.feasibleIndex;
            nodes2 = obj.vademecums{2}.feasibleIndex;
            obj.communNodes = intersect(nodes1,nodes2);
        end
        
        function computeTxiRhoInCommunNodes(obj)
            obj.txi1 = obj.vademecums{1}.chi(obj.communNodes);
            obj.rho1 = obj.vademecums{1}.volume(obj.communNodes);
            obj.txi2 = obj.vademecums{2}.chi(obj.communNodes);
            obj.rho2 = obj.vademecums{2}.volume(obj.communNodes);                        
        end
        
        function computeVademecumDiff(obj)
            obj.vademecumDiff.mxV  =  obj.vademecums{1}.mxV;
            obj.vademecumDiff.myV  =  obj.vademecums{1}.myV;
            obj.vademecumDiff.C    = obj.computeHomogenizedTensorDiff();
            obj.vademecumDiff.invP = obj.computeAmplificatorTensorDiff(); 
            obj.vademecumDiff.volume = obj.vademecums{2}.volume;            
        end
        
        function computeTxiVariable(obj)
            for ivad = 1:numel(obj.vademecums)
                mxV = obj.vademecums{ivad}.mxV;
                myV = obj.vademecums{ivad}.myV;
                for i = 1:length(mxV)
                    for j = 1:length(myV)
                        mx = mxV(i);
                        my = myV(j);
                        chi = atan(mx/my);
                        obj.vademecums{ivad}.chi(i,j) = chi;
                    end
                end
            end
        end
        
        function Cdif = computeHomogenizedTensorDiff(obj)
            obj.tensor1 = obj.vademecums{1}.C;
            obj.tensor2 = obj.vademecums{2}.C;
            Cdif = obj.computeDifference();
        end        
        
        function invPdif = computeAmplificatorTensorDiff(obj)
            obj.tensor1 = obj.vademecums{1}.invP;
            obj.tensor2 = obj.vademecums{2}.invP;
            invPdif = obj.computeDifference();
        end
        
        function tDif = computeDifference(obj)
            nComponent = length(obj.index);
            for icomponent = 1:nComponent
                obj.obtainIJindex(icomponent);
                obj.computeDiff();
            end
            tDif = obj.tensorDiff;
        end
        
        function obtainIJindex(obj,iplot)
            obj.iIndex = obj.index(iplot,1);
            obj.jIndex = obj.index(iplot,2);
        end        
        
        function zDif = computeDiff(obj)
            i = obj.iIndex;
            j = obj.jIndex;  
            t1 = obj.tensor1(i,j,:,:);
            t2 = obj.tensor2(i,j,:,:);            
            z1 = t1(obj.communNodes);
            z2 = t2(obj.communNodes);            
            z1Int = griddata(obj.txi2,obj.rho2,z2,obj.txi1,obj.rho1);
            zDif = (z1 - z1Int);
            obj.tensorDiff(i,j,:,:) = zDif;            
        end
        
    end
    
end