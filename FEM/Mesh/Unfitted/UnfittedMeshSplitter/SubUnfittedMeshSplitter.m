classdef SubUnfittedMeshSplitter < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        subMesh
    end
       
    methods (Access = public, Static)
        
        function obj = create(cParams)            
            f = SubUnfittedMeshSplitterFactory;
            obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
        
        function mI = split(obj)
            m         = obj.subMesh;
            faces     = m.mesh.connec;
            compLabel = obj.computeComponentLabel(faces);
            nComp = unique(compLabel);
            mI    = cell(length(nComp),1);
            for iComp = 1:length(nComp)
                isComp = compLabel == nComp(iComp);
                sM = obj.createSubMeshComponent(m,isComp);
                mI{iComp} = sM;
            end
        end            
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.subMesh = cParams.subMesh;
        end    
           
    end
    
    methods (Access = private, Static)
        
       function compID = computeComponentLabel(faces)
            s.faces = faces;
            sp = SplitterInConnectedComponents(s);
            mComp = sp.split();
            compID = mComp(faces(:,1))';
        end                
        
    end
    
    methods (Access = protected, Static)
        createSubMeshComponent(obj)        
    end
    
end