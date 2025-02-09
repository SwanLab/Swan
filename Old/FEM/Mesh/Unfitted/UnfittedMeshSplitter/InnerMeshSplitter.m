classdef InnerMeshSplitter < SubUnfittedMeshSplitter
       
    methods (Access = public)
        
        function obj = InnerMeshSplitter(cParams)
            obj.init(cParams)
        end
        
    end
    
    methods (Access = protected, Static)
        
        function sM = createSubMeshComponent(m,nodes)
            s.fullCells      = m.fullCells(nodes);
            s.backgroundMesh = m.backgroundMesh;
            sM = InnerMesh(s);            
        end
        
    end
    
end