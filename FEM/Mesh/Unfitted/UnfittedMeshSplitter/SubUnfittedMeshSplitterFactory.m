classdef SubUnfittedMeshSplitterFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(s)
            switch class(s.subMesh)
                case 'InnerMesh'
                    obj = InnerMeshSplitter(s);
                case 'InnerCutMesh'
                    obj = InnerCutMeshSplitter(s);
            end
            
        end
        
    end    
    
end