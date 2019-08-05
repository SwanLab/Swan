classdef Mesh_Unfitted < handle
    
    properties (GetAccess = public, SetAccess = protected)
%         meshBackground
        levelSet_background
        levelSet_unfitted
    end
    
    methods (Access = public, Abstract)
        
        computeMesh(obj,levelSet)
        computeMass(obj)
        plot(obj)
        add2plot(obj)
        
    end
    
    methods (Access = public, Static)
       
        function obj = create2(cParams)
            f   = Mesh_Unfitted_Factory();
            obj = f.create(cParams);            
        end
        
    end
    
    methods (Access = public)
        
        function setLevelSetUnfitted(obj,LS)
            obj.levelSet_unfitted = LS;
        end
        
        function setLevelSetBackground(obj,LS)
            obj.levelSet_background = LS;
        end
        
    end
    
end