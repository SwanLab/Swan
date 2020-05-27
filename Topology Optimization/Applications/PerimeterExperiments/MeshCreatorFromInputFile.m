classdef MeshCreatorFromInputFile < handle
    
    properties (Access = public)
        mesh
        domainLength
    end
    
    properties (Access = private)
        inputFile
    end
    
    methods (Access = public)
        
        function obj = MeshCreatorFromInputFile(cParams)
            obj.init(cParams);
        end
        
        function createMesh(obj)
           [connec,coord,borderNodes,borderElements] = obj.loadSquareMeshParams();
            s.coord  = coord;
            s.connec = connec;
            s.borderNodes    = borderNodes;
            s.borderElements = borderElements;
            obj.mesh = Mesh_Total(s);            
        end        
    
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.inputFile = cParams.inputFile; 
        end
        
        function [connec,coord,borderNodes,borderElements] = loadSquareMeshParams(obj)
            eval(obj.inputFile);
            coord  = coord(:,2:3);
            connec = connec(:,2:end);
            borderNodes = [];
            borderElements = [];
            if exist('External_border_nodes')
                borderNodes    = External_border_nodes;
                borderElements = External_border_elements;
            end
            obj.domainLength = max(coord(:,1)) - min(coord(:,1));
        end                
        
    end
    
    
end