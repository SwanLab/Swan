classdef DesignVarMonitor_LevelSet < DesignVarMonitor_Abstract
    
    properties (Access = protected, Abstract)
        unfittedType
        meshIncludeBoxContour
    end
    
    properties (Access = protected)
        designVarName = 'Level Set - \phi';
        meshUnfitted
    end
    
    methods (Access = public)
        
        function obj = DesignVarMonitor_LevelSet(mesh)
            obj@DesignVarMonitor_Abstract(mesh);
            obj.createUnfittedMesh();
        end
        
        function plot(obj,phi)
            obj.updateMesh(phi);
            obj.refreshFigure();
        end
        
    end
    
    methods (Access = protected)
        
        function initPlotting(obj)
            obj.patchHandle = patch('Faces',obj.mesh.connec,'Vertices',obj.mesh.coord,...
                'FaceColor',obj.getColor(),'EdgeColor',[0 0 0],'EdgeAlpha',0.5,...
                'FaceLighting','flat','EdgeLighting','flat');
            
            set(obj.axes,'CLim',[0, 1],'XTick',[],'YTick',[]);
        end
        
    end
    
    methods (Access = protected, Static)
        
        function color = getColor()
            color = [1 0 0];
        end
        
    end
    
    methods (Access = private)
        
        function createUnfittedMesh(obj)
            interpolation = Interpolation.create(obj.mesh,'LINEAR');
            obj.meshUnfitted = Mesh_Unfitted_Factory.create(obj.unfittedType,obj.mesh,interpolation,...
                'includeBoxContour',obj.meshIncludeBoxContour);
        end
        
        function updateMesh(obj,phi)
            obj.meshUnfitted.computeMesh(phi);
        end
        
        function refreshFigure(obj)
            cla(obj.axes);
            obj.meshUnfitted.add2plot(obj.axes);
            light(obj.axes)
            obj.BCplotter.plot();
        end
        
    end
    
end