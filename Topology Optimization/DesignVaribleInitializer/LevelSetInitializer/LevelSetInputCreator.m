classdef LevelSetInputCreator < handle
    
    properties (Access = private)
        input
    end
    
    methods (Access = public)
        
        function obj = LevelSetInputCreator(settings,mesh)
            input.nHoles        = settings.N_holes;
            input.rHoles        = settings.R_holes;
            input.phaseHoles    = settings.phase_holes;
            input.warningHoleBC = settings.warningHoleBC;
            input.m             = settings.widthSquare;
            input.m1            = settings.widthH;
            input.m2            = settings.widthV;
            input.levFib        = settings.levFib;
            input.yn            = settings.yn;
            input.initialCase   = settings.initial_case;            
            input.ndim          = mesh.ndim;
            input.coord         = mesh.coord;            
            input.bc            = obj.createBcInput(mesh);
            obj.input = input;
        end
        
        function i = getValue(obj)
            i = obj.input;            
        end
        
    end
    
    methods (Access = private, Static)
        
        function bc = createBcInput(mesh)
            bc = [];
            if ~isempty(mesh.dirichlet) && ~isempty(mesh.pointload)
                bc = unique([mesh.dirichlet(:,1); mesh.pointload(:,1)]);
            end
        end
        
    end
    
    
    
end

