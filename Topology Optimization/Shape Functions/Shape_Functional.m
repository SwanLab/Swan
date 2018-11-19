classdef Shape_Functional < handle
    properties
        value
        gradient
        target_parameters = struct;
        filter
        Msmooth
        dvolu
    end
    
    methods (Access = public)
        function obj = Shape_Functional(settings)
            obj.filter = Filter.create(settings);
            obj.filter.setupFromGiDFile(settings.filename,settings.ptype);
            diffReacProb = DiffReact_Problem;
            diffReacProb.setupFromGiDFile(settings.filename);
            diffReacProb.preProcess;
            obj.Msmooth = diffReacProb.element.M;
            obj.dvolu = diffReacProb.geometry.dvolu;
        end
    end
    
    methods (Abstract, Access = public)
        computeCostAndGradient(obj, x)
    end
end
