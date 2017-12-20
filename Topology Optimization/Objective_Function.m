classdef Objective_Function < handle
    properties
        value
        gradient
    end
    methods
        computeFunction(obj)
        computeGradient(obj)
    end
end
