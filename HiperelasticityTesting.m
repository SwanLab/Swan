classdef HiperelasticityTesting < handle 

   properties (Access = private)
        fileName
        nsteps
    end

   methods (Access = public)

        function obj = HiperelasticityTesting()
            close all
            obj.init()
            s.nsteps = obj.nsteps;
            s.bcCase = obj.fileName;
            s.printing = false;
            h = HyperelasticProblem(s);
            

        end

   end

   methods (Access = private)
       function init(obj)
           obj.fileName = 'HoleDirich'; %'Metamaterial'
           obj.nsteps = 20;
       end

   end

end