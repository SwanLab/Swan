classdef HiperelasticityTesting < handle 

   properties (Access = private)
        fileName
        nsteps
        bcCase
    end

   methods (Access = public)

        function obj = HiperelasticityTesting()
            close all
            obj.init()
            s.nsteps = obj.nsteps;
            s.bcCase = obj.fileName;
            s.printing = true;
            s.bcCase = obj.bcCase;
            s.fileName = obj.fileName;
            s.meshGen  = 'EIFEMMesh';
            s.nSubdomains = [5,5];
%             h = HyperelasticProblem(s);
            h2 = HyperelasticProblem_refactoring(s);
            

        end

   end

   methods (Access = private)
       function init(obj)
%            obj.fileName = 'HoleDirich';
%            obj.nsteps = 20; 
%             obj.fileName = 'Metamaterial'; 
%            obj.nsteps = 75; %          
           obj.fileName = 'DEF_Q4auxL_1.mat';
            obj.nsteps = 50; 
           obj.bcCase = 'Traction';
       end

   end

end