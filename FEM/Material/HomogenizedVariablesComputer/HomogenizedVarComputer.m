classdef HomogenizedVarComputer < handle


    methods (Access = public, Static)
        
        function obj = create(cParams)
           f = HomogenizedVarComputerFactory();
           obj = f.create(cParams); 
        end        
        
    end

end