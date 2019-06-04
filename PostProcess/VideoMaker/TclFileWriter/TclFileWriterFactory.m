classdef TclFileWriterFactory < handle
    
    methods (Access = public, Static)
        
        function obj = create(cParams)
            switch cParams.type 
                case 'LevelSet'
                    obj = TclFileWriter_LevelSet(cParams);
                case {'Density','RegularizedDensity'}
                    obj = TclFileWriter_Density(cParams);
            end            
        end
        
        
    end
end