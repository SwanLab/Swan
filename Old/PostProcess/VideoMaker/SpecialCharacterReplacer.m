classdef SpecialCharacterReplacer < handle
    
    methods (Access = public, Static)
        
        function oStr = replace(iStr)
            if ispc
                oStr = replace(iStr,'\','\\\\');
            elseif isunix
                oStr = iStr;
            elseif ismac
                oStr = iStr;
            end            
        end
        
    end

end

