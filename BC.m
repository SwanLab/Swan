classdef BC
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        displacements
        Fpointload
        iN
        iD
    end
    
    methods
        % Constructor
        function obj = BC(nunkn,filename)
            if nargin ~= 0 
                [obj.displacements,obj.Fpointload] = Preprocess.getBC(filename);
                
                for i = 1:length(obj.displacements(:,1))
                    obj.iD(i) = obj.displacements(i,1)*nunkn - nunkn + obj.displacements(i,2);
                end
                
                for i = 1:length(obj.Fpointload(:,1))
                    obj.iN(i)= obj.Fpointload(i,1)*nunkn - nunkn + obj.Fpointload(i,2);
                end
                
            end
        end
    end
    
end

