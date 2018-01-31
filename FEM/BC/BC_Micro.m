classdef BC_Micro < BC
    %BC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Physical_Problem,?Element,?Filter}, SetAccess = private)
        pnodes
    end
    
    methods (Access = ?Physical_Problem)
        % Constructor
        function obj = BC_Micro(nunkn,filename,coords,ptype,ndim)
            obj@BC(nunkn,filename);
            obj.computeFixedNodesValues(ptype,ndim);
            obj.computeiDiN(nunkn);
            obj.pnodes = Preprocess.getPeriodicBC(coords);
        end
    end
    
    methods (Access = private)
        function obj=computeFixedNodesValues(obj,ptype,ndim)
            switch ptype
                case 'ELASTIC'
                    ifix=1;
                    for i=1:size(obj.fixnodes_perimeter,1)
                        for j=1:ndim
                            obj.fixnodes(ifix,1)=obj.fixnodes_perimeter(i,1); % node
                            obj.fixnodes(ifix,2)=j; % idim
                            obj.fixnodes(ifix,3)=0; % U_imp     
                            ifix=ifix+1;
                        end
                    end
%                 case 'THERMAL'          %HAS TO BE REVISED FOR THERMAL
%                     ifix=1;
%                     for i=1:size(obj.fixnodes_perimeter,1)
%                         obj.fixnodes(ifix,1)=obj.fixnodes_perimeter(i); % node
%                         obj.fixnodes(ifix,2)=1; % idim
%                         obj.fixnodes(ifix,3)=0; % U_imp
%                         ifix=ifix+1;
%                     end
            end
        end
    end
end


