classdef DIM
    %DIM Summary of this class goes here
    %   Detailed explanation goes here
    properties (GetAccess = public, SetAccess = private)
        nstre
        %       nnode
        nunkn
    end
    properties (GetAccess = {?Physical_Problem,?PhysicalVariables,?Postprocess}, SetAccess = private)
        ndim
        %       ngaus
    end
    
    methods (Access = ?Physical_Problem)
        function obj = DIM(ptype,pdim)
            switch ptype
                case {'ELASTIC','HYPERELASTIC'}
                    switch pdim
                        case '2D'
                            obj.ndim = 2;
                            obj.nunkn = 2;
                            obj.nstre = 3;
                        case '3D'
                            obj.ndim = 3;
                            obj.nunkn = 3;
                            obj.nstre = 6;
                    end
                case 'THERMAL'
                    error('Still not implemented.')
                case 'Stokes'
                    switch pdim
                        case '2D'
                            obj.ndim = 2;
                            nunkn_u = 2;
                            nunkn_p = 1;
                            obj.nunkn = [nunkn_u nunkn_p];
                            obj.nstre = 0;
                        case '3D'
                            obj.ndim = 3;
                            obj.nunkn = 3;
                            obj.nstre = 6;
                    end
            end
        end
    end
end



