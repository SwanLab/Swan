classdef Material
    %Material Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = {?Material_Elastic,?Material_Stokes}, SetAccess = private)
        nelem
    end
    
    methods (Access = protected)
        function obj = Material(nelem)
            obj.nelem = nelem;
        end
    end
    
    methods (Static) %(Access = ?Physical_Problem, Static)
        function material = create(geometry,mesh)
            ptype=mesh.ptype;
            pdim=mesh.pdim;
            nelem=geometry(1).interpolation.nelem;
            connec=mesh.connec;
            cartd=geometry(1).cartd;
            nnode=geometry(1).interpolation.nnode;
            coord=mesh.coord;
            switch ptype
                case 'ELASTIC'
                    % !! IT HAS BEEN ASSUMED THAT THERE'S ONLY ISOTROPIC MATERIALS.
                    %    THIS HAS TO BE CHANGED FOR THE OPT TOP PROBLEM !!
                    switch pdim
                        case '2D'
                            material = Material_Elastic_ISO_2D(nelem);
                        case '3D'
                            material = Material_Elastic_ISO_3D(nelem);
                    end
                    
                case 'HYPERELASTIC'
                    switch pdim
                        case '2D'
                            material = Material_Hyperelastic_2D(nelem,connec,cartd,nnode,coord);
                    end
                    
                case {'THERMAL', 'DIFF-REACT'}
                    error('Still not implemented.')
                case 'Stokes'
                    material = Material_Stokes(nelem);
                otherwise
                    error('Invalid ptype.')
            end
        end
    end
    
end