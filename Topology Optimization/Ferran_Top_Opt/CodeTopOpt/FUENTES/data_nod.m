function [npnod,nndof,ndime,nunkn,nstre]=data_nod(coordinates,type,element,problembsc)
% Esta funci�n entrega informaci�n a la rutina principal de MatFEM acerca
% del numero total de nodos del problema y del numero de grados de libertad
% dependiendo de si la formulaci�n es estandar o mixta.
ptype = problembsc.problemtype;
ftype = problembsc.phisical_type;


switch ptype
    case '1D'
        ndime = 1;  
        nstre = 1;    
    case '2D'
        ndime = 2;
        switch ftype
            case {'ELASTIC'}
                switch element.material.subtype
                    case {'PLANESTRAIN','PLANESTRES'}
                        nstre = 3;
                    case 'AXI'
                        nstre = 4;
                end
                
            case {'THERMAL'}
                nstre = 2;
        end
    case '3D'
        ndime = 3;
        switch ftype
            case {'ELASTIC'}
                nstre = 6;
            case {'THERMAL'}
                nstre = 3;
        end
        
end


switch ftype
    case {'ELASTIC'}
        nunkn = ndime;
        switch type
            case 'BAR'
                npnod  = size(coordinates,1);        % Number of nodes
                nndof  = npnod;                    % Number of total DOF
            case 'TRIANGLE'
                npnod  = size(coordinates,1);        % Number of nodes
                nndof  = 2*npnod;                    % Number of total DOF
            case 'LINEAR_TRIANGLE_MIX'
                npnod  = size(coordinates,1);        % Number of nodes
                nndof  = 3*npnod;                    % Number of total DOF
            case 'LINEAR_TRIANGLE_MIX_COUPLED'
                npnod  = size(coordinates,1);        % Number of nodes
                nndof  = 3*npnod;                    % Number of total DOF
            case {'QUAD'}
                npnod  = size(coordinates,1);        % Number of nodes
                nndof  = 2*npnod;                    % Number of total DOF
            case 'HEXAHEDRA'
                npnod  = size(coordinates,1);        % Number of nodes
                nndof  = 3*npnod;                    % Number of total DOF
            otherwise
                error('Formulaci�n no ha sido implementada')
        end
        
        
        
    case {'THERMAL'}
        nunkn = 1;
        switch type
            case 'BAR'
                npnod  = size(coordinates,1);        % Number of nodes
                nndof  = npnod;                    % Number of total DOF
            case 'TRIANGLE'
                npnod  = size(coordinates,1);        % Number of nodes
                nndof  = npnod;                    % Number of total DOF
            case {'QUAD'}
                npnod  = size(coordinates,1);        % Number of nodes
                nndof  = npnod;                    % Number of total DOF
            case 'HEXAHEDRA'
                npnod  = size(coordinates,1);        % Number of nodes
                nndof  = npnod;                    % Number of total DOF
            otherwise
                error('Formulaci�n no ha sido implementada')
        end
end




end