function [phifunct_n,imesh,fext,element,fixnodes,problembsc,coordinates,dim,Msmooth,Stiff_smooth,emass] = remeshing(remesh,phifunct_n,imesh,fext,element,fixnodes,problembsc,coordinates,dim,Msmooth,Stiff_smooth,emass,file_name,Vfrac,lagrange_type,function_type,file_write,conectivities,TYPE)

switch remesh
    case 1
        
        [fext,element2,fixnodes,problembsc,new_coordinates,~] = call_read_data(file_name,Vfrac,imesh+1,lagrange_type,function_type,file_write,TYPE);
        element = replace_struct(element,element2);
                                                               
        [phifunct_n] = transfer_info(conectivities,coordinates,new_coordinates,phifunct_n);
        coordinates = new_coordinates;
                
        [dim.npnod,dim.nndof,dim.ndime] = data_nod(coordinates,element.type,element,problembsc);
        [dim.nelem,dim.nnode,dim.neleq] = data_elem(element.conectivities, element.type);
        [coordinatesn,coordinatesa] = init_coord(coordinates);
        [Msmooth,emass] = mass_matrix(dim,problembsc,element,coordinatesn,coordinatesa);
        [Stiff_smooth] = stiff_unitary_triang( dim,element,problembsc,coordinatesn,coordinatesa,2,1,dim.npnod*1);
                                               
        imesh = imesh + 1;
    end
end

function struct = replace_struct(struct,struct2)

fields = fieldnames(struct2);
for ifield = 1:length(fields)
    struct.(fields{ifield}) = struct2.(fields{ifield});
end
end


