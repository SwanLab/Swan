function ToGiD_post_OPT_internal_loop(dim,mstep,file_name,ngaus,element,problembsc,phifunct,young)

etype = element.type;
ptype = problembsc.problemtype;
ndime=dim.ndime; npnod=dim.npnod;nelem=dim.nelem; nnode=dim.nnode;
switch  etype
    case 'TRIANGLE'
        gtype = 'Triangle'; %gid type
    case 'QUAD'
        gtype = 'Quadrilateral';
    case 'HEXAHEDRA'
        gtype = 'Hexahedra';
end

% Leer tipo de post proceso (step,time)--------------
if isfield(problembsc,'kindpost')
    switch problembsc.kindpost
        case 'TIME'
            timestep = tp_time;
        case 'STEP'
            timestep = mstep;
        otherwise
            timestep = mstep;
            fprintf(1,'asignado postproceso "No. step".........\n')
    end
else
    timestep = mstep;   % por default no. de paso
end
% Escribe el fichero de resultados
res_file = strcat(file_name,'_',num2str(mstep),'.flavia.res');
fid = fopen(res_file,'w');

switch  etype
    case 'TRIANGLE'
        if nnode==3
            idxgp = [1 2 3]; job=2;
            gid_write_headerpost(fid,gtype,ngaus,job)
        elseif nnode==6
            idxgp = [];
        end
    case 'QUAD'
        if nnode==4
            idxgp = [1 2 3 4 ]; job =3;
            gid_write_headerpost(fid,gtype,ngaus,job)
        elseif nnode==8
            idxgp = [1 7 9 3 4 8 6 2 5]; job=1;
            gid_write_headerpost(fid,gtype,ngaus,job)
        end
end


if (problembsc.ppinfo(31)==1)
    nameres = 'PHI FUNCTION';
    gid_write_nodal_sclfield(fid,nameres,timestep,phifunct);
end

if (problembsc.ppinfo(30)==1)
    nameres = 'GAUSS YOUNG';
    gid_write_gauss_sclfield(fid,nameres,timestep,young,idxgp);
end


end
