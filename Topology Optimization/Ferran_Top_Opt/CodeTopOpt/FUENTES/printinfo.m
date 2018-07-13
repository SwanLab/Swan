function [  ] = printinfo(problembsc,file_name,iter,coordinatesn,element,...
    dim,gamma_gp,gamma_print,g_gp,phi,structural_values)
% Postproceso

d_u = structural_values.d_u;

ndime=dim.ndime; npnod=dim.npnod; nndof=npnod*ndime;
nelem=dim.nelem; nnode=dim.nnode;
[posgp,weigp,ngaus] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);

switch problembsc.phisical_type
    case {'ELASTIC'}
        
        ptype = problembsc.problemtype;
        switch ptype
            case '1D'
                tdisp(1:npnod,1)=d_u;
            case '2D'
                tdisp(1:npnod,1)=d_u(1:2:nndof);
                tdisp(1:npnod,2)=d_u(2:2:nndof);
            case '3D'
                tdisp(1:npnod,1)=d_u(1:3:nndof);
                tdisp(1:npnod,2)=d_u(2:3:nndof);
                tdisp(1:npnod,3)=d_u(3:3:nndof);
        end
        
    case {'THERMAL'}
        tdisp(1:npnod,1)=d_u;
end

contactb=[];
if mod(iter,problembsc.frqpostpros) ==0
    ToGID (file_name,iter,coordinatesn,element,contactb,problembsc,nnode);
    ToGiD_post_OPT(dim,iter,file_name,ngaus,element,tdisp,coordinatesn,gamma_gp,gamma_print,g_gp,phi,structural_values,problembsc);
end


end

