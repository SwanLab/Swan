function [  ] = printinfo_internal_loop(OK,problembsc,file_name,iter,coordinatesn,element,...
    dim,phifunct,young)
% Postproceso

ndime=dim.ndime; npnod=dim.npnod; nndof=npnod*ndime;
nelem=dim.nelem; nnode=dim.nnode;
[~,~,ngaus] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);
if (OK==1)
     contactb=[];
    if mod(iter,problembsc.frqpostpros) ==0
        ToGID (file_name,iter,coordinatesn,element,contactb,problembsc,nnode);
        ToGiD_post_OPT_internal_loop(dim,iter,file_name,ngaus,element,problembsc,phifunct,young)

    end
end

end

