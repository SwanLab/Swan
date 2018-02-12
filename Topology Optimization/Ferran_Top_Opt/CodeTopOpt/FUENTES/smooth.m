function [nf] = smooth(nelem,npnod,ndime,nnode,gf,Msmooth,element,coordinatesn,coordinatesa,problembsc)
% gf : gauss funtion f
% nf : nodal function f
% job = 1: Msmooth is the lumped matrix, nf = int(gf,)/M
% job = 2: Msmooth is th full mass matrix

job = problembsc.smoothing_proc;
ptype = problembsc.problemtype;
for i=1:nnode
    lnods(i,:)= element.conectivities(:,i);
end
etype = element.type;
neres = 0;
% compute rhs int(f*Na)
rhs = zeros(npnod,1);
erhs=zeros(nnode,nelem);
[posgp,weigp,ngaus] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);
for igaus=1:ngaus    
    [cartd,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,ptype);
    [shape,deriv,heslo] = shape_deriv_functions(igaus,posgp,ptype,etype,nnode,neres);
    dvolu = weigp(igaus)*djacb;
    for inode=1:nnode
        erhs(inode,:) = erhs(inode,:) + shape(inode)*gf(igaus,:).*dvolu';
    end
end
for inode=1:nnode
    rhs = rhs + sparse(lnods(inode,:),1,erhs(inode,:),npnod,1);
end

%compute nodal function
if (job==1)
    nf = rhs./Msmooth;
    if (nnode>4) 
        fprintf(1,'WARNING: LUMPED MATRIX SHOULD NOT BE USE WITH THIS COMBINATION OF ELEM AND NNODE');
    end
elseif (job==2)
    nf = Msmooth\rhs;
end

end

