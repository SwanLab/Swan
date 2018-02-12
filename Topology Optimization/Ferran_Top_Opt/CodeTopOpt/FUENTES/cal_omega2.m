function [volume] = cal_omega2(phifunct,dim,element,problembsc,coordinates)

ptype = problembsc.problemtype;
% basic dimensions
npnod=dim.npnod; ndime=dim.ndime; nelem=dim.nelem;
nnode=dim.nnode; 
% basic variables
[coordinatesn,coordinatesa] = init_coord(coordinates);


% calculo del volumen, omega
neres = 0;
lnods=zeros(nnode,nelem);
ephi=zeros(nnode,nelem);
for i=1:nnode
    lnods(i,:)= element.conectivities(:,i);
    ephi(i,:)= phifunct(lnods(i,:));
end
etype = element.type;
[posgp,weigp,ngaus] = cal_posgp_weigp(etype,ndime,nnode,element.ngaus);


evol=zeros(ngaus,nelem);
geometric_volum=zeros(ngaus,nelem);
%vol=zeros(npnod,1);

for igaus=1:ngaus
    [shape,~,~] = shape_deriv_functions(igaus,posgp,ptype,etype,nnode,neres);
    phigp=zeros(1,nelem);
    for inode=1:nnode
        phigp = phigp + shape(inode)*ephi(inode,:);
    end
    e=(phigp<=0);
    [~,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,ptype);
    evol(igaus,e) = weigp(igaus)*djacb(e)';
    geometric_volum(igaus,:) = weigp(igaus)*djacb';
end
% for inode=1:nnode
%     vol = vol + sparse(lnods(inode,:),1,evol(inode,:),npnod,1);
% end
% volume=sum(vol);

volume = sum(sum(evol))/sum(sum(geometric_volum));


end

