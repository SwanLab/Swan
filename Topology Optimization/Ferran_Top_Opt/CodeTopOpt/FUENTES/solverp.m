function [d_u,free] = solverp(d_u,K,pnods,fextlod,fixed,dim)

n = size(d_u,1);
nunkn = dim.nunkn;
nlib = size(pnods(1,:),2);
pglib = zeros(nlib*nunkn,1);
qglib = zeros(nlib*nunkn,1);
for iunkn = 1:nunkn
    index_glib = nlib*(iunkn - 1) + [1:nlib];
    pglib(index_glib,1) = (pnods(1,:)-1)*nunkn + iunkn;
    qglib(index_glib,1) = (pnods(2,:)-1)*nunkn + iunkn;
end

% pglib = [nunkn*pnods(1,:)-1,nunkn*pnods(1,:)]';
% qglib = [nunkn*pnods(2,:)-1,nunkn*pnods(2,:)]';

list = [pglib; qglib; fixed];
free = setdiff([1:1:n],list);

newK = [K(free,free) , K(free,pglib)+K(free,qglib);...
  K(pglib,free)+K(qglib,free), K(pglib,pglib)+K(pglib,qglib)+K(qglib,pglib)+ K(qglib,qglib)];

rhs1 = fextlod(pglib)+fextlod(qglib);
rhs = [fextlod(free) ; rhs1];
usol = newK \rhs;
d_u(free)=usol(1:1:size(free,2));
d_u(pglib)=usol(size(free,2)+1:1:size(newK,1));
d_u(qglib)=d_u(pglib);

%fprintf(1,'COND NUMBER %d \n',condest(newK));

end

