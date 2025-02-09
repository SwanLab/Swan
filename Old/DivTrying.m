function DivTrying





[du,uN] = Derivative('u')
[dv,vN] = Derivative('v')


divu = sum(diag(du))
divv = sum(diag(dv))

e = divu*divv

%du = 0.5*(du+du');
%dv = 0.5*(dv+dv');
%e = du(:)'*dv(:)


dV = 0.5;
nDim = 2;
nnode = 3;

for iDim = 1:nDim
    for jDim = 1:nDim
        for iNode = 1:nnode
            for jNode = 1:nnode
                idof = nnode*(iDim-1)+iNode;%nDim*(iNode-1)+iDim;
                jdof = nnode*(jDim-1)+jNode;                               
                Kij = diff(diff(e,uN(jDim,jNode)),vN(iDim,iNode));
                DD(idof,jdof) = Kij*dV;
            end
        end
    end
end
DD
eigDD = eig(DD)


end


function [du,uS] = Derivative(s)

x = sym('x','real');
y = sym('y','real');

uS(1,1) = sym([s,'x1'],'real');
uS(1,2) = sym([s,'x2'],'real');
uS(1,3) = sym([s,'x3'],'real');
uS(2,1) = sym([s,'y1'],'real');
uS(2,2) = sym([s,'y2'],'real');
uS(2,3) = sym([s,'y3'],'real');



N(1,1) = 1 - x - y;
N(2,1) = x;
N(3,1) = y;

u = uS*N;

%ux = ux1*N1 + ux2*N2 + ux3*N3;
%uy = uy1*N1 + uy2*N2 + uy3*N3;

duxdx = diff(u(1),x);
duxdy = diff(u(1),y);
duydx = diff(u(2),x);
duydy = diff(u(2),y);

du(1,1) = duxdx;
du(1,2) = duydx;
du(2,1) = duxdy;
du(2,2) = duydy;

end