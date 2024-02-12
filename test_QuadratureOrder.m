o1.ngaus = 1;
o1.posgp(:,1) = [0,0];
o1.weigp = 4;

o2.ngaus = 4;
a =  0.577350269189626;
o2.posgp(:,1) = [-a,-a];
o2.posgp(:,2) = [+a,-a];
o2.posgp(:,3) = [-a,+a];
o2.posgp(:,4) = [+a,+a];
o2.weigp =  [1,1,1,1];

posgl(1) =-0.774596669241483;
posgl(2) = 0.0;
posgl(3) = 0.774596669241483;
weigl(1) = 0.555555555555556;
weigl(2) = 0.888888888888889;
weigl(3) = 0.555555555555556;
o3.ngaus = 9;
igaus = 0;
nlocs = 3;
for ilocs = 1:nlocs
    for jlocs = 1:nlocs
        igaus = igaus+1;
        o3.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
        o3.posgp(1,igaus) = posgl(ilocs);
        o3.posgp(2,igaus) = posgl(jlocs);
    end
end

posgl(1) = sqrt(3/7-2/7*sqrt(6/5));
posgl(2) =-sqrt(3/7-2/7*sqrt(6/5));
posgl(3) = sqrt(3/7+2/7*sqrt(6/5));
posgl(4) =-sqrt(3/7+2/7*sqrt(6/5));
weigl(1) = (18+sqrt(30))/36;
weigl(2) = (18+sqrt(30))/36;
weigl(3) = (18-sqrt(30))/36;
weigl(4) = (18-sqrt(30))/36;
o4.ngaus = 16;
igaus = 0;
for ilocs = 1:4
    for jlocs = 1:4
        igaus = igaus+1;
        o4.weigp(  igaus) = weigl(ilocs)*weigl(jlocs);
        o4.posgp(1,igaus) = posgl(ilocs);
        o4.posgp(2,igaus) = posgl(jlocs);
    end
end


syms x y
f(x,y) = (x+1)^8*(y+1)^7;

intAn = int(int(f,x,-1,1),y,-1,1);

oArray = [o1 o2 o3 o4];
intNum = zeros(size(oArray));
for o = 1:length(oArray)
    for i = 1:oArray(o).ngaus
        intNum(o) = intNum(o) + oArray(o).weigp(i)*subs(subs(f,'x',oArray(o).posgp(1,i)),'y',oArray(o).posgp(2,i));
    end
end

display(intAn)
display(intNum)


mesh = UnitQuadMesh(2,2);

sAF.fHandle = @(x) x(1,:,:);
sAF.ndimf   = 1;
sAF.mesh    = mesh;
xFun = AnalyticalFunction(sAF);

p1fun = xFun.project('P1');
p2fun = xFun.project('P2');
p3fun = xFun.project('P3');