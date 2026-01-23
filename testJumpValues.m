mesh = CohesiveMesh();
u = LagrangianFunction.create(mesh,2,1);

fValuesNew = [1 1 1 1 1 1; 1 1 11 1]
u.setFValues(fValuesNew)




function computeJump(mesh,u)

end