s.testName = 'test_dehomogenizingSingularities';            
test = DehomogenizingSingularitiesTest(s);

% 
% % 
% m = obj.backgroundMesh;
% s.type    = m.type;
% s.connec  = m.connec;
% s.fValues(1,:,:) = reshape(ls,3,[]);
% f         = P1DiscontinuousFunction(s);
% f.plot(m); shading flat
% 
% 
% m = obj.mesh;
% s.type    = m.type;
% s.connec  = m.connec;
% s.fValues(1,:,:) = (phi');
% f         = P1DiscontinuousFunction(s);
% f.plot(m); shading flat
% 
% 
% m = obj.mesh;
% s.type    = m.type;
% s.connec  = m.connec;
% s.fValues(1,:,:) = reshape(obj.phi(:,1),3,[]);
% f         = P1DiscontinuousFunction(s);
% f.plot(m); shading flat