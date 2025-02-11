function ExperimentingGraph
nmax = 280;
G = graph();
G = addnode(G,nmax);
plot(G)

allEdges = nchoosek(1:nmax,2);
allEdgesN = allEdges;
nEdgeMax = size(allEdges,1);
sparsity = 0.05;
nEdge = round(sparsity*nEdgeMax);
for i = 1:nEdge
    nPosibleEdges = size(allEdgesN,1);
    newEdge = randi(nPosibleEdges);
    edge    = allEdgesN(newEdge,:);
    allEdgesN(edge,:) = [];
    G = addedge(G,edge(1),edge(2));
    
    if mod(i,round(nEdge/100)) == 0
        figure(1)
        clf
        subplot(4,1,1)
        plot(G)
        d = centrality(G,'degree');
        subplot(4,1,2)
        histogram(d)
        xlim([0 0.02*nEdge])
        subplot(4,1,3)
        t = computeTriangularValues(G);
        histogram(t)
        C = computeTransition(d,t);
        subplot(4,1,4)
        histogram(C)
        drawnow
    end
end

end

function t = computeTriangularValues(G)
A = G.adjacency;
A3 = A^3;
t = (diag(A3)/2);
end


function C = computeTransition(d,t)
C = 2*t./(d.*(d-1));
mean(C)
isCNan = isnan(C);
C(isCNan) = 0;
end