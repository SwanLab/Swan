function PrintResults(PrintInputs)
    subMesh = PrintInputs.subMesh;
    U = PrintInputs.U;
    type = PrintInputs.type;

    subdomains = size(subMesh,1);
    for i=1:subdomains
        subMesh(i).coord(:,1) = subMesh(i).coord(:,1) + (i-1)*0.002;
        s = [];
        s.fValues = full(U{i});
        s.fValues = [s.fValues zeros(size(s.fValues,1),1)]; % In order to plot deformed (z=0)
        s.mesh    = subMesh(i);
        p.filename = ['domain',char(string(i))];
        ResultSubDom = P1Function(s);
        ResultSubDom.print(p);
    end
    
%     if type == "Three"
%         interMesh = PrintInputs.interMesh;
%         interU = PrintInputs.interU;
%         interfaces = size(interMesh,1);
%         for j=1:interfaces
%             interMesh(j).mesh.coord(:,1) = interMesh(j).mesh.coord(:,1) + (j-1)*0.002;
%             s = [];
%             s.fValues = full(interU{j});
%             s.fValues = [s.fValues zeros(size(s.fValues,1),1)]; % In order to plot deformed (z=0)
%             s.mesh    = interMesh(j).mesh;
%             p.filename = ['interface',char(string(j))];
%             ResultSubDom = P1Function(s);
%             ResultSubDom.print(p);
%         end
%     end


    

end