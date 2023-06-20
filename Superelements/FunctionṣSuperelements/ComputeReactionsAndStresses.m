function  [totalReact, stress] = ComputeReactionsAndStresses(stressInputs)
    lambda = stressInputs.Fields.lambda;
    subC = stressInputs.subMatrices.C;
    subBoundMesh = stressInputs.meshDecomposed.subBoundMesh;
    plot = stressInputs.plot;
    type = stressInputs.type;

    
    stress = cell(size(length(lambda)));
    % Stresses
    j = 1;
    for i=1:length(lambda)
        j = j + 0.5;
        s.fValues = full(lambda{i});
        s.mesh = subBoundMesh(floor(j),1).mesh;
        if size(subC{1,1},1)==size(subC{1,1},2)
            stress(i) = {P1Function(s)};
        else
            stressFun = P0Function(s);
            stress(i) = {stressFun.project('P1')};
        end

    end

    % Reactions  
    totalReact = zeros(length(lambda),3);
    for i=1:length(lambda)
        VectorLambda = zeros(size(lambda{i},1)*size(lambda{i},2),1);
        VectorLambda(1:2:end-1) = lambda{i}(:,1);
        VectorLambda(2:2:end)   = lambda{i}(:,2);
        if type == "Two"
            React = subC{i,1}*VectorLambda;
            React = [React(1:2:end-1), React(2:2:end)];
            if plot
               plotReactions(React,subBoundMesh(i,1),2);
            end
            Center = subBoundMesh(i,1).mesh.coord(end,2)/2;
            RootMoment = sum(React(:,1)'*(subBoundMesh(i,1).mesh.coord(:,2) - Center));
        elseif type == "Three"
            if rem(i,2)==0
               React = subC{round(i/2),2}*VectorLambda;
               React = [React(1:2:end-1), React(2:2:end)];
               if plot
               plotReactions(React,subBoundMesh(round(i/2),2),2);
               end
               Center = subBoundMesh(round(i/2),2).mesh.coord(end,2)/2;
               RootMoment = sum(React(:,1)'*(subBoundMesh(round(i/2),2).mesh.coord(:,2) - Center)); 
            else
               React = subC{round(i/2),1}*VectorLambda;
               React = [React(1:2:end-1), React(2:2:end)];
               if plot
               plotReactions(React,subBoundMesh(round(i/2),1),2);
               end
               Center = subBoundMesh(round(i/2),1).mesh.coord(end,2)/2;
               RootMoment = sum(React(:,1)'*(subBoundMesh(round(i/2),1).mesh.coord(:,2) - Center));
            end    
        end

        totalReact(i,:) = [sum(React(:,1)) sum(React(:,2)) RootMoment];

    end
end
