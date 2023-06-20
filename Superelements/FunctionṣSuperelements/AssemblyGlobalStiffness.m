function [globalK] = AssemblyGlobalStiffness(AKinputs)
   subK = AKinputs.subK;
   subdomains = AKinputs.subdomains;

   globalK = subK{1};
   for i=2:subdomains
        globalK = [                globalK                 sparse(size(globalK,1),size(subK{i},2));
                   sparse(size(subK{i},1),size(globalK,2))                   subK{i}              ];
   end
end

