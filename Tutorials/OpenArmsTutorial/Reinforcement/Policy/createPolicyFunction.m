function policyFunction = createPolicyFunction(policyType)
    switch lower(policyType)
        case 'epsilongreedy'
            policyFunction = @(state, w, epsilon, params, activeTiles) ...
                             epsilonGreedy(state, w, epsilon, params, activeTiles);
        case 'softmax'
            policyFunction = @(state, w, params, activeTiles) ...
                             softmaxPolicy(state, w, params, activeTiles); 
        otherwise
            error('Unknown policy type');
    end
end