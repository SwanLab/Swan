function policyFunction = createPolicyFunction(policyType)
    switch lower(policyType)
        case 'epsilongreedy'
            policyFunction = @(state, w, epsilon, params, getActiveTiles) ...
                             epsilonGreedy(state, w, epsilon, params, getActiveTiles);
        case 'softmax'
            policyFunction = @(state, w, params, getActiveTiles) ...
                             softmaxPolicy(state, w, params, getActiveTiles); 
            % Note: softmaxPolicy does not need epsilon, so we pass ~
        otherwise
            error('Unknown policy type');
    end
end