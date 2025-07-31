clear; clc;

params.state_range = [-1, 1];
params.actions = [-1, 0, 1];
params.n_actions = length(params.actions);
params.walk_step = 0.05;
params.arbitrariness = params.walk_step / 5;

params.tiles_per_dim = 100;     % or [100] for 1D
params.num_tilings = 8;
params.n_features = params.tiles_per_dim * params.num_tilings * params.n_actions;

params.alpha = 0.1 / params.num_tilings;
params.gamma = 0.99;
params.lambda = 0.9;
params.epsilon = 0.8;
params.n_episodes = 1000;
params.offsets = generateDiagonalOffsets(1, params.num_tilings);

% Epsilon decay strategy
params.epsilon_strategy = 'decay';   % or 'linear_decay'
params.epsilon_decay = 0.995;
%params.min_epsilon = 0.001;

env = create1DEnv(params);

% Define your policy (ε-greedy)
policy = @(s, w, epsilon, params) epsilonGreedy(s, w, epsilon, params, @getActiveTiles);

% Run SARSA(λ)
weights = qlearning_lambda(env, policy, @getActiveTiles, params);

%% Visualization
plotValueFunction(weights, params)
