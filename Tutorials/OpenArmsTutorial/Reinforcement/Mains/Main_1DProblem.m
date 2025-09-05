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
params.gamma = 1;
params.lambda = 0.8;
params.epsilon = 0.3;
params.n_episodes = 2000;
params.offsets = generateDiagonalOffsets(1, params.num_tilings);

% Actor-critic hyperparameters
params.alpha_c = 0.05/params.num_tilings;    % critic learning rate
params.alpha_a = 0.05/params.num_tilings;    % actor learning rate

% Epsilon decay strategy
params.epsilon_strategy = 'decay';   % or 'linear_decay'
params.epsilon_decay = 0.999;
%params.min_epsilon = 0.001;

params.initial_state = [-1];

% --- Environment ---
env = OneDEnv(params);

% --- Active Tiles ---
activeTiles = ActiveTiles(params);

% --- Policy ---
policyType = 'epsilonGreedy';
policyFunction = createPolicyFunction(policyType);

% --- Agent ---
%agent = Sarsa(params.n_features);
agent = Q_learning(params.n_features, activeTiles);

% --- Training ---
weights = TD_lambda(env, policyFunction, activeTiles, params, agent);


%% Visualization
plotValueFunction(weights, params, activeTiles)
