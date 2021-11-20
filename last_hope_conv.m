% convergence plot

NC = 4; % #Cells which equals #BSs
NU = 10; % #USers in each cell.
P = 16;
nvar = 1.9905e-08; % Noise Variance ?????
epsilon = 1e-5; % For convergence test.
inner_radius = 500; 
minR_ratio = 0.01;
numIter = 2000;
num_reals = 1000; %# of channel realizations
alpha_rng = 1:10; % the range of the alpha computation indices
seed = 1;

% system model
% Load the presaved channel parameters.
fileName = sprintf('channels_for_NU_fixed/Channels%dx%dpower%d.mat', NC, NU, P);
load(fileName,'H', 'in', 'D'); 
last_hope

plot_conv(alpha_rng, conv, WR_vs_iter);
