% This file generates an ER graph, computes s, and plot the distribution of
% s.

% Clear result of last computation
clear;
clc;

% Value assignment
N = 500;
d_av = 12;
p = d_av/(N-1);
ITERATION_TIME = 500;
iteration = 0;
s_vector = zeros(ITERATION_TIME, 1);

while iteration < ITERATION_TIME
    iteration = iteration + 1;

    % Generate the ER graph
    A = erdos_reyni(N, p);

    % Compute the degree vector
    u = ones(500, 1); % an all-one vector with 500 rows and 1 column
    Deg = A * u; % Deg is also 500-by-1

    % Compute the Laplacian matrix and its eigenvalues
    sorted_Deg = sort(Deg);
    diag_matrix = diag(Deg); % denoted as Delta in the report
    Q = diag_matrix - A;
    eigen_Q = eig(Q);

    % Compute s
    % Ref: http://nl.mathworks.com/matlabcentral/newsreader/view_thread/12610
    S_temp = eigen_Q - sorted_Deg;
    S_temp_positive = find(S_temp >= 0);
    s = length(S_temp_positive);
    s_vector(iteration) = s;
end

% Simulation
s_vector_bin = unique(s_vector);
s_vector_hist = hist(s_vector, s_vector_bin);
plot(s_vector_bin, s_vector_hist/ITERATION_TIME); % divided by ITERATION_TIME to show probability

figure
semilogy(s_vector_bin, s_vector_hist/ITERATION_TIME); % divided by ITERATION_TIME to show probability

% Fitting
figure
pd = fitdist(s_vector,'Normal');
y = pdf(pd, s_vector_bin);
semilogy(s_vector_bin, y, 'LineWidth', 2);