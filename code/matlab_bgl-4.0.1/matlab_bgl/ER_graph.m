% This file generates an ER graph and does some subsequent computation.

% clear result of last computation
clear;
clc;

% Value assignment
N = 500;
d_av = 12;
p = d_av/(N-1);

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

%plot
plot(sorted_Deg);
hold on;
plot(eigen_Q);
hold off;
