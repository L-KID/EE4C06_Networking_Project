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
xlabel('k');
ylabel('Degree and the Laplacian eigenvalues');
title('The degree vector and the Laplacian eigenvalues of a graph');
legend('Laplacian eigenvalues u_{(k)}','Ordered degree d_{(k)}');
hold off;

% The following code of plotting the distribution of the degrees might be
% INCORRECT!
% Ref: http://stackoverflow.com/questions/7907017/count-occurrences-on-a-array-using-matlab
% Besides, normalization is not considered.
Deg_bin = unique(Deg);
Deg_hist = hist(Deg, Deg_bin);

figure
plot(Deg_bin, Deg_hist);
hold on;

% Laplacian eigenvalues distribution
rounded_eigen_Q = round(eigen_Q);
rounded_eigen_Q_bin = unique(rounded_eigen_Q);
rounded_eigen_Q_hist = hist(rounded_eigen_Q, rounded_eigen_Q_bin);

plot(rounded_eigen_Q_bin, rounded_eigen_Q_hist);
hold off;

