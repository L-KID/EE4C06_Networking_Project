% This file generates an ER graph and does some subsequent computation.

% Clear result of last computation
clear;
clc;

%% Parameters declaration
% Value assignment
N = 500;
d_av = 12;
p = d_av/(N-1);
u = ones(500, 1); % create an all-one vector with 500 rows and 1 column
num_simulation = 500; % Number of simulation times (should be 100000)

% Define 3 cells to store arrays
Deg_bin = cell(num_simulation,1);
Deg_org = cell(num_simulation,1);
eigen_Q = cell(num_simulation,1);

%% Computation with 100000 random graph for both Degree and Eigenvalue
%Store all data in the cell, including unique array, and the original value
for i = 1:1:num_simulation
    % Generate the ER graph
    A = erdos_reyni(N, p);
    Deg = A * u; % Deg is also 500-by-1
    % Plot distribution of degree and eigenvalues for 1 ER graph
    if i == 1 % for Q1, only one simulation is sufficient
        % Compute the Laplacian matrix and its eigenvalues
        sorted_Deg = sort(Deg);
        diag_matrix = diag(Deg);
        Q1 = diag_matrix - A;
        eigen_Q1 = eig(Q1);
        
        %plot
        plot(sorted_Deg);
        hold on;
        plot(eigen_Q1);
        xlabel('k');
        ylabel('Degree and the Laplacian eigenvalues');
        title('The degree vector and the Laplacian eigenvalues of a graph');
        legend('Ordered degree d_{(k)}','Laplacian eigenvalues u_{(k)}');
        hold off;
    end
    Diag_matrix = diag(Deg);
    Q = Diag_matrix - A;
    eigen_Q(i,1) = {eig(Q)};
    Deg_org(i,1) = {Deg};
    Deg_bin(i,1) = {unique(Deg)};
end

% compute all unique value from all data, and hist data
Deg_all = unique(cell2mat(Deg_bin));
Deg_hist = hist(cell2mat(Deg_org), Deg_all);

%% Plots
figure
plot(Deg_all, Deg_hist/(N*num_simulation),'-o'); % divided by N to show probability
hold on;

% Laplacian eigenvalues distribution
rounded_eigen_Q = round(cell2mat(eigen_Q));
rounded_eigen_Q_bin = unique(rounded_eigen_Q);
rounded_eigen_Q_hist = hist(rounded_eigen_Q, rounded_eigen_Q_bin);

plot(rounded_eigen_Q_bin, rounded_eigen_Q_hist/(N*num_simulation),'-*' ); % divided by N to show probability
xlabel('x');
ylabel('f_u(x)');
title('The distribution of degree and Laplacian eigenvalues');
legend('degree','Laplacian eigenvalues');
hold off;
