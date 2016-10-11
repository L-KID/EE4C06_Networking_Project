% This file generates an ER graph and does some subsequent computation.

% TO DO: conduct multiple experiments and plot the average value to
% eliminate uncertainty

% Clear result of last computation
clear;
clc;

%% Parameters declaration
% Value assignment
N = 500;
m = 6;
m0 = 7;
n = 100; % Number of ER graphs(should be 100000)

% Generate the ER graph
%A = erdos_reyni(N, p);

% Compute the degree vector
u = ones(500, 1); % an all-one vector with 500 rows and 1 column

% Define 3 cells to store arrays
Deg_bin = cell(n,1);
Deg_org = cell(n,1);
eigen_Q = cell(n,1);

%% Computation with 100000 random graph for both Degree and Eigenvalue
%Store all data in the cell, including unique array, and the original value
for i = 1:1:n
    A = scalefree(N, m0, m);
    %Deg = A * u; % Deg is also 500-by-1
    Deg = A*u;
    % Plot distribution of degree and eigenvalues for 1 ER graph
    if i == 1
        % Compute the Laplacian matrix and its eigenvalues
        sorted_Deg = sort(Deg);
        diag_matrix = diag(Deg); % denoted as Delta in the report
        Q1 = diag_matrix - A;
        eigen_Q1 = eig(Q1);
        
        %plot
        plot(sorted_Deg);
        hold on;
        plot(eigen_Q1);
        xlabel('k');
        ylabel('Degree and the Laplacian eigenvalues');
        title('The degree vector and the Laplacian eigenvalues of a graph');
        legend('Laplacian eigenvalues u_{(k)}','Ordered degree d_{(k)}');
        hold off;
    end
    Diag_matrix = diag(Deg);
    Q = Diag_matrix -A;
    eigen_Q(i,1) = {eig(Q)};
    Deg_org(i,1) = {Deg};
    Deg_bin(i,1) = {unique(Deg)};
end

% compute all unique value from all data, and hist data
Deg_all = unique(cell2mat(Deg_bin));
Deg_hist = hist(cell2mat(Deg_org),Deg_all);

%% Plots
figure
plot(Deg_all, Deg_hist/N); % divided by N to show probability
hold on;

% Laplacian eigenvalues distribution
rounded_eigen_Q = round(cell2mat(eigen_Q));
rounded_eigen_Q_bin = unique(rounded_eigen_Q);
rounded_eigen_Q_hist = hist(rounded_eigen_Q, rounded_eigen_Q_bin);

plot(rounded_eigen_Q_bin, rounded_eigen_Q_hist/N ); % divided by N to show probability
xlabel('x');
ylabel('f_u(x)');
title('The distribution of degree and Laplacian eigenvalues');
legend('degree','Laplacian eigenvalues');
hold off;
