% This file generates an ER graph and does some subsequent computation.

% Clear result of last computation
clear;
close all;
clc;

%% Parameters declaration
% Value assignment
N = 500;
d_av = 12;
p = d_av/(N-1);
u = ones(500, 1);
num_simulation = 500; % Number of simulation times (should be 100000)

% Define 3 cells to store arrays
Deg_bin = cell(num_simulation,1);
Deg_org = cell(num_simulation,1);
eigen_Q = cell(num_simulation,1);

%% Computation with 100000 random graph for both Degrees and Eigenvalues
% Store all data in the cell, including unique array, and the original value
for i = 1:1:num_simulation
    % Generate the ER graph
    A = erdos_reyni(N, p);
    Deg = A * u;
    % Plot distribution of degree and eigenvalues for one ER graph
    if i == 1 % for question 1, only one simulation is sufficient
        % Compute the Laplacian matrix and its eigenvalues
        sorted_Deg = sort(Deg);
        diag_matrix = diag(Deg);
        Q1 = diag_matrix - A;
        eigen_Q1 = eig(Q1);
        
        % Plot
        plot(sorted_Deg)
        hold on
        plot(eigen_Q1)
        xlabel('k')
        ylabel('Degrees and Laplacian eigenvalues')
        title('The degree vector and the Laplacian eigenvalues of the ER graph')
        legend('degrees d_{(k)}','Laplacian eigenvalues u_{(k)}')
        hold off
    end
    Diag_matrix = diag(Deg);
    Q = Diag_matrix - A;
    eigen_Q(i,1) = {eig(Q)};
    Deg_org(i,1) = {Deg};
    Deg_bin(i,1) = {unique(Deg)};
end

% Compute all unique value from all data, and hist data
Deg_all = unique(cell2mat(Deg_bin));
Deg_hist = hist(cell2mat(Deg_org), Deg_all);

%% Plots
figure
plot(Deg_all, Deg_hist/(N*num_simulation),'-o')
hold on

% Laplacian eigenvalues distribution
rounded_eigen_Q = round(cell2mat(eigen_Q));
rounded_eigen_Q_bin = unique(rounded_eigen_Q);
rounded_eigen_Q_hist = hist(rounded_eigen_Q, rounded_eigen_Q_bin)/(N*num_simulation);

plot(rounded_eigen_Q_bin, rounded_eigen_Q_hist,'-*' )
xlabel('x')
ylabel('f_?(x)')
title('The distribution of degrees and Laplacian eigenvalues (ER)')
legend('degrees','Laplacian eigenvalues')
hold off

% Fitting Laplacian eigenvalues
figure
plot(rounded_eigen_Q_bin, rounded_eigen_Q_hist,'r.','MarkerSize',25 )
hold on
pd = fitdist(rounded_eigen_Q,'Kernel'); % fitting use Kernel distribution
y = pdf(pd, rounded_eigen_Q_bin);
plot(rounded_eigen_Q_bin, y, 'LineWidth', 2) % fitting figure
xlabel('k')
ylabel('f_?(x)')
legend('Distribution','Fitting')
title('Fitting Laplacian eigenvalues distribution by Kernel function (ER)')
hold off
