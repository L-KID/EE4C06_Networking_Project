% This file generates an ER graph and does some subsequent computation.

%% Clear result of last computation
clear;
close all;
clc;

%% Parameters declaration
% Value assignment
N = 500;
d_av = 12;
p = d_av/(N-1);
u = ones(500, 1);
NUM_SIMULATION = 100000; % Number of simulation times (should be 100000)

% Define 3 cells to store arrays
Deg_bin = cell(NUM_SIMULATION,1);
Deg_org = cell(NUM_SIMULATION,1);
eigen_Q = cell(NUM_SIMULATION,1);

%% Plot degrees and eigenvalues for one ER graph
A = erdos_reyni(N, p);
Deg = A * u;

% Compute the Laplacian matrix and its eigenvalues
sorted_Deg = sort(Deg);
diag_matrix = diag(Deg);
Q = diag_matrix - A;
eig_Q = eig(Q);
        
% Plot
plot(sorted_Deg)
hold on
plot(eig_Q)
ylim([0,inf])
xlabel('k')
ylabel('Degrees and the Laplacian eigenvalues')
title('The degree vector and the Laplacian eigenvalues of one ER graph')
legend('degrees d_{(k)}','Laplacian eigenvalues u_{(k)}')
hold off
savefig('../../figures/ER/fig/ER_deg_and_eig_single.fig');
saveas(gcf, '../../figures/ER/png/ER_deg_and_eig_single.png');

%% Plot degrees, eigenvalues and their distribution for multiple ER graphs
total_eigen = zeros(500,1);
total_Deg = zeros(500,1);
for i = 1:1:NUM_SIMULATION
    % Generate the ER graph
    A = erdos_reyni(N, p);
    Deg = A * u;
    Diag_matrix = diag(Deg);
    Q = Diag_matrix - A;    
    eig_Q = eig(Q);
    sorted_Deg = sort(Deg);
    total_Deg = total_Deg + sorted_Deg;
    total_eigen = total_eigen + eig_Q;  
    % Store all data in the cell, including unique array, and the original value
    eigen_Q(i,1) = {eig(Q)};
    Deg_org(i,1) = {Deg};
    Deg_bin(i,1) = {unique(Deg)};
end

%Plot degrees and eigenvalue for multiple ER graphs
figure
plot(total_Deg/NUM_SIMULATION)
hold on
plot(total_eigen/NUM_SIMULATION)
ylim([0,inf])
xlabel('k')
ylabel('Degrees and the Laplacian eigenvalues')
title('The average degrees and Laplacian eigenvalues of multiple ER graphs')
legend('degrees d_{(k)}','Laplacian eigenvalues u_{(k)}')
hold off
savefig('../../figures/ER/fig/ER_deg_and_eig_multiple.fig');
saveas(gcf, '../../figures/ER/png/ER_deg_and_eig_multiple.png');

% Compute all unique value from all data, and hist data
Deg_all = unique(cell2mat(Deg_bin));
Deg_hist = hist(cell2mat(Deg_org), Deg_all);

% Plot distribution of degree and eigenvalues for the ER graph
figure
plot(Deg_all, Deg_hist/(N*NUM_SIMULATION),'-o')
hold on

% Laplacian eigenvalues distribution
rounded_eigen_Q = round(cell2mat(eigen_Q));
rounded_eigen_Q_bin = unique(rounded_eigen_Q);
rounded_eigen_Q_hist = hist(rounded_eigen_Q, rounded_eigen_Q_bin)/(N*NUM_SIMULATION);

plot(rounded_eigen_Q_bin, rounded_eigen_Q_hist,'-*' )
xlabel('x')
ylabel('f_\mu(x)')
title('The distribution of degrees and Laplacian eigenvalues (ER)')
legend('degrees','Laplacian eigenvalues')
hold off
savefig('../../figures/ER/fig/ER_distribution.fig');
saveas(gcf, '../../figures/ER/png/ER_distribution.png');

%% Fitting distribution of Laplacian eigenvalues
figure
plot(rounded_eigen_Q_bin, rounded_eigen_Q_hist,'r.','MarkerSize',25 )
hold on
% pd = fitdist(rounded_eigen_Q,''); % fitting use Kernel distribution
% y = pdf(pd, rounded_eigen_Q_bin);
[f, xi] = ksdensity(rounded_eigen_Q,'Width',0.8);
plot(xi, f, 'LineWidth', 2)
%plot(rounded_eigen_Q_bin, y, 'LineWidth', 2) % fitting figure
xlabel('k')
ylabel('f_\mu(x)')
legend('Distribution','Fitting')
title('Fitting Laplacian eigenvalues distribution by Kernel function (ER)')
hold off
savefig('../../figures/ER/fig/ER_fitting.fig');
saveas(gcf, '../../figures/ER/png/ER_fitting.png');
