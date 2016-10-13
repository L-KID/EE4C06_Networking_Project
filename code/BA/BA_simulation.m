% This file generates a BA graph and does some subsequent computation.

% Clear result of last computation
clear;
close all;
clc;

%% Parameters declaration
% Value assignment
N = 500;
m = 6;
m0 = 7; % m0 > m
u = ones(500, 1);
num_simulation = 100; % Number of simulation times (should be 100000)

% Define 3 cells to store arrays
Deg_bin = cell(num_simulation,1);
Deg_org = cell(num_simulation,1);
eigen_Q = cell(num_simulation,1);

%% Plot distribution of degree and eigenvalues for 1 BA graph
A = scalefree(N, m0, m);
Deg = A * u;
sorted_Deg = sort(Deg);
diag_matrix = diag(Deg);
Q1 = diag_matrix - A;
eigen_Q1 = eig(Q1);
% Plot 
plot(sorted_Deg)
hold on
plot(eigen_Q1)
xlabel('k')
ylabel('Degrees and the Laplacian eigenvalues')
title('The degree vector and the Laplacian eigenvalues of the BA graph')
legend('degrees d_{(k)}','Laplacian eigenvalues u_{(k)}')
hold off

%% Computation with 100000 random graph for both Degree and Eigenvalue
% Store all data in the cell, including unique array, and the original value
for i = 1:1:num_simulation
    % Generate the BA graph
    A = scalefree(N, m0, m);
    Deg = A * u;
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
loglog(Deg_all, Deg_hist/(N*num_simulation))
xlim([0 200]) % limit x axis
hold on

% Laplacian eigenvalues distribution
rounded_eigen_Q = round(cell2mat(eigen_Q));
rounded_eigen_Q_bin = unique(rounded_eigen_Q);
rounded_eigen_Q_hist = hist(rounded_eigen_Q, rounded_eigen_Q_bin)/(N*num_simulation);

loglog(rounded_eigen_Q_bin, rounded_eigen_Q_hist )
xlim([0 200])
xlabel('x')
ylabel('f_?x)')
title('The distribution of degrees and Laplacian eigenvalues (BA)')
legend('degrees','Laplacian eigenvalues')
hold off

% Fitting Laplacian eigenvalues
figure
loglog(rounded_eigen_Q_bin, rounded_eigen_Q_hist,'r.','MarkerSize',25 )
hold on
pd = fitdist(rounded_eigen_Q,'Kernel'); % fitting use Kernel distribution
y = pdf(pd, rounded_eigen_Q_bin);
loglog(rounded_eigen_Q_bin, y, 'LineWidth', 2) % fitting figure
xlabel('k')
ylabel('f_?x)')
legend('Distribution','Fitting')
title('Fitting Laplacian eigenvalues distribution by Kernel function (BA)')
hold off


% Fitting by Gaussian distribution
% For BA, maybe it is not a Gaussian distribution, maybe a exp 
% f = fittype('a*exp(-((x-b)/c)^2)');
% loglog(rounded_eigen_Q_bin, rounded_eigen_Q_hist,'r.','MarkerSize',25 )
% startPoints = [0.15 12 3];
% [cfun, gof] = fit(rounded_eigen_Q_bin(:), rounded_eigen_Q_hist(:), f, 'Start', startPoints);
% yy = cfun.a*exp(-((rounded_eigen_Q_bin-cfun.b)/cfun.c).^2);
% hold on
% loglog(rounded_eigen_Q_bin, yy,'LineWidth', 2)
% xlabel('k')
% ylabel('f_?x)')
% legend('Distribution','Fitting')
% title('Fitting Laplacian eigenvalues distribution by Kernel function (BA)')
% hold off