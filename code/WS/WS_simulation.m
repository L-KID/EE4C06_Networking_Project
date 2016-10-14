% This file generates a WS graph and does some subsequent computation.

% Clear result of last computation
clear;
close all;
clc;
%%
% Value assignment
N = 500;
k= 12;
p = 0.1 ;

u = ones(500, 1);
num_simulation = 2000;
% Define 3 cells to store arrays
Deg_bin = cell(num_simulation,1);
Deg_org = cell(num_simulation,1);
eigen_Q = cell(num_simulation,1);

%% Plot distribution of degree and eigenvalues for 1 WS graph
A = small_world(N, k, p);
Deg = A * u;
   
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
ylabel('Degrees and the Laplacian eigenvalues')
title('The degree vector and the Laplacian eigenvalues of the WS graph')
legend('degrees d_{(k)}','Laplacian eigenvalues u_{(k)}')
hold off


%% Computation with 100000 random graph for both Degree and Eigenvalue
% Store all data in the cell, including unique array, and the original value
for i = 1:1:num_simulation
    % Generate the WS graph
    A = small_world(N, k, p);
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
plot(Deg_all, Deg_hist/(N*num_simulation),'-o')
hold on;

% Laplacian eigenvalues distribution
rounded_eigen_Q = round(cell2mat(eigen_Q));
rounded_eigen_Q_bin = unique(rounded_eigen_Q);
rounded_eigen_Q_hist = hist(rounded_eigen_Q, rounded_eigen_Q_bin)/(N*num_simulation);

plot(rounded_eigen_Q_bin, rounded_eigen_Q_hist,'-*' )
xlabel('x')
ylabel('f_u(x)')
title('The distribution of degrees and Laplacian eigenvalues (WS)')
legend('degrees','Laplacian eigenvalues')
hold off

% Fitting Laplacian eigenvalues by Kernel
% Not so good
% figure
% plot(rounded_eigen_Q_bin, rounded_eigen_Q_hist,'r.','MarkerSize',25 )
% hold on
% pd = fitdist(rounded_eigen_Q,'Beta'); % fitting use Kernel distribution
% y = pdf(pd, rounded_eigen_Q_bin);
% plot(rounded_eigen_Q_bin, y, 'LineWidth', 2) % fitting figure
% xlabel('k')
% ylabel('f_?x)')
% legend('Distribution','Fitting')
% title('Fitting Laplacian eigenvalues distribution by Kernel function (WS)')
% hold off

% Fitting by Gaussian distribution
figure
f = fittype('a*exp(-((x-b)/c)^2)');
plot(rounded_eigen_Q_bin, rounded_eigen_Q_hist,'r.','MarkerSize',25 )
startPoints = [0.15 12 3];
[cfun, gof] = fit(rounded_eigen_Q_bin(:), rounded_eigen_Q_hist(:), f, 'Start', startPoints);
yy = cfun.a*exp(-((rounded_eigen_Q_bin-cfun.b)/cfun.c).^2);
hold on
plot(rounded_eigen_Q_bin, yy,'LineWidth', 2)
xlabel('k')
ylabel('f_u(x)')
legend('Distribution','Fitting')
title('Fitting Laplacian eigenvalues distribution by Kernel function (WS)')
hold off