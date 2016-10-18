% This function is to generate figures of degree, eigenvalue and their 
% distribution for different graphs

function y = graph_fig(G, P, N)
% G is the type of graph, which should be 'ER', 'BA', or 'WS'
% P is a array for the parameters of graph
% N is the number of simulation

u = ones(500, 1);
% Define 3 cells to store arrays
Deg_bin = cell(N,1);
Deg_org = cell(N,1);
eigen_Q = cell(N,1);

total_eigen = zeros(P(1),1);
total_Deg = zeros(P(1),1);

%% Compute degrees and eigenvalues for N times and store in cells
if G=='ER'
    for i = 1:1:N
            % Generate the ER graph
            A = erdos_reyni(P(1), P(2));
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
            
elseif G=='BA'
            for i = 1:1:N
            % Generate the BA graph
            A = scalefree(P(1),P(2),P(3));
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
    
elseif G=='WS'
    for i = 1:1:N
            % Generate the WS graph
            A = small_world(P(1),P(2),P(3));
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
    
else
        disp('Wrong graph type');
        return 
end

%% Plot degrees and eigenvalues for one graph
Deg = A * u;
sorted_Deg = sort(Deg);
diag_matrix = diag(Deg);
Q = diag_matrix - A;
eig_Q = eig(Q);
% Plot 
figure
plot(sorted_Deg)
hold on
plot(eig_Q)
ylim([0,inf])
xlabel('k')
ylabel('Degrees and the Laplacian eigenvalues')
title(['The degree vector and the Laplacian eigenvalues of one ' G ' graph'])
legend('degrees d_{(k)}','Laplacian eigenvalues u_{(k)}')
hold off
savefig(['../../figures/' G '/fig/' G '_deg_and_eig_single.fig']);
saveas(gcf, ['../../figures/' G '/png/' G '_deg_and_eig_single.png']);

%% Plot degrees and eigenvalues for multiple graphs
figure
plot(total_Deg/N)
hold on
plot(total_eigen/N)
ylim([0,inf])
xlabel('k')
ylabel('Degrees and the Laplacian eigenvalues')
title(['The average degrees and Laplacian eigenvalues of multiple ' G ' graphs'])
legend('degrees d_{(k)}','Laplacian eigenvalues u_{(k)}')
hold off
savefig(['../../figures/' G '/fig/' G '_deg_and_eig_multiple.fig']);
saveas(gcf, ['../../figures/' G '/png/' G '_deg_and_eig_multiple.png']);
              
%% Compute all unique value from all data, and hist data
Deg_all = unique(cell2mat(Deg_bin));
Deg_hist = hist(cell2mat(Deg_org), Deg_all);

% Laplacian eigenvalues distribution
rounded_eigen_Q = round(cell2mat(eigen_Q));
rounded_eigen_Q_bin = unique(rounded_eigen_Q);
rounded_eigen_Q_hist = hist(rounded_eigen_Q, rounded_eigen_Q_bin)/(P(1)*N);

%% Plot distribution of degree and eigenvalues
if G=='BA'
figure
loglog(Deg_all, Deg_hist/(P(1)*N))
hold on
loglog(rounded_eigen_Q_bin, rounded_eigen_Q_hist )
xlabel('x')
ylabel('f_\mu(x)')   
title('The distribution of degrees and Laplacian eigenvalues (BA)')
legend('degrees','Laplacian eigenvalues')
hold off
savefig('../../figures/BA/fig/BA_distribution.fig');
saveas(gcf, '../../figures/BA/png/BA_distribution.png');
else
figure
plot(Deg_all, Deg_hist/(P(1)*N),'-o')
hold on
plot(rounded_eigen_Q_bin, rounded_eigen_Q_hist,'-*' )
xlabel('x')
ylabel('f_\mu(x)')
title(['The distribution of degrees and Laplacian eigenvalues (' G ')'])
legend('degrees','Laplacian eigenvalues')
hold off
savefig(['../../figures/' G '/fig/' G '_distribution.fig']);
saveas(gcf, ['../../figures/' G '/png/' G '_distribution.png']);
end

%% Fitting
if G=='ER'
figure
plot(rounded_eigen_Q_bin, rounded_eigen_Q_hist,'r.','MarkerSize',25 )
hold on
[f, xi] = ksdensity(rounded_eigen_Q,'Width',0.8);
plot(xi, f, 'LineWidth', 2)
xlabel('k')
ylabel('f_\mu(x)')
legend('Distribution','Fitting')
title('Fitting Laplacian eigenvalues distribution by Kernel function (ER)')
hold off
savefig('../../figures/ER/fig/ER_fitting.fig');
saveas(gcf, '../../figures/ER/png/ER_fitting.png');
end

if G=='BA'
figure
loglog(rounded_eigen_Q_bin, rounded_eigen_Q_hist,'r.','MarkerSize',25 )
hold on
pd = fitdist(rounded_eigen_Q,'Kernel'); % fitting use Kernel distribution
y = pdf(pd, rounded_eigen_Q_bin);
loglog(rounded_eigen_Q_bin, y, 'LineWidth', 2) % fitting figure
xlabel('k')
ylabel('f_\mu(x)')
legend('Distribution','Fitting')
title('Fitting Laplacian eigenvalues distribution by Kernel function (BA)')
hold off
savefig('../../figures/BA/fig/BA_fitting.fig');
saveas(gcf, '../../figures/BA/png/BA_fitting.png');
end

if G=='WS'
figure
f = fittype('a*exp(-((x-b)/c)^2)');
plot(rounded_eigen_Q_bin, rounded_eigen_Q_hist,'r.','MarkerSize',25 )
startPoints = [0.15 12 3];
[cfun, gof] = fit(rounded_eigen_Q_bin(:), rounded_eigen_Q_hist(:), f, 'Start', startPoints);
yy = cfun.a*exp(-((rounded_eigen_Q_bin-cfun.b)/cfun.c).^2);
hold on
plot(rounded_eigen_Q_bin, yy,'LineWidth', 2)
xlabel('k')
ylabel('f_\mu(x)')
legend('Distribution','Fitting')
title('Fitting Laplacian eigenvalues distribution by Kernel function (WS)')
hold off
savefig('../../figures/WS/fig/WS_fitting.fig');
saveas(gcf, '../../figures/WS/png/WS_fitting.png');
end

end


