%% Clear result of last computation
clear;
close all;
clc;

%% Parameters declaration
% Value assignment
N = 500;
d_av = 12;
p1 = d_av/(N-1);
m = 6;
m0 = 7; % m0 > m
k= 12;
p2 = 0.1;

%% Plot all the figures
graph_fig('ER', [N p1], 100);
graph_fig('BA', [N m0 m], 100);
graph_fig('WS', [N k p2], 100);