% This function generates the adjacency matrix of the given graph

function A = graph_gen(G, P)
addpath(genpath('./'));
switch G
    case 'BA'
        A = scalefree(P(1), P(2), P(3));
    case 'ER'
        A = erdos_reyni(P(1), P(2));
    case 'WS'
        A = small_world(P(1), P(2), P(3));
    otherwise
        disp('Wrong graph type');
        return
end
end
