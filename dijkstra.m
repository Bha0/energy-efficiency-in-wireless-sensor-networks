% Define a function to find the shortest path using Dijkstra's algorithm
function shortestPath = dijkstra(graph, startNode, endNode)
    numNodes = size(graph, 1);
    visited = false(1, numNodes);
    dist = Inf(1, numNodes);
    prev = zeros(1, numNodes);
    
    dist(startNode) = 0;
    
    for i = 1:numNodes
        % Find the unvisited node with the shortest distance
        u = -1;
        minDist = Inf;
        for j = 1:numNodes
            if ~visited(j) && dist(j) < minDist
                u = j;
                minDist = dist(j);
            end
        end
        
        if u == -1
            break;
        end
        
        visited(u) = true;
        
        % Update distances to neighbors
        for v = 1:numNodes
            if ~visited(v) && graph(u, v) > 0
                alt = dist(u) + graph(u, v);
                if alt < dist(v)
                    dist(v) = alt;
                    prev(v) = u;
                end
            end
        end
    end
    
    % Reconstruct the shortest path
    u = endNode;
    shortestPath = [];
    while u ~= startNode
        shortestPath = [u, shortestPath];
        u = prev(u);
    end
    shortestPath = [startNode, shortestPath];
end