# Biobjective-ShortestPath-ParetoFront
Find all Pareto optimal solutions of the biobjective shortest path problem in a network.

A network with two objective functions is considered:
1. A linear function that denotes the minimum distance.
2. A bottleneck function that denotes the maximum capacity.

The algorithm finds the complete list of all Pareto optimal solutions of the biobjective shortest path problem in the network. The method is based on two modifications of the Dijkstra's algorithm, and a general procedure that calculates all the elements of the Pareto front. The modifications of the Dijkstra's algorithm uses a Fibonacci heap for a priority queue operations that ensures $O(n log n + m)$ complexity, where $n$ is the number of vertices of the network and $m$ is the number of edges. The entire method has complexity $k_0 O(n log n + m)$, where $k_0$ is the number of classes of Pareto equivalent paths.
