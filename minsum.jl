# Compute a list of all (1, n)-paths with minimum length

# include
#include("network.jl")                       # Representation of network
#include("fibonacci-heap.jl")                # Fibonacci heap data structure

# Record to store in the Fibonacci heap
# for Dijkstra's algorithm for shortest path
mutable struct DistRec

    # data fields
    id::Int64                               # the id of the vertex
    dist_est::Float64                       # shortest distance estimate is used as 
                                            # the key value for the Fibonacci heap

    # Constructor
    # id            id of the vertex
    # dist          shortest distance estimate
    # return        the constructed record
    function DistRec(id::Int64, dist_est::Float64)
        vr = new(id, dist_est)

        return vr
    end
end

# overload less than operator for the Fibonacci heap key value
function Base.:<(vr_frst::DistRec, vr_scnd::DistRec)
    return vr_frst.dist_est < vr_scnd.dist_est
end

# Initialize the shortest path estimates and predecessors vectors
# numb          number of vertices in the graph
# srce          source vertex
# return        dist   vector of shortest path estimate
#               prnt   shortest path subnetwork predecessors adjacency list
function initdist(numb::Int64, srce::Int64)
    dist = Vector{Union{Nothing, Float64}}(nothing, numb)
    prnt = AdjLst{Int64}(numb)

    for i in eachindex(dist)
        dist[i] = Inf
    end
    dist[srce] = 0

    return dist, prnt
end

# Modification of the Dijkstra's algorithm that solves the MINSUM problem
# alnet         adjacency list of the input network
# srce          source vertex
function minsum(alnet::AdjLst{NtwrkRec}, srce::Int64)
    numb = size(alnet)                      # number of vertices in the network
    dist, prnt = initdist(numb, srce)    # initialize SP estimate, SP predecessor
    fheap = FibonacciHeap{DistRec}()        # key: shortest path estimate
    table =                                 # table of heap nodes by vertex id
        Vector{Union{Nothing, FHeapNode{DistRec}}}(nothing, numb)

    # fill in Fibonacci heap and table of nodes
    for i = 1:numb
        vr = DistRec(i, dist[i])
        table[i] = insert!(fheap, vr)
    end

    # main loop of the algorithm
    while !isempty(fheap)
        node = extractmin!(fheap)           # extract the minimum node from heap
        sta_vrtx = node.key.id              # the id of the current min vertex

        # visit each adjacent vertex of the current
        for adj in alnet.list[sta_vrtx]
            end_vrtx = adj.vrtx_end         # vertex id
            wght = adj.wght_one             # edge weight: distance

            # relax procedure
            new_dist =                      # sum of distance so far and edge weight 
                dist[sta_vrtx] + wght
            if dist[end_vrtx] > new_dist
                # update estimate in the vector
                dist[end_vrtx] = new_dist 

                # reorder the Fibonacci heap by decreasekey function
                new_key = DistRec(end_vrtx, new_dist)
                decreasekey!(fheap, table[end_vrtx], new_key)

                # remove the previous adjacency list and set new adjacent vertex
                newadj!(prnt, end_vrtx, sta_vrtx)
            elseif dist[end_vrtx] == new_dist
                # append the new adjacent vertex, and preserve the previous
                pushadj!(prnt, end_vrtx, sta_vrtx)
            end
        end
    end

    return dist, prnt
end

# Create a list of all (1, n)-paths from the shortest paths predecessor adjacency list
# following a BFS algorithmic scheme
# prnt          shortest path subnetwork predecessors adjacency list
# return        the list of paths
function listpaths(prnt::AdjLst{Int64})
    paths_lst = AdjLst{Int64}()         # list that contains the paths
    pushlst!(paths_lst)                 # crate the first path

    # inset the last vertex as first vertex in the first path
    prependadj!(paths_lst, 1, size(prnt))

    more = true
    while more

        more = false
        numb_paths = size(paths_lst)    # number of paths in the paths list

        # visit each path in the path list
        for indx_paths = 1:numb_paths

            # current is the last vertex inserted in the path
            curr = paths_lst.list[indx_paths][1]

            # get the predecessor of the current and prepend it the path
            if curr != 1
                prependadj!(paths_lst, indx_paths, prnt.list[curr][1])

                if prnt.list[curr][1] != 1
                    more = true
                end
            end

            # if the adj list of predecessors contains more than one predecessor,
            # copy the list or each of them and prepend them
            for indx_pred = 2:size(prnt.list[curr], 1)
                # copy the list
                pushlst!(paths_lst)                 # crate the path
                for i = 2:size(paths_lst.list[indx_paths], 1)
                    pushadj!(
                             paths_lst, 
                             size(paths_lst), 
                             paths_lst.list[indx_paths][i])
                end
                prependadj!(paths_lst, size(paths_lst), prnt.list[curr][indx_pred])

                if prnt.list[curr][indx_pred] != 1
                    more = true
                end
            end
        end
    end

    return paths_lst
end

# Represent the shortest paths subnetwork with the outgoing adjacency list
# alnet         adjacency list of the input network
# prnt          (1, n)-shortest path subnetwork predecessors adjacency list
# return        shortest paths subnetwork adjacency list
function outadj(alnet::AdjLst{NtwrkRec}, prnt::AdjLst{Int64})
    numb = size(alnet)                      # number of vertices
    spnet = AdjLst{NtwrkRec}(numb)          # shortest paths subnetwork adjacency list
    adj_mtrx = toadjmtrx(alnet)             # network as adjacency matrix

    for i = 2:numb
        for j = 1:size(prnt.list[i], 1)
            nr = NtwrkRec(
                i,                          # ending edge of the vertex
                # weight value of the first criteria
                adj_mtrx[prnt.list[i][j], i].first,
                # weight value of the second criteria
                adj_mtrx[prnt.list[i][j], i].second) 
            pushadj!(
                spnet,                      # adjacency list 
                prnt.list[i][j],            # starting vertex of the edge 
                nr)
        end
    end

    return spnet
end
