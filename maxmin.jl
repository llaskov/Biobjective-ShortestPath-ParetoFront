# Compute a list of all (1, n)-paths with maximum capacity

# include
#include("network.jl")                       # Representation of network
#include("fibonacci-heap.jl")                # Fibonacci heap data structure

# Record to store in the Fibonacci heap
# for Dijkstra's algorithm for maximal capacity
mutable struct CpctRec

    # data fields
    id::Int64                               # the id of the vertex
    cpct_est::Float64                       # capacity estimate is used as 
                                            # the key value for the Fibonacci heap

    # Constructor
    # id            id of the vertex
    # cpct_est      shortest distance estimate
    # return        the constructed record
    function CpctRec(id::Int64, cpct_est::Float64)
        vr = new(id, cpct_est)

        return vr
    end
end

# overload less than operator for the Fibonacci heap key value
function Base.:<(vr_frst::CpctRec, vr_scnd::CpctRec)
    # we look for maximal capacity so reverse the inequality
    return vr_frst.cpct_est > vr_scnd.cpct_est
end

# Initialize the max capacity path estimates
# numb          number of vertices in the graph
# srce          source vertex
# return        cpct   vector of max capacity path estimates
function initcpct(numb::Int64, srce::Int64)
    cpct = Vector{Union{Nothing, Float64}}(nothing, numb)

    for i in eachindex(cpct)
        cpct[i] = -Inf
    end
    cpct[srce] = Inf

    return cpct
end

# Modification of the Dijkstra's algorithm that finds the capacity of the n-th vertex
# alnet         adjacency list of the input network
# return        capacity of the n-th vertex
function capacity(alnet::AdjLst{NtwrkRec})
    numb = size(alnet)                      # number of vertices in the network
    cpct = initcpct(numb, 1)                # vector of capacity estimates
    fheap = FibonacciHeap{CpctRec}()        # key: capacity
    table =                                 # table of heap nodes by vertex id
        Vector{Union{Nothing, FHeapNode{CpctRec}}}(nothing, numb)

    # fill in Fibonacci heap and table of nodes
    for i = 1:numb
        vr = CpctRec(i, cpct[i])
        table[i] = insert!(fheap, vr)
    end

    # main loop of the algorithm
    while !isempty(fheap)
        node = extractmin!(fheap)           # extract the minimum node from heap
        sta_vrtx = node.key.id              # the id of the current min vertex

        # visit each adjacent vertex of the current
        for adj in alnet.list[sta_vrtx]
            end_vrtx = adj.vrtx_end         # vertex id
            wght = adj.wght_two             # edge weight: capacity

            # relax procedure
            new_cpct =                      # min of capacity so far and edge weight
                min(cpct[sta_vrtx], wght)
            if cpct[end_vrtx] < new_cpct
                # update estimate in the vector
                cpct[end_vrtx] = new_cpct 

                # reorder the Fibonacci heap by decreasekey function
                new_key = CpctRec(end_vrtx, new_cpct)
                decreasekey!(fheap, table[end_vrtx], new_key)
            end
        end
    end

    return cpct[numb]
end

# Compose the maximum capacity subnetwork
# alnet         adjacency list of the input network
# return        n_cpct          capacity of the n-th vertex
#               al_max_cpcpt    adjacency list of the maximum capacity subnetwork
function maxmin(alnet::AdjLst{NtwrkRec})
    n_cpct = capacity(alnet)                # capacity of the n-th vertex
    al_max_cpcpt = AdjLst{NtwrkRec}(size(alnet)) 
    for i in eachindex(alnet.list), j in eachindex(alnet.list[i])
        if alnet.list[i][j].wght_two >= n_cpct
            pushadj!(
                al_max_cpcpt,                   # adjacency list 
                i,                              # starting vertex of the edge 
                alnet.list[i][j])               # record to insert in adjacency list
        end
    end

    return n_cpct, al_max_cpcpt
end
