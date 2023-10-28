# Compute the list of all Pareto optimal paths

# include
include("network.jl")                       # Representation of network
include("fibonacci-heap.jl")                # Fibonacci heap data structure
include("minsum.jl")                        # List of all paths with minimum length 
include("maxmin.jl")                        # List of all paths with maximum capacity

# Compose a subnetwork of those edges for which the capacity is g(i, j) > c
# al_in_net     adjacency list of the input network
# cpct_bnd      the capacity bound c
# return        adjacency list of the subnetwork with restricted capacity edges
function restrict(al_in_net::AdjLst{NtwrkRec}, cpct_bnd::Float64)
    # adjacency list of the subnetwork with restricted capacity edges
    al_rt_net = AdjLst{NtwrkRec}(size(al_in_net))
    for i in eachindex(al_in_net.list), j in eachindex(al_in_net.list[i])
        if al_in_net.list[i][j].wght_two > cpct_bnd
            pushadj!(
                al_rt_net,                  # adjacency list 
                i,                          # starting vertex of the edge 
                al_in_net.list[i][j])       # record to insert in adjacency list
        end
    end

    return al_rt_net
end

# Compose a list of Pareto optimal solutions
# al_in_net     adjacency list of the input network
# return        list of Pareto optimal solutions
function lpop(al_in_net::AdjLst{NtwrkRec})
    par_list =                              # resulting list 
        Vector{AdjLst{NtwrkRec}}(undef, 0)  
    al_rs_net = al_in_net                   # copy of the adj. list to be restricted
    n_cpct_0 = capacity(al_in_net)
    more = true
    while more
        dist, prnt = minsum(al_in_net, 1)
        al_in_net = outadj(al_in_net, prnt)
        if dist == Inf
            more = false
        else
            n_cpct_1, al_max_cpct = maxmin(al_in_net)
            push!(par_list, al_max_cpct)
            if n_cpct_0 == n_cpct_1
                more = false
            else
                al_rs_net = restrict(al_rs_net, n_cpct_1)
            end
            al_in_net = al_rs_net
        end
    end

    return par_list
end

# test minimum length
println("--> Test minimum length:")
al_in_net = readadj("net11.txt")            # input network adjacency list
println("Size (number of vertics): ", size(al_in_net))
dist, prnt = minsum(al_in_net, 1)
display(prnt)
paths_lst = listpaths(prnt)
display(paths_lst)
spnet = outadj(al_in_net, prnt)
display(spnet)

# test maximum capacity
println("--> Test maximum capacity:")
n_cpct, al_max_cpct = maxmin(al_in_net)
println("Capacity of n: ", n_cpct)
display(al_max_cpct)

# test list of Pareto optimal paths
println("--> Test Pareto optimal paths:")
par_list = lpop(al_in_net)
display(par_list)
