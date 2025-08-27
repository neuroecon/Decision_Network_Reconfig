function Ci_modified = revise_community(Ci, threshold)
Ci_modified = Ci;  % make a copy to modify
unique_communities = unique(Ci);

for i = 1:length(unique_communities)
    comm = unique_communities(i);
    nodes_in_comm = find(Ci == comm);
    if numel(nodes_in_comm) <= threshold
        Ci_modified(nodes_in_comm) = 100;
    end
end