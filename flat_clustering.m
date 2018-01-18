function [delt,s_hat] = flat_clustering(S,parent)

    % get the parent-child node pair and preallocate indicator vectors
    uniqueParent = fliplr(unique(parent(parent>0)));
    s_hat = S;
    
    delt = true(1,numel(S));
    for clust = uniqueParent
        % sum stability of children of current parent
        children = (parent == clust);
        s_child = sum(s_hat(children));

        [s_hat(clust),ind] = max([S(clust),s_child]);

        if ind == 1
            % this parent node is stable
            children_all = subnodes(clust,parent);
            delt(children_all) = false;
        else
            delt(clust) = false;
        end
    end
    delt = find(delt);
end