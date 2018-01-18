function res = fit_hdbscan(X,min_points, min_clust_size)
    % inputs:
    % X (observations x dimensions)
    % min_points - min points
    % min_clust_size - minimum number of points in a cluster

    %% create connected tree 
    % compute core distances and the mutual reachability
    d = sum(X.^2,2); 
    D = real(sqrt(bsxfun(@plus,d,d')-2*(X*X')));

    dCore = sort(D,1);
    dCore = dCore(min_points,:);

    % mutual reachability
    [n,~] = size(X);
    mr = D;
    maxCore = bsxfun(@max,dCore,dCore');
    idx = maxCore > mr;
    mr(idx) = maxCore(idx);
    
    
    % create a minimum spanning tree (matlab build in)
    mst = minspantree(graph(mr),'Method','dense','Type','tree');
    
    nodes = mst.Edges.EndNodes;
    weights = mst.Edges.Weight;

    % sort weights (descended)
    epsilon = sort(unique( weights),'descend');
    epsilon = epsilon(1:1:end);
    num_epsilon = numel(epsilon);

    noise_lambda = zeros(n,1);
    cluster_last = zeros(n,1);
    parent_clust = zeros(1,100); % 100 is max clusters
    lambda_mininum = zeros(1,100);
    maximum_lambda = zeros(n,100);
    current_max = 1;
    lambda_mininum(1) =1./epsilon(1);
    new_ide = ones(n,1,'int8');

    %% hierarchical search
    for i = 2:num_epsilon
        
        old_ide = new_ide;
        
        % find edges greater than current epsilon value
        idx = weights > epsilon(i);
        if ~any( idx )
            continue
        end
        
        % find the nodes of the cut edges and remove bad ones (those previously labeled as noise)
        nodes_end = nodes(idx,:);
        nodes(idx,:) = [];
        weights(idx) = []; 
        mst = mst.rmedge( nodes_end(:,1),nodes_end(:,2) );

        
        % remove noise
        self_cuts = (nodes_end(:,1) == nodes_end(:,2));
        if any( self_cuts )
            new_ide(nodes_end(self_cuts,1)) = 0;
            nodes_end(self_cuts,:) = [];
            if isempty( nodes_end )
                continue
            end
        end
        
        % remove noisy nodes and skip loop if no remaining nodes 
        uniqueCuts = unique( nodes_end );
        idx_bad = (new_ide(uniqueCuts) == 0);     
        if any( idx_bad )
            if all( idx_bad )
                continue
            end
            
            % if only some nodes noisy, remove those
            nodes_bad = uniqueCuts(idx_bad);
            nodes_end(any(ismember(nodes_end',nodes_bad))',:) = []; % remove
        end

        
        %% find sub-trees    
        for k = 1:size(nodes_end,1)

            % get the connected components from the end nodes
            parent = old_ide(nodes_end(k,1));
            sub_t = bfsearch(mst,nodes_end(k,1));
            sub_t2 = bfsearch(mst,nodes_end(k,2));
            n_t1 = length(sub_t);
            n_t2 = length(sub_t2);
            val_tree = [n_t1,n_t2] >= min_clust_size;
            
            %  check for noise
            if val_tree
                new_max = current_max + 2;
                tem = new_max - 1;
                new_ide(sub_t) = tem;
                new_ide(sub_t2) = new_max; 
                parent_clust(tem:new_max) = parent;
                current_max = new_max;
                maximum_lambda([sub_t;sub_t2],parent) = epsilon(i-1);
                lambda_mininum([tem,new_max]) = epsilon(i);
            else % or split
                if val_tree(1)   
                    is_noises = sub_t2;
                elseif val_tree(2)
                    is_noises = sub_t;
                else
                    is_noises = [sub_t;sub_t2];
                end
                
                cluster_last(is_noises) = old_ide(is_noises);
                noise_lambda(is_noises) = epsilon(i-1);
                maximum_lambda(is_noises,parent) = epsilon(i-1);
                new_ide(is_noises) = 0;
                
                if ~any( new_ide )
                    break
                end
            end
        end
    end

    % calculate the condensed tree
    lambda_mininum = 1./lambda_mininum(1:current_max);
    maximum_lambda = 1./maximum_lambda(:,1:current_max);
    maximum_lambda(isinf(maximum_lambda))=0;
    parent_clust = parent_clust(1:current_max);

    % compute stabilities for the clusters
    in_cluster = ones( n,current_max );
    in_cluster(maximum_lambda == 0) = 0; 
    S = sum( maximum_lambda - bsxfun( @times,lambda_mininum,in_cluster ) );
    clusterTree = struct( 'clusters',1:current_max,'parents',parent_clust,'lambdaMin',lambda_mininum,'stability',S );    
    res = struct( 'lambda',1./epsilon(1:i),'clusterTree',clusterTree,'dCore',dCore,'lambdaNoise',1./noise_lambda,'lastClust',cluster_last,'lambdaMax',sparse( maximum_lambda ) );
end