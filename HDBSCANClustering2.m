function [idx,C] = HDBSCANClustering2(X, Y, min_pts,min_clust_size,outlier_thresh)
% Colormap for random cluster color.
cmap = colormap(colorcube); 
X = [Y; X]';

% OLD KMEANS
%[idx,C] = kmeans(X,numSC,'Distance','cityblock','Replicates',5,'Options',opts);

% NEW (HDBSCAN) %%%%%%%%%%%%%%%%%%%%%
% the hierarchical cluster structure

res = fit_hdbscan(X,min_pts,min_clust_size);


% extract best clusters
tree = res.clusterTree;
% Get The best / optimal clustering (max lambdas)
best_clusters = flat_clustering( tree.stability,tree.parents );
lambdaMax = full(res.lambdaMax);
p = length( best_clusters );
corePoints = cell( 1,p );
coreLambda = zeros( 1,p );
for ii = 1:numel(best_clusters)
    thisClust = best_clusters(ii);
    % find subnodes
    children = get_all_subnodes(thisClust,tree.parents);
    % get maximum lambda value of this cluster or its children
    maxLambda = max( lambdaMax(:,[thisClust,children]),[],2 );
    coreLambda(ii) = max(maxLambda);
    corePoints{ii} = find( maxLambda == coreLambda(ii) );
end



% compute the outlier scores
tree = res.clusterTree;
score = outliner_score( best_clusters,tree.parents,res.lastClust,res.lambdaNoise );
% compute labels and probability of cluster membership
lambda_max = full(res.lambdaMax);
n = size(lambda_max,1);
P = zeros(n,1);
idx = zeros( n,1,'int8' );
% loop clusters
for k = 1:numel( best_clusters )
    pts = lambda_max(:,best_clusters(k)) > 0;
    idx(pts) = k;
    P(pts) = lambda_max(pts,best_clusters(k)) / coreLambda(k);
end



% set labels with outlier scores to zero
idx( score > outlier_thresh ) = 0;
% Calculate mid points
C = [];
for i = 1:max(idx)
    C(i,:) = mean(X(idx==i,:));
end



end
