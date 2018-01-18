function [scatteringCenters,num_clusters] = HDBSCANClustering(X, Y, SNR, min_pts,min_clust_size,outlier_thresh)
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
    children = subnodes(thisClust,tree.parents);
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

% END HDBSCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numSC = max(idx);
num_clusters = numSC;
% Find out which scattering center has higher average intensity and find out if the distance between is too close.


%% Create cell with cluster results

for scInd = 1:numSC
    scatteringCenters.scX(scInd) = C(scInd,2);
    scatteringCenters.scY(scInd) = C(scInd,1);
    scatteringCenters.scXstddev(scInd) = std(X(idx==scInd,2));
    scatteringCenters.scYstddev(scInd) = std(X(idx==scInd,1));
    scatteringCenters.scSNR(scInd) = mean(SNR(idx==scInd));
    scatteringCenters.scSNRstddev(scInd) = std(SNR(idx==scInd));
    
    clusterColor{scInd} = cmap(scInd*4,:);
    
    % Hardcode the first 3 colors as distinguishable for debugging purposes.
    clusterColor{1} = [102 152 255]/255;
    clusterColor{2} = [140 255 102]/255;
    clusterColor{3} = [209 122 189]/255;
    clusterColor{4} = [50 122 20]/255;
    clusterColor{5} = [50 2 200]/255;
    
    scatter(X(idx==scInd,1),X(idx==scInd,2),20,clusterColor{scInd},'filled')
    error_ellipse(X(idx==scInd,:))
    
    
    plot(C(scInd,1),C(scInd,2),'rx','MarkerSize',18,'LineWidth',3) % Plot calculated centerpoint of the clusters.
end
scatter(X(idx==0,1),X(idx==0,2),20,[0.5 0.5 0.5],'filled')


title 'HDBSCAN clustering'
axis([-4 4 -4 4])
xlabel('Lateral position [m]')
ylabel('Relative range [m]')
set(gca,'YDir','normal');

end