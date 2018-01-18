function score = outliner_score(clusters,parent,last_cluster,lambda_noise)
    % Calculate the outliner scores
    score = zeros(size(last_cluster));
    has_processed = false(1,max(clusters));
    for clust = clusters
        if has_processed(clust)
            continue
        end
        children = subnodes(clust,parent);
        pts = ismember(last_cluster,[clust,children]);
        r_lambda = max(lambda_noise(pts));
        score(pts) = 1 - (lambda_noise(pts) ./ r_lambda);
        has_processed([clust,children]) = true;
    end
end