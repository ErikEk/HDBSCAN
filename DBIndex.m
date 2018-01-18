function DB = DBIndex(X,idx,C)
numClusters = max(unique(idx));
N = length(X);
%Piecewise distance for each datapoint to clustering center
distCluster = zeros(1,numClusters);
numElements  = zeros(1,numClusters);
for i=1:numClusters
    for j=1:N
        if idx(j) == i
            distCluster(i) = distCluster(i) + sqrt( (C(i,1)-X(j,1))^2 + (C(i,2)-X(j,2))^2  );
            numElements(i) = numElements(i) +1;
        end
    end
end
aveDist = distCluster./numElements;

temp = zeros(1,numClusters);
sum = 0.0;
for i=1:numClusters
    for j=1:numClusters
        if i ~= j
            temp(j) =(aveDist(i) + aveDist(j)) / sqrt( (C(i,1)-C(j,1))^2 + (C(i,2)-C(j,2))^2  );
        end
    end
    sum = sum + max(temp);
    temp = zeros(1,numClusters);
end

DB = sum/double(numClusters);