function [minPts,minClusterSize,minThresholdSize] = BestHDBSCANTuning(x,y)

currentBestDBI = 1;

%epsilon = 0.1:0.1:1.5;
minPts = 3:1:6;
minClusterSize = 3:1:6;
tresholdSize = 0.7:0.1:1.0;
%list = zeros(length(epsilon),length(minPts));

for i=1:length(minPts)
    for j=1:length(minClusterSize)
        for k=1:length(tresholdSize)
            [idx,C] = HDBSCANClustering2(x, y, minPts(i), minClusterSize(j),tresholdSize(k));
            data = [y;x]';
            DB = DBIndex(data,idx,C);

            if DB == 0
                DB = 1;
            end
            %list(i,j) = DB;  %For plotting
            if DB < currentBestDBI
                currentBestDBI = DB;
                currentBestMinPts = minPts(i);
                currentBestMinClusterSize = minClusterSize(j);
                currentBestThresholdSize = tresholdSize(k);
            end
        end
    end
end

% surf(list)
% xlabel('minpts')
% ylabel('eps')
% zlabel('DBI')



%Output
minPts = currentBestMinPts;
minClusterSize = currentBestMinClusterSize;
minThresholdSize = currentBestThresholdSize;
end