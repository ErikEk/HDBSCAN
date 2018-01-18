clear all; close all; clc;

%% Parameters
modelRangeInterval = '10m'; % Select either '5m' or '10m'
numSC = 2; % Only for K-means simple model!
minRCS = -30; %[dBSm]
maxRCS = 30; %[dBdraw_circleSm]
plotMarkerSize = 5; % Size of the plot points.
targetLength = 4.369; % meters 4.269
targetWidth = 2.041; % meters (w/ mirrors)
im=imread('images/v40.png');

%% Variables
modelOutput = [];
angleDict = [0; 30; 60; 90; 120; 150; 180; 210; 240; 270; 300; 330];
cmap = colormap(jet); % RCS colormap for colorbar.
close; % Close the figure window that opens when using cmap.

%% Load radar detection data
load(['data/v40data_',  modelRangeInterval]);

%% String dicts for plot titles and user 
if strcmp(modelRangeInterval, '10m')
    rangeSegments = {'<20 m', '20m to 30m', '30m to 40m', '>40m'};
    rangeStr = {'_range_10m-20m_', '_range_20m-30m_', '_range_30m-40m_', '_range_40m-50m_','0'};
elseif strcmp(modelRangeInterval, '5m')
    rangeSegments = {'<10 m', '10m to 15m', '15m to 20m', '20m to 25m', '25m to 30m', '30m to 35m', '35m to 40m', '40m to 45m', '45m to 50m', '>50m'};
    rangeStr = {'_range_5-10m', '_range_10m-15m_', '_range_15m-20m_', '_range_20-25m_', '_range_25m-30m_', '_range_30-35m_', '_range_35m-40m_', '_range_40-45m_', '_range_45m-50m_', '_range_50-55m_'};
end

%angleSegments = {'345deg to 15deg','15deg to 45deg', '45deg to 75deg', '75deg to 105deg', '105deg to 135deg', '135deg to 165deg', '165deg to 195deg', '195deg to 225deg', '225deg to 255deg', '255deg to 285deg', '285deg to 315deg', '315deg to 345deg'};

%% Figure output folder
figureFolder = strcat('figures/', mfilename);

if ~exist(figureFolder,'dir')
    % Folder does not exist yet, so create it.
    mkdir(figureFolder);
end

figHandle = figure('Position', [500, 70, 1000, 800]); % 

%% Variables for radar position plot

for rangeInd = 1:length(rangeSegments)
    
    for angleInd = 1:length(angleDict)
        
        angleDeg = angleDict(angleInd);
        angle = angleDeg*(pi/180);
        
        % Load the model struct into variables.
        x = V40backprojections.x{rangeInd,angleInd};
        y = V40backprojections.y{rangeInd,angleInd};
        RCS = V40backprojections.RCS{rangeInd,angleInd};
        SNR = V40backprojections.SNR{rangeInd,angleInd};
        RCSCol = V40backprojections.RCSCol{rangeInd,angleInd};
        
        %% Plot the current range and angle segment being clustered for clarity (can be skipped)
        
        clf;
        drawModelBins(angleInd, rangeInd, modelRangeInterval)
        pause(0.5)
        clf;
        hold on;
        
        %% Plot detections in bin with RCS colormap
        
        % Plot direction arrow
        arrow1 = [0 -6];
        arrow2 = [0 -5.5];
        
        arrow1x = (arrow1(1)).*cos(angle) - (arrow1(2)).*sin(angle);
        arrow1y = (arrow1(1)).*sin(angle) + (arrow1(2)).*cos(angle);
        arrow2x = (arrow2(1)).*cos(angle) - (arrow2(2)).*sin(angle);
        arrow2y = (arrow2(1)).*sin(angle) + (arrow2(2)).*cos(angle);
        
        arrow1 = [arrow1x arrow1y];
        arrow2 = [arrow2x arrow2y];
        
        darrow = arrow2-arrow1;
        quiver(arrow1(1),arrow1(2),darrow(1),darrow(2),6, 'color', 'k', 'linewidth', 2) % Draw the direction arrow
        
        % Plot image of the target
        image([-targetWidth, targetWidth]/2,[targetLength, -targetLength]/2,im);
        
        % Plot the detections from the current bin with RCS colormap
        scatter(y, x, plotMarkerSize*3, RCSCol, 'filled');
        
        title(['Detections for ', num2str(angleDeg), ' deg and ', rangeSegments{rangeInd}])
        xlabel('Lateral position [m]')
        ylabel('Relative range [m]')
        set(gca,'YDir','normal');
        axis([-4 4 -4 4])
        
        c1 = colorbar;
        colormap(jet)
        caxis([minRCS maxRCS])
        ylabel(c1,'RCS (dBsm)')
        
        % Save as png
        plotType = 'backproj_all';
        figPrintName = strcat(figureFolder, '/', rangeStr{rangeInd}, '_', num2str(angleDeg), 'deg', plotType);
        print(figHandle,'-dpng',figPrintName);
        clf;
        hold on;
        
        %% Clustering of detections to create scattering centers.
       
        image([-targetWidth, targetWidth]/2,[targetLength, -targetLength]/2,im); % Draw the target location
        quiver(arrow1(1),arrow1(2),darrow(1),darrow(2),6, 'color', 'k', 'linewidth', 2) % Draw the direction arrow        
        fprintf('\n');
        disp(['Clustering detections for angle: ', num2str(angleDeg), ' deg and range: ', rangeSegments{rangeInd}]);

        [scatteringCenters,numSC] = HDBSCANClustering(x, y, SNR,3,4,0.9); %%% <-------------- HDBSCAN clustering function call
        
        numSCmatrix(rangeInd,angleInd) = length(scatteringCenters.scSNR); % Save model parameters to matrices.
        for scInd = 1:numSC

            scX(rangeInd,angleInd,scInd) = scatteringCenters.scX(scInd);
            scY(rangeInd,angleInd,scInd) = -scatteringCenters.scY(scInd); % ATT! Inverted y, this is done to match the coordinate system used in the simulations
            scXstddev(rangeInd,angleInd,scInd) = scatteringCenters.scXstddev(scInd);
            scYstddev(rangeInd,angleInd,scInd) = scatteringCenters.scYstddev(scInd);
            scSNR(rangeInd,angleInd,scInd) = scatteringCenters.scSNR(scInd);
            scSNRstddev(rangeInd,angleInd,scInd) = scatteringCenters.scSNRstddev(scInd);
        end
        
        % Save as png
        plotType = 'clusteringResult';
        figPrintName = strcat(figureFolder, '/', rangeStr{rangeInd}, '_', num2str(angleDeg), 'deg', plotType);
        print(figHandle,'-dpng',figPrintName);
                       
    end % end of for angleInd
end % of for rangeInd

for rangeIter = 1:length(rangeSegments) % Put the number of scattering centers first in the model vector
    for angleIter = 1:length(angleDict)
        modelOutput = [modelOutput numSCmatrix(rangeIter,angleIter)];
    end
end

for rangeIter = 1:length(rangeSegments) % Put the model values after the number of corresponding scattering centeres, in the order: X, Y, Xstddev, Ystddev, SNR, SNRstddev
    for angleIter = 1:length(angleDict)
        for scIter = 1:numSCmatrix(rangeIter,angleIter)
            modelOutput = [modelOutput scX(rangeIter,angleIter,scIter) scY(rangeIter,angleIter,scIter) scXstddev(rangeIter,angleIter,scIter) scYstddev(rangeIter,angleIter,scIter) scSNR(rangeIter,angleIter,scIter) scSNRstddev(rangeIter,angleIter,scIter)];
        end
    end
end

dlmwrite('scatteringCenters.csv',modelOutput, ',')
% save('modelValues.mat', 'modelOutput'); % Optional .mat save for debugging

disp('Done!')