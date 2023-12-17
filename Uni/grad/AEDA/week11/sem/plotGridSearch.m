function peaks = plotGridSearch(X,Y,Z)
% AED - Experimental Data Analysis, Supplementary course script
% -------------------------------------------------------------------------
% Function that plots optimization grids for SVM binary classifier, along
% with the optimal hyperparameter peak points.

% INPUT:
% X is a 1xN vector of all possible values for parameter 1 (e.g. BoxConstraint)
% Y is a 1xN vector of all possible values for parameter 2 (e.g. KernelScale)
% Z(i,j,k) is a NxNx3 matrix, where:
%   > i = 1,2,...,N
%   > j = 1,2,...,N
%   > k = 1,2,3; where 1==Sensitivity, 2==Specificity, 3==Accuracy
%   > and Z(i,j,1/2/3) is the sensitivity/specificity/accuracy for the i-th
%     value of BoxConstraint (C) in X and j-th value of KernelScale (σ) in Y

% OUTPUT:
% Figure with plotted optimization surfaces and peak points.
% peaks = table with the optimal values.
% -------------------------------------------------------------------------
peak_points = nan(3); % Output table init

figure('Name','Grid Search','Units','Normalized','OuterPosition',[0, 0, 1, 1]); % Figure and layout init
TL = tiledlayout(1,3,'TileSpacing','compact');
title(TL,'Grid Search - Hyperparameter optimization on SVM with RBF kernel')
out_names = {'Sensitivity','Specificity','Accuracy'};
for i = 1:3
    nexttile
    % Plotting the optimization surface    
    surf(X,Y,Z(:,:,i)','FaceLighting','gouraud'); hold on;
    
    % Locating and plotting the surface peak
    [max_from_col, idx_max_BoxConstraints_rows] = max(Z(:,:,i)); % Column maxima
    [max_peak, idx_max_KernelScale] = max(max_from_col); % Maxima from column maxima
    
    idx_max_BoxConstraint = idx_max_BoxConstraints_rows(idx_max_KernelScale);
    maxBoxConstraint = X(idx_max_BoxConstraint);
    maxKernelScale = Y(idx_max_KernelScale);
    
    % Save the peak point value and location
    peak_points(i,:) = [maxBoxConstraint, maxKernelScale, max_peak];
    % Plot the peak point using peak value and location
    H_peak = scatter3(maxBoxConstraint, maxKernelScale, max_peak,...
        120,'o','Filled','MarkerEdgeColor','black','MarkerFaceColor','red');
    
    % Axis labeling and configuration
    xlabel('BoxConstraint (-)');
    ylabel('KernelScale (-)');
    zlabel(sprintf('%s (-)',out_names{i}));
    title({sprintf('%s optimization surface',out_names{i}),...
        sprintf('Optimal point: C = %.2f, σ = %.2f, z = %.2f',...
        maxBoxConstraint, maxKernelScale, max_peak)})
    ax = gca;
    ax.Projection = 'perspective';
    ax.Box = 'on';
    uistack(H_peak,'top');
end
lgd = legend(H_peak,'Optimal parameter peak','Orientation','horizontal');
lgd.Layout.Tile = 'south'; % Places legend to shared layout space

% Formation of output table
peaks = array2table(peak_points,...
    'RowNames',out_names,...
    'VariableNames',{'BoxConstraint','KernelScale','Output'});
end