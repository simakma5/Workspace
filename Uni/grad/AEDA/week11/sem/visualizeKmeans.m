function [F1, F2] = visualizeKmeans(IDX,CENTROIDS,DATA)
narginchk(3,3)

% Input checks
if size(DATA,1) ~= size(IDX,1)
    error("Index and data row length doesn't match! %d vs %d",...
        size(IDX,1),size(DATA,1))
end
if size(CENTROIDS,1) ~= 2
    error("Invalid number of centroid rows! 2 needed.")
elseif size(CENTROIDS,2) ~= 4
    error("Invalid number of centroid columns! 4 needed.")
end
if size(DATA,2) ~= 4
    error("Invalid number of data columns! 4 needed.")
end
    
% Scatter 3D visualization
F1 = figure('Name','AED11: Data Visualisation 1', ...
    'Units','Normalized','OuterPosition',[0, 0, 1, 1]);
TL = tiledlayout(2,2,'TileSpacing','compact');
title(TL,{'K-Means: 2-Cluster data visualization',...
    'Δ = Untreated Score - Treated Score'})

perms = nchoosek(1:4,3); % All possible permutations
SETS = {DATA(:,1) DATA(:,2) DATA(:,3) DATA(:,4)};

axL = {'UPDRS Axial Score Δ (-)',...
    'UPDRS Tremor Score Δ (-)',...
    'UPDRS Bradykinesia Score Δ (-)',...
    'UPDRS Rigidity Score Δ (-)'};
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
titL = {'Axial','Tremor','Bradykinesia','Rigidity'};
for n = 1:length(SETS)
    D = [SETS{perms(n,:)}];
    L = axL(perms(n,:));
    nexttile
    grid on
    hold on
    
    H = gobjects(1,2);
    H_CNTR = gobjects(1,2);
    for i = 1:2 % Over both clusters
        gix = (IDX == i);
        % Datapoints scatter
        H(i) = scatter3(D(gix,1),D(gix,2),D(gix,3),...
            80,'filled','MarkerEdgeColor',...
            'black','MarkerFaceColor',colors{i});
        % Centroids scatter
        H_CNTR(i) = scatter3(CENTROIDS(i,perms(n,1)),...
            CENTROIDS(i,perms(n,2)),...
            CENTROIDS(i,perms(n,3)),...
            80,'filled','Marker','*','MarkerEdgeColor',colors{i});
    end
    
    hold off
    xlabel(L{1}),ylabel(L{2}),zlabel(L{3});
    title(sprintf('%s + %s + %s',...
        titL{perms(n,1)},titL{perms(n,2)},titL{perms(n,3)}))
    ax = gca;
    ax.View = [-17.4514 13.4453];
    ax.Projection = 'perspective';
    ax.Box = 'on';
end

lgd = legend([H,H_CNTR],'Cluster 1','Cluster 2',...
    'Centroid 1','Centroid 2',...
    'Orientation','horizontal');
title(lgd,'Clusters')
lgd.Layout.Tile = 'South';

% Plot Matrix
titL = {'Δ Axial','Δ Tremor','Δ Bradykinesia','Δ Rigidity'};
F2 = figure('Name','AED11: Data Visualisation 2', ...
    'Units','Normalized','OuterPosition',[0, 0, 1, 1]);
[~,ax] = gplotmatrix(DATA,[],IDX,...
    [0 0.4470 0.7410;0.8500 0.3250 0.0980],[],[],'off','variable',titL);
title({'K-Means: 2-Cluster data plot matrix',...
    'Δ = Untreated Score - Treated Score'})

% Turn grid on in all subplots
for i = 1:4 
    for j = 1:4
        if i ~= j
            ax(i,j).XGrid = 'on';
            ax(i,j).YGrid = 'on';
        end
    end
end

end