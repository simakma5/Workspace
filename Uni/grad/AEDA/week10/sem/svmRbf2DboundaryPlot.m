%% svmRbf2DboundaryPlot.m - Plot RBF SVM Decision Boundary in 2D
% Use this function to plot input 2D data from two groups as scatter points
% on a 3D functional surface in RBF kernel space, along with a decision 
% boundary (non-linear curve) of the provided trained support vector 
% machine model (intersection of the functional surface and a z=0 plane).
%
% Input arguments:
% - data_norm: Normalized data. Data must be normalized to zscores and must
%              be in a form of matrix NxP, where N is number of data points
%              and P is the number of parameters, in this case, P = 2.
% - group: Vector of group labels, length Nx1. Must contain two distinct
%          labels, "HC" and "PD".
% - SVMmodel3DLin: Trained SVM model with a RBF kernel, created by the 
%                  fitcsvm function.
% - pNames: (Optional argument) String array of parameter names.
%
% Output arguments:
% - Farr: Structure array containing the figure and graphic object handles.
%
% Created 8. 12. 2023 by Petr Krýže, FEE CTU Prague, for the purposes of
% the course Experimental Data Analysis.
function Farr = svmRbf2DboundaryPlot(data_norm, group, SVMmodel2DRBF, pNames)
% Input check
narginchk(3,4);

if nargin <= 3
    pNames = ["x1","x2"];
end

assert(size(data_norm,2) == 2,"Data must have 2 columns.")
assert(size(data_norm,1) == length(group),"Group must be the same length as the number of rows in data.")
assert(length(pNames) == 2,"pNames must be a string array of length 2.")

if ~isa(group,"categorical")
    group = categorical(group);
    group = group(:);
end

% Group specific data
D_HC_norm = data_norm(group == "HC",:);
D_PD_norm = data_norm(group == "PD",:);
NHC = sum(group == "HC");
NPD = sum(group == "PD");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D SVM RBF Kernel

% Prediction
[clasif_labels_rbf2D,~] = predict(SVMmodel2DRBF,data_norm);

% Grid
xrange = linspace(-5,5,100);
yrange = linspace(-5,5,100);
[xxrbf,yyrbf] = meshgrid(xrange,yrange);

% Model parameters
gamma = SVMmodel2DRBF.KernelParameters.Scale;
suppVs = SVMmodel2DRBF.SupportVectors;
NsuppV = length(suppVs);
alphas = SVMmodel2DRBF.Alpha;
SVMlabels = SVMmodel2DRBF.SupportVectorLabels;
bias = SVMmodel2DRBF.Bias;

% RBF functional surface calculation
% z-values initialization
z_rbf = ones(length(xrange),length(yrange))*bias;
z_points_HC = ones(NHC,1)*bias;
z_points_PD = ones(NPD,1)*bias;

% Loop over all support vectors
for i = 1:NsuppV
    sV = suppVs(i,:);
    for j = 1:length(xrange)
        for k = 1:length(yrange)
            xp = [xrange(j),yrange(k)];
            z_rbf(k,j) = z_rbf(k,j) + alphas(i)*SVMlabels(i)*RBFkernel(gamma, xp, sV);
        end
    end

    % Calculation of the z-coordinate for all data points
    for j = 1:NHC
        xp = D_HC_norm(j,1:2);
        z_points_HC(j) = z_points_HC(j) + alphas(i)*SVMlabels(i)*RBFkernel(gamma, xp, sV);
    end

    for j = 1:NPD
        xp = D_PD_norm(j,1:2);
        z_points_PD(j) = z_points_PD(j) + alphas(i)*SVMlabels(i)*RBFkernel(gamma, xp, sV);
    end
end

% Calculation of the prediction score for the whole grid (for contour)
xyGrid = [xxrbf(:),yyrbf(:)];
[~,scoresRBF] = predict(SVMmodel2DRBF,xyGrid);

% Plotting
F1 = figure('Name',"AED10: 2D SVM w/ RBF Kernel Decision Boundary");
TL = tiledlayout(1,2);
title(TL, "2D SVM w/ RBF Kernel Decision Boundary")

% 1) Plot with true labels
nexttile
H1 = scatter3(D_HC_norm(:,1),D_HC_norm(:,2),z_points_HC,100,...
    'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0],"LineWidth",1.5);
hold on
H2 = scatter3(D_PD_norm(:,1),D_PD_norm(:,2),z_points_PD,100,...
    'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0],"LineWidth",1.5);

RBFH1 = surf(xxrbf,yyrbf,z_rbf,"EdgeAlpha",0.2,"FaceAlpha",0.35,"EdgeColor",[0 0 0]);

[~,DB1] = contour(xxrbf,yyrbf,reshape(scoresRBF(:,2),size(xxrbf)),[0 0],...
    "LineWidth",2,"Color","black","LineStyle","-","ZLocation",0);

title("True Labels")
xlabel(sprintf("Normalized %s",pNames(1)));
ylabel(sprintf("Normalized %s",pNames(2)));
zlabel("RBF Kernel Space");
grid on
ax = gca;
ax.Box = 'on';
legend([H1,H2,RBFH1,DB1],...
    "True HC",...
    "True PD",...
    "RBF Surface",...
    "Decision Boundary",...
    'Location','northeast');

% 2) Plot with classification labels
nexttile
z_points_all = [z_points_HC;z_points_PD];

H1 = scatter3(data_norm(clasif_labels_rbf2D == "HC",1),...
    data_norm(clasif_labels_rbf2D == "HC",2),...
    z_points_all(clasif_labels_rbf2D == "HC"),100,...
    'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0],"LineWidth",1.5);
hold on
H2 = scatter3(data_norm(clasif_labels_rbf2D == "PD",1),...
    data_norm(clasif_labels_rbf2D == "PD",2),...
    z_points_all(clasif_labels_rbf2D == "PD"),100,...
    'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0],"LineWidth",1.5);

RBFH2 = surf(xxrbf,yyrbf,z_rbf,"EdgeAlpha",0.1,"FaceAlpha",0.1,"EdgeColor",[0 0 0]);

[~,DB2] = contour(xxrbf,yyrbf,reshape(scoresRBF(:,2),size(xxrbf)),[0 0],...
    "LineWidth",2,"Color","black","LineStyle","-","ZLocation",0);
contour3(xxrbf,yyrbf,reshape(scoresRBF(:,2),size(xxrbf)),30,...
    "LineWidth",1.5,"Color","flat","ZLocation",0);

SVS = SVMmodel2DRBF.SupportVectors; 
SVH = plot3(SVS(:,1),SVS(:,2),z_points_all(SVMmodel2DRBF.IsSupportVector),...
    "+g","MarkerSize",10);

title('Classification Labels')
xlabel(sprintf("Normalized %s",pNames(1)));
ylabel(sprintf("Normalized %s",pNames(2)));
zlabel("RBF Kernel Space");
grid on
ax = gca;
ax.Box = 'on';
legend([H1,H2,RBFH2,DB2,SVH],...
    "Classified as HC",...
    "Classified as PD",...
    "RBF Surface",...
    "Decision Boundary",...
    "Support Vectors",...
    'Location','northeast');

% Confusion matrix
F2 = figure("Name","AED10: 2D SVM w/ RBF Kernel Classification results");
C_rbf2D = confusionmat(group,clasif_labels_rbf2D);
confusionchart(C_rbf2D,{'HC','PD'})
title("2D SVM w/ RBF Kernel Classification results")

% Return Handles
GRAPH1.FigureHandle = F1;
GRAPH1.PlotHandles = [H1,H2,RBFH1,DB1];
GRAPH2.FigureHandle = F2;
GRAPH2.PlotHandles = [H1,H2,RBFH2,DB2,SVH];
Farr = [GRAPH1,GRAPH2];

%% Functions
function ker = RBFkernel(gamma, x, x_i)
    ker = exp(-(norm(x/gamma - x_i/gamma)^2));
end

end