%% svmLim3DboundaryPlot.m - Plot Linear SVM Decision Boundary in 3D
% Use this function to plot input 3D data from two groups as scatter points
% along with a decision boundary (2D plane) of the provided trained support
% vector machine model.
%
% Input arguments:
% - data_norm: Normalized data. Data must be normalized to zscores and must
%              be in a form of matrix NxP, where N is number of data points
%              and P is the number of parameters, in this case, P = 3.
% - group: Vector of group labels, length Nx1. Must contain two distinct
%          labels, "HC" and "PD".
% - SVMmodel3DLin: Trained SVM model with a linear kernel, created by the 
%                  fitcsvm function.
% - pNames: (Optional argument) String array of parameter names.
%
% Output arguments:
% - Farr: Structure array containing the figure and graphic object handles.
%
% Created 8. 12. 2023 by Petr Krýže, FEE CTU Prague, for the purposes of
% the course Experimental Data Analysis.
function [Farr] = svmLin3DboundaryPlot(data_norm, group, SVMmodel3DLin, pNames)
% Input check
narginchk(3,4);

if nargin <= 3
    pNames = ["x1","x2","x3"];
end

assert(size(data_norm,2) == 3,"Data must have 3 columns.")
assert(size(data_norm,1) == length(group),"Group must be the same length as the number of rows in data.")
assert(length(pNames) == 3,"pNames must be a string array of length 3.")

if ~isa(group,"categorical")
    group = categorical(group);
    group = group(:);
end

% Group specific data
D_HC_norm = data_norm(group == "HC",:);
D_PD_norm = data_norm(group == "PD",:);

% Model parameters
b0 = SVMmodel3DLin.Bias;
b1 = SVMmodel3DLin.Beta(1); 
b2 = SVMmodel3DLin.Beta(2);
b3 = SVMmodel3DLin.Beta(3);

% Grid
xrange = linspace(-3,3,25);
yrange = linspace(-3,3,25);
[xxlin,yylin] = meshgrid(xrange,yrange);

% Linear Kernel SVM Decision boundary (plane) equation
% 0 = b0 + b1*x + b2*y + b3*z
% z = (-b0-b1*x-b2*y)/b3
z_lin = (-b0-b1*xxlin-b2*yylin)/b3;

[clasif_labels_lin,~] = predict(SVMmodel3DLin,data_norm);

% Plotting
F1 = figure("Name","AED10: 3D SVM w/ Linear Kernel Decision Boundary");
TL = tiledlayout(1,2);
title(TL, "3D SVM w/ Linear Kernel Decision Boundary")

% 1) Linear plot with real labels
nexttile
H1 = scatter3(D_HC_norm(:,1),D_HC_norm(:,2),D_HC_norm(:,3),100,...
    "MarkerFaceColor",[0 0 1],"MarkerEdgeColor",[0 0 0],"LineWidth",1.5);
hold on
H2 = scatter3(D_PD_norm(:,1),D_PD_norm(:,2),D_PD_norm(:,3),100,...
    "MarkerFaceColor",[1 0 0],"MarkerEdgeColor",[0 0 0],"LineWidth",1.5);

DB1 = surf(xxlin,yylin,z_lin,"FaceColor",[0 0 0],"FaceAlpha",0.3,"EdgeAlpha",0.5);

title("True Labels")
xlabel(sprintf("Normalized %s",pNames{1}));
ylabel(sprintf("Normalized %s",pNames{2}));
zlabel(sprintf("Normalized %s",pNames{3}));
grid on
ax = gca;
ax.Box = "on";
zlim([-3, 7])
legend([H1,H2,DB1],...
    "True HC",...
    "True PD",...
    "Decision boundary (Plane)",...
    "Location","northeast")

% 2) Linear plot with classification labels
nexttile
H3 = scatter3(data_norm(clasif_labels_lin=="HC",1),...
    data_norm(clasif_labels_lin=="HC",2),...
    data_norm(clasif_labels_lin=="HC",3),...
    100,...
    "MarkerFaceColor",[0 0 1],"MarkerEdgeColor",[0 0 0],"LineWidth",1.5);
hold on
H4 = scatter3(data_norm(clasif_labels_lin=="PD",1),...
    data_norm(clasif_labels_lin=="PD",2),...
    data_norm(clasif_labels_lin=="PD",3),...
    100,...
    "MarkerFaceColor",[1 0 0],"MarkerEdgeColor",[0 0 0],"LineWidth",1.5);

DB2 = surf(xxlin,yylin,z_lin,"FaceColor",[0 0 0],"FaceAlpha",0.3,"EdgeAlpha",0.5);

SVS = SVMmodel3DLin.SupportVectors; 
SVH = plot3(SVS(:,1),SVS(:,2),SVS(:,3),"+k","MarkerSize",10);

title("Classification Labels")
xlabel(sprintf("Normalized %s",pNames{1}));
ylabel(sprintf("Normalized %s",pNames{2}));
zlabel(sprintf("Normalized %s",pNames{3}));
grid on
ax = gca;
ax.Box = "on";
zlim([-3, 7])
legend([H3,H4,DB2,SVH],...
    "Classified as HC",...
    "Classified as PD",...
    "Decision boundary (Plane)",...
    "Support Vectors",...
    "Location","northeast")

% Confusion matrix
F2 = figure("Name","AED10: 3D SVM w/ Linear Kernel Classification results");
C_lin = confusionmat(group,clasif_labels_lin);
confusionchart(C_lin,["HC","PD"])
title("3D SVM w/ Linear Kernel Classification results")

% Return Handles
GRAPH1.FigureHandle = F1;
GRAPH1.PlotHandles = [H1,H2,DB1];
GRAPH2.FigureHandle = F2;
GRAPH2.PlotHandles = [H3,H4,DB2,SVH];
Farr = [GRAPH1,GRAPH2];

end