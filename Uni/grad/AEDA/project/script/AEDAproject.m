%% Initialization
close all; clear; clc
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";
T = readtable('data.xlsx', 'ReadVariableNames', true);
HC = T(ismember(T.Group, "HC"),:);
MS = T(ismember(T.Group, "MS"),:);
speechParameters = ["DDKR", "DDKI", "stdF0", "jitter", "HNR", "DUS", "RFA", "AR", "IntSD", "F0SD"];
speechParametersValues = table2array(T(:,16:end));
alpha = 0.05;

%% Speech parameters correlation
normality = reportNormality(MS, speechParameters, alpha);
figure('Name', 'Table of speech parameter correlations', 'Units', 'normalized', 'Position', [.25 .3 .5 .4]);
reportCorrelation(MS, speechParameters, normality, alpha);

%% Resulting parameters correlation
result_parameters = ["WhiteMatter", "GrayMatter", "CerebellarWhiteMatter", "CerebellarGrayMatter", "WholeBrainTissue"];
normality_results = reportNormality(MS, result_parameters, alpha);
figure('Name', 'Table of result parameter correlations', 'Units', 'normalized', 'Position', [.25 .3 .5 .4]);
reportCorrelation(MS, result_parameters, normality_results, alpha);

%% Ttest2/ranksum for speech parameters
iiter = 1;
for parameter = speechParameters
    if normality.(parameter)
        [h, p, ci, stats] = ttest2(MS.(parameter), HC.(parameter), 'Tail', 'both', 'Alpha', alpha);
    else
        [p, h , stats] = ranksum(MS.(parameter), HC.(parameter), 'Tail', 'both', 'Alpha', alpha);
    end
    hs_same(iiter) = h; % 1 - jine rozdeleni, 0 - stejne rozdeleni
    % parametr IntSD = vychází stejně pro zdravé a nemocné
        % nebudeme ho používat 
    iiter = iiter + 1;
end

%% Vizualizace

figure('Name', 'KDE')
t1 = tiledlayout(2, 5, "TileSpacing","tight");
for names = speechParameters
    [f, xf] = kde(MS.(names));
    nexttile
    area(xf, f, 'LineWidth',1.5)
    ylabel('KDE (-)')
    title([names, ' Normality: ', normality.(names)])
    grid on
end

%% N-way ANOVA
for parameter = result_parameters
    [p, tbl, stats, terms] = anovan(MS.(parameter),{MS.Sex,MS.EDSS, MS.Age}, 'varnames', {'Sex', 'EDSS', 'Age'}, 'model','interaction', 'continuous', 2:3, 'display', 'off');
    figure('Name', 'N-way ANOVA for parameter ' + parameter, 'Units', 'normalized', 'Position', [.25 .4 .5 .3])
    uitable('Data', tbl(2:end,2:end), 'RowName', tbl(2:end,1), 'ColumnName', tbl(1,2:end), 'Units', 'Normalized', 'Position', [.05, .05, .9, .9]);
    
    table2latex(tbl, 'table.tex');
end

%% LASSO
Lasso_Param = zeros(10,5);
for i = 11:15
    [B, FitInfo] = lasso(table2array(MS(:,16:end)), table2array(MS(:,i)), 'CV', 10);
    % lassoPlot(B, FitInfo, 'PlotType', 'CV');
    lambda_optimal = FitInfo.LambdaMinMSE;
    Lasso_Param(:,i-10) = B(:,FitInfo.Lambda==lambda_optimal);
end

%% Linear regression model 
fitlm(cat(2,MS.stdF0, MS.DUS, MS.AR),MS.WhiteMatter, 'RobustOpts','on');
fitlm(MS.AR,MS.GrayMatter, 'RobustOpts','on');
fitlm(cat(2,MS.DDKR,MS.jitter,MS.HNR, MS.AR),MS.CerebellarWhiteMatter, 'RobustOpts','on');
fitlm(cat(2,MS.DDKR, MS.jitter),MS.CerebellarGrayMatter, 'RobustOpts','on');
fitlm(cat(2,MS.stdF0, MS.jitter,MS.AR),MS.WholeBrainTissue, 'RobustOpts','on');

fitlm(cat(2,MS.EDSS, MS.SDMT, MS.MSFC_PASAT3, MS.MSFC_T25FW, MS.MSFC_9HPT, MS.BDI),MS.WhiteMatter, 'RobustOpts','on');
fitlm(cat(2,MS.EDSS, MS.SDMT, MS.MSFC_PASAT3, MS.MSFC_T25FW, MS.MSFC_9HPT, MS.BDI),MS.GrayMatter, 'RobustOpts','on');
fitlm(cat(2,MS.EDSS, MS.SDMT, MS.MSFC_PASAT3, MS.MSFC_T25FW, MS.MSFC_9HPT, MS.BDI),MS.CerebellarWhiteMatter, 'RobustOpts','on');
fitlm(cat(2,MS.EDSS, MS.SDMT, MS.MSFC_PASAT3, MS.MSFC_T25FW, MS.MSFC_9HPT, MS.BDI),MS.CerebellarGrayMatter,'RobustOpts','on');
fitlm(cat(2,MS.EDSS, MS.SDMT, MS.MSFC_PASAT3, MS.MSFC_T25FW, MS.MSFC_9HPT, MS.BDI),MS.WholeBrainTissue, 'RobustOpts','on');

%% Functions
% Report normality
function normality = reportNormality(table, parameters, alpha)
    for sample = parameters
        [H, pValue, SWstatistic] = swtest(table.(sample), alpha);
        normality.(sample) = ~H;
    
        if pValue < 0.001
            pReport = "p < 0.001";
        elseif pValue < 0.01
            pReport = sprintf('p = %.3f', pValue);
        else
            pReport = sprintf('p = %.2f', pValue);
        end
    
        if H == 0
            fprintf(['\tShapiro-Wilk test confirms that the null hypothesis of composite normality is a reasonable assumption\n' ...
                'regarding the population distribution of a random sample %s (M = %.2f, SD = %.2f)\n' ...
                'at significance level α = %.2f, W(%d) = %.4f, %s.\n\n'], ...
                sample, mean(table.(sample)), std(table.(sample)), alpha, length(table.(sample))-1, SWstatistic, pReport)
        else
            fprintf(['\tShapiro-Wilk test rejects the null hypothesis of composite normality regarding the population\n' ...
                'distribution of a random sample %s (M = %.2f, SD = %.2f)\n' ...
                'at significance level α = %.2f, W(%d) = %.4f, %s.\n\n'], ...
                sample, mean(table.(sample)), std(table.(sample)), alpha, length(table.(sample))-1, SWstatistic, pReport)
        end
    end
end

% Report correlation
function reportCorrelation(table, parameters, normality, alpha)
    Nparameters = length(parameters);
    resultsTable = array2table(nan(Nparameters), 'RowNames', parameters, 'VariableNames', parameters);
    corrTypes = ["Spearman", "Pearson"];
    for sample1 = parameters
        for sample2 = parameters(1:find(parameters==sample1))
            type = corrTypes(1 + (normality.(sample1) && normality.(sample2)));
            [rho, pValue] = corr(table.(sample1), table.(sample2), 'Type', type);
    
            if pValue < 0.001
                pReport = "p < 0.001";
            elseif pValue < 0.01
                pReport = sprintf('p = %.3f', pValue);
            else
                pReport = sprintf('p = %.2f', pValue);
            end
    
            if pValue < alpha
                resultsTable{sample1, sample2} = rho;
                fprintf("A %s correlation test was conducted to assess the relationship between two parameters:\n" + ...
                    "\t%s (M=%.3f, SD=%.3f),\n" + ...
                    "\t%s (M=%.3f, SD=%.3f).\n" + ...
                    "Considering the level of significance α = %.2f, there was a significant relationship between\n" + ...
                    "these parameters, r(%d) = %.3f, %s.\n\n", ...
                type, sample1, mean(table.(sample1)), std(table.(sample1)), sample2, mean(table.(sample2)), std(table.(sample2)), ...
                alpha, height(table)-2, rho, pReport);
            else
                resultsTable{sample1, sample2} = 0;
            end
        end
    end

    uit = uitable("Data", resultsTable{:,:}, 'RowName', parameters, 'ColumnName', parameters, 'Units', 'Normalized', 'Position', [.05, .05, .9, .9]);
    for row = 1:Nparameters
        for col = 1:Nparameters
            if isnan(resultsTable{row,col})
                style = uistyle('BackgroundColor', "#A9A9A9");
            elseif abs(resultsTable{row,col}) >= 0.8
                style = uistyle('BackgroundColor', "#FFAC1C");
            elseif abs(resultsTable{row,col}) >= 0.6
                style = uistyle('BackgroundColor', "#FFBF00");
            elseif abs(resultsTable{row,col}) >= 0.4
                style = uistyle('BackgroundColor', "#FFC000");
            elseif abs(resultsTable{row,col}) >= 0.2
                style = uistyle('BackgroundColor', "#FFD580");
            else
                style = uistyle('BackgroundColor', "#FFDEAD");
            end
            addStyle(uit, style, 'cell', [row,col]);
        end
    end
end


