function showXPortParameters(freq, Spar, varargin)
%% showXPortParameters plots data in par to figure hFig with text in legend.
%
% INPUTS:
%   freq: vector of frequency points in Hz, cell
%   par: 3D matrix of parameters obtained from SXPParse, cell
%   legendText: legend text. Cell with title to all parameters.
%               Angle and abs value signs will be added automatically.
%               Part @@ will be removed by number of parameter (11, 21, ...)
%   dB: if absolute value data are in dB or linear, logical
%   LineStyle: like plot use, cell
%   LineWidth: double
%   YAxis: limits of Y axis. Struct with field trans and refl. This parameter is
%          not mandatory.
%
%  SYNTAX:
%
%   legendText.abs = {'|S_{@@}|', '|S2_{@@}|'};
%   legendText.angle = {'\angle(S_{@@})', '\angle(S2_{@@})'};
%   Yaxis.trans = [0 -50];
%   Yaxis.refl = [0 -20];
%   showTwoPortParameters({freq; freq}, {S_SMA; S_SMA2}, true, ...
%      legendText, {'-'; '--'}, [2, 1], Yaxis)
%
% © 2015, Viktor Adler, CTU in Prague, adlervik@fel.cvut.cz


% number of S-parameters
nPar = length(Spar);
% if just single frequency vector is provided, replicate it to all S-parameters
if length(freq) == 1
   freq = repmat(freq, 1, nPar);
end

% nFreq = length(freq{iPar});
if ismatrix(Spar{1}) % just vector [n x 1] with S11 probably, all S-parameters
   % in cell have to be vectors
   nPorts = 1;
   for iPar = 1:nPar % reshape a vector to 1 x 1 x nFreq
      %       nFreq = length(freq{iPar});
      Spar{iPar} = reshape(Spar{iPar}, 1, 1, []);
   end
else
   nPorts = size(Spar{1}, 1);
end

YLabelDefault = repmat({'S_{@@}'}, 1, nPar);

p = inputParser();
p.addParameter('dB', true, @(x)validateattributes(x, {'logical'}, {'scalar'}));

legendTextDefault = repmat({'S_{@@}'}, 1, nPar);
p.addParameter('LegendText', legendTextDefault, @(x)validateattributes(x, {'cell'}, {'nonempty'}));

defaultLineStyles = {'-', '--', '-.'};
nLinestyles = length(defaultLineStyles);
p.addParameter('LineStyle', defaultLineStyles, @(x)validateattributes(x, {'cell'}, {'nonempty'}));

p.addParameter('LineWidth', ones(1, nPar), @(x)validateattributes(x, {'double'}, {'nonempty'}));
% YAxis is struct with refl and trans fields
p.addParameter('YAxis', [], @(x)validateattributes(x, {'struct'}, {'nonempty'}));

p.addParameter('Title', '', @(x)validateattributes(x, {'char'}, {'row'}));

% defaultColorMap = @jet;
% p.addParameter('ColorMap', defaultColorMap, @(x)validateattributes(x, {'function_handle'}, {'scalar'}));

if nPar <= 6
   allColors = {[0 0 1], [1 0 0], [0 1 0], [1 0 1], [0 1 1], [0 0 0]};
   colorDefault = allColors(1:nPar);
else
   colorDefault = jet(nPar);
   colorDefault = mat2cell(colorDefault, ones(nPar, 1), 3);
end

nColors = length(colorDefault);
p.addParameter('Color', colorDefault, @(x)validateattributes(x, {'cell'}, {'nonempty'}));

p.addParameter('XLim', ...
   [min(cellfun(@min, freq)), max(cellfun(@max, freq))], ...
   @(x)validateattributes(x, {'double'}, {'row'}));

p.addParameter('Phase', true, @(x)validateattributes(x, {'logical'}, {'scalar'}));

p.addParameter('Mag', true, @(x)validateattributes(x, {'logical'}, {'scalar'}));

p.addParameter('Symmetric', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));

% p.addParameter('AxisFontSize', 10, @(x)validateattributes(x, {'double'}, {'scalar', 'positive'}));

p.addParameter('LegendFontSize', 9, @(x)validateattributes(x, {'double'}, {'scalar', 'positive'}));

p.addParameter('ShowLegend', true, @(x)validateattributes(x, {'logical'}, {'scalar'}));

p.addParameter('XLabelInterpreter', 'tex', @(x)validateattributes(x, {'char'}, {'row'}));

p.addParameter('YLabelInterpreter', 'tex', @(x)validateattributes(x, {'char'}, {'row'}));

p.addParameter('XLabel', 'f (GHz)', @(x)validateattributes(x, {'char'}, {'row'}));

p.addParameter('AxesLineWidth', 1, @(x)validateattributes(x, {'double'}, {'scalar'}));

p.addParameter('AxesFontSize', 12, @(x)validateattributes(x, {'double'}, {'scalar', 'positive'}));

p.addParameter('GridLineStyle', '-', @(x)validateattributes(x, {'char'}, {'scalar'}));

p.addParameter('FontName', 'Times New Roman', @(x)validateattributes(x, {'char'}, {'scalar'}));

p.addParameter('XTick', [], @(x)validateattributes(x, {'double'}, {'row'}));

p.addParameter('Polar', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));

p.parse(varargin{:});
dB = p.Results.dB;
legendText = p.Results.LegendText;
lineStyle = p.Results.LineStyle;
lineWidth = p.Results.LineWidth;
titleText = p.Results.Title;
showPhase = p.Results.Phase;
showMag = p.Results.Mag;
isSymmetric = p.Results.Symmetric;
legendFontSize = p.Results.LegendFontSize;
lineColor = p.Results.Color;
showLegend = p.Results.ShowLegend;
showPolarGraph = p.Results.Polar;

AxesLineWidth = p.Results.AxesLineWidth;
AxesFontSize = p.Results.AxesFontSize;
GridLineStyle = p.Results.GridLineStyle;
FontName = p.Results.FontName;

YAxis = p.Results.YAxis;
YLabelInterpreter = p.Results.YLabelInterpreter;

XLabelInterpreter = p.Results.XLabelInterpreter;
XLim = p.Results.XLim;
XLabel = p.Results.XLabel;
XTick = p.Results.XTick;

%%

if ~showPolarGraph && nPorts == 1 && showPhase && showMag
   orderOfAxes = [1, 2]; % order of axes in subplot command
   %    x = [1;1]; % the first number of Sxy parameter
   x = repmat(sort(repmat(1:nPorts, 1, nPorts)), 1, 2);
   %    y = [1;1]; % the second number
   y = repmat(repmat(1:nPorts, 1, nPorts), 1, 2);
   axisOrderMag = 1;
   axisOrderPhase = 2;
   nRowAx = nPorts; % 1 rows of axes
   nCollAx = 2*nPorts; % 2 colls of axes
elseif ~showPolarGraph && showPhase && showMag && ~isSymmetric
   nRowAx = nPorts; % 2 rows of axes
   nCollAx = 2*nPorts; % 4 colls of axes
   
   combinationsOfAxes = reshape(1:nRowAx*nCollAx, nCollAx, nRowAx).';
   half1OfAxes = combinationsOfAxes(:, 1:nRowAx);
   half2OfAxes = combinationsOfAxes(:, nRowAx+1:nCollAx);
   
   orderOfAxes = [reshape(half1OfAxes.', 1, []), reshape(half2OfAxes.', 1, [])];
   
   %    orderOfAxes = [1, 2, 5, 6, 3, 4, 7, 8]; % order of axes in subplot command
   %    x = [1;1;2;2;1;1;2;2]; % the first number of Sxy parameter
   x = repmat(sort(repmat(1:nPorts, 1, nPorts)), 1, 2);
   y = repmat(repmat(1:nPorts, 1, nPorts), 1, 2);
   %    y = [1;2;1;2;1;2;1;2]; % the second number
   %    axisOrderMag = 1:4;
   axisOrderMag = 1:nPorts^2;
   %    axisOrderPhase = 5:8;
   axisOrderPhase = max(axisOrderMag) + axisOrderMag;
   
elseif ~showPolarGraph && showPhase && showMag && isSymmetric
   orderOfAxes = [1, 2, 3, 4]; % order of axes in subplot command
   x = [1;2;1;2]; % the first number of Sxy parameter
   y = [1;1;1;1]; % the second number
   axisOrderMag = 1:2;
   axisOrderPhase = 3:4;
   nRowAx = 2; % 2 rows of axes
   nCollAx = 2; % 4 colls of axes
elseif ~showPolarGraph && ~showPhase && showMag && isSymmetric % show magnitude but no phases
   orderOfAxes = [1, 2]; % order of axes in subplot command
   x = [1;2]; % the first number of Sxy parameter
   y = [1;1]; % the second number
   axisOrderMag = 1:2;
   axisOrderPhase = [];
   nRowAx = 1; % 1 rows of axes
   nCollAx = 2; % 2 colls of axes
elseif ~showPolarGraph && ~showPhase && showMag && ~isSymmetric% show magnitude but no phases
   nRowAx = nPorts; % 2 rows of axes
   nCollAx = nPorts; % 2 colls of axes
   orderOfAxes = 1:nPorts^2; % order of axes in subplot command
   x = repelem(1:nPorts, nPorts).';
   %    x = [1;1;2;2]; % the first number of Sxy parameter
   y = repmat((1:nPorts).', nPorts, 1);
   %    y = [1;2;1;2]; % the second number
   axisOrderMag = orderOfAxes;
   axisOrderPhase = [];
   
elseif ~showPolarGraph && showPhase && ~showMag && isSymmetric% show phases but no magnitudes
   orderOfAxes = [1, 2]; % order of axes in subplot command
   x = [1;2]; % the first number of Sxy parameter
   y = [1;1]; % the second number
   axisOrderMag = [];
   axisOrderPhase = 1:2;
   nRowAx = 1; % 1 rows of axes
   nCollAx = 2; % 2 colls of axes
elseif ~showPolarGraph && showPhase && ~showMag && ~isSymmetric% show phases but no magnitudes
   %    orderOfAxes = [1, 2, 3, 4]; % order of axes in subplot command
   orderOfAxes = 1:nPorts^2; % order of axes in subplot command
   %    x = [1;1;2;2]; % the first number of Sxy parameter
   x = repelem(1:nPorts, nPorts).';
   %    y = [1;2;1;2]; % the second number
   y = repmat((1:nPorts).', nPorts, 1);
   axisOrderMag = [];
   axisOrderPhase = orderOfAxes;
   nRowAx = nPorts; % 2 rows of axes
   nCollAx = nPorts; % 2 colls of axes
elseif showPolarGraph
   nRowAx = nPorts; % 2 rows of axes
   nCollAx = nPorts; % 2 colls of axes
   orderOfAxes = 1:nPorts^2; % order of axes in subplot command
   x = repelem(1:nPorts, nPorts).';
   %    x = [1;1;2;2]; % the first number of Sxy parameter
   y = repmat((1:nPorts).', nPorts, 1);
   %    y = [1;2;1;2]; % the second number
   axisOrderMag = orderOfAxes;
   axisOrderPhase = [];
end

legendAllText = cell(nPar, 1);
if nPar > nColors
   lineColor = repmat(lineColor, 1, ceil(nPar/nColors));
end

if nPar > nLinestyles
   lineStyle = repmat(lineStyle, 1, ceil(nPar/nLinestyles));
end

% open new figure, of get the last one
hFig = gcf;
set(hFig, 'Color', [1 1 1]);
% create tiledlayout
hTiledL = tiledlayout(hFig, nRowAx, nCollAx, 'Padding', 'tight', 'TileSpacing', 'tight');

% if showPolarGraph
%    magPlotFcn = @polarplot;
% else
%    magPlotFcn = @plot;
% end

if showPolarGraph
   for iAx = axisOrderMag
      nexttile(hTiledL, orderOfAxes(iAx)); % probably can't return polaraxes object
      for iPar = 1:nPar
         %          polaraxes(hAx);
         polarplot(angle(squeeze(Spar{iPar}(x(iAx), y(iAx), :))), ...
            abs(squeeze(Spar{iPar}(x(iAx), y(iAx), :))), ...
            lineStyle{iPar}, ...
            'LineWidth', lineWidth(iPar), 'Color', lineColor{iPar});
         hAx = gca;
         hold on
         legendAllText{iPar} = replaceDelimiter(legendText{iPar}, ...
            num2str(x(iAx)*10 + y(iAx)));
      end
      showLegendFcn();
      YLabel = replaceDelimiter(YLabelDefault{iPar}, num2str(x(iAx)*10 + y(iAx)));
      hAx.Title.String = YLabel;
   end
   
else % show rectangular graph
   % first 4 axes are for absolute values
   for iAx = axisOrderMag
      % make subplot
      hAx = nexttile(hTiledL, orderOfAxes(iAx));
      %    hAx = subplot(nRowAx, nCollAx, orderOfAxes(iAx));
      for iPar = 1:nPar
         data = abs(squeeze(Spar{iPar}(x(iAx), y(iAx), :)));
         if dB
            data = 20*log10(data);
         end
         plot(hAx, freq{iPar}/1e9, data, ...
            lineStyle{iPar}, ...
            'LineWidth', lineWidth(iPar), 'Color', lineColor{iPar})
         hold on
         grid on
         
         legendAllText{iPar} = replaceDelimiter(['|', legendText{iPar}, '|'], ...
            num2str(x(iAx)*10 + y(iAx)));
         YLabel = replaceDelimiter(YLabelDefault{iPar}, num2str(x(iAx)*10 + y(iAx)));
         if dB
            hAx.YLabel.String = [YLabel, ' (dB)'];
         else
            hAx.YLabel.String = [YLabel, ' (-)'];
         end
      end
      
      showLegendFcn();
      %    set(hAx, 'XLim', XLim/1e9, 'FontSize', axisFontSize);
      if ~isempty(YAxis)
         if x(iAx) == y(iAx) % reflection
            %       YaxisLim = sort([Yaxis.refl; Yaxis.trans; Yaxis.trans; Yaxis.refl], 2);
            set(hAx, 'YLim', sort(YAxis.refl))
         else
            set(hAx, 'YLim', sort(YAxis.trans))
         end
      end
      hAx.XLim = XLim/1e9;
      hAx.XLabel.Interpreter = XLabelInterpreter;
      hAx.XLabel.String = XLabel;
      hAx.LineWidth = AxesLineWidth;
      hAx.FontSize = AxesFontSize;
      hAx.GridLineStyle = GridLineStyle;
      hAx.FontName = FontName;
      hAx.YLabel.Interpreter = YLabelInterpreter;
      if ~isempty(XTick)
         hAx.XTick = XTick/1e9;
      end
   end
end

% next 4 graphs are for angle
for iAx = axisOrderPhase
   hAx = nexttile(hTiledL, orderOfAxes(iAx));
   %    hAx = subplot(nRowAx, nCollAx, orderOfAxes(iAx));
   for iPar = 1:nPar
      
      plot(hAx, freq{iPar}/1e9, angle(squeeze(Spar{iPar}(x(iAx), y(iAx), :)))/pi*180, ...
         lineStyle{iPar}, ...
         'LineWidth', lineWidth(iPar), 'Color', lineColor{iPar})
      hold on
      grid on
      
      legendAllText{iPar} = replaceDelimiter(['\angle(', legendText{iPar}, ')'], ...
         num2str(x(iAx)*10+y(iAx)));
   end
   %    set(hAx, 'XLim', XLim/1e9, 'FontSize', axisFontSize);
   YLabel = replaceDelimiter(YLabelDefault{iPar}, num2str(x(iAx)*10 + y(iAx)));
   hAx.YLabel.String = [YLabel, ' (deg.)'];
%    if showLegend
%       legend(hAx, legendAllText, 'Location', 'best', 'FontSize', legendFontSize)
%    end
   showLegendFcn();
   hAx.XLim = XLim/1e9;
   hAx.YLim = [-200 200];
   hAx.XLabel.Interpreter = XLabelInterpreter;
   hAx.XLabel.String = XLabel;
   hAx.LineWidth = AxesLineWidth;
   hAx.FontSize = AxesFontSize;
   hAx.GridLineStyle = GridLineStyle;
   hAx.FontName = FontName;
   hAx.YLabel.Interpreter = YLabelInterpreter;
end



hTiledL.Title.String = titleText;

   function showLegendFcn()
      if showLegend
         legend(hAx, legendAllText, 'Location', 'best', 'FontSize', legendFontSize)
      end
   end

end

function output = replaceDelimiter(legendText, newText)
% replace '@@' with 12, 22, ...
output = strsplit(legendText, '@@');
output = strjoin(output, newText);
end