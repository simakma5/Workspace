function [freq, Port1Meas, Port2Meas, SThruMeas, SLineMeas] = loadData(folderName)

[freq, SOpen] = SXPParse(fullfile(folderName, 'Open.s2p'));
[~, SShort] = SXPParse(fullfile(folderName, 'Short_A1.s2p'));
[~, SMatch] = SXPParse(fullfile(folderName, 'Match_A1.s2p'));
[~, SThruMeas] = SXPParse(fullfile(folderName, 'Thru_A1.s2p'));
[~, SLine1Meas] = SXPParse(fullfile(folderName, 'Line_A1.s2p'));
[~, SLine2Meas] = SXPParse(fullfile(folderName, 'Line_A2.s2p'));
[~, SLine3Meas] = SXPParse(fullfile(folderName, 'Line_A3.s2p'));
[~, SLine4Meas] = SXPParse(fullfile(folderName, 'Line_A4.s2p'));
[~, SLine5Meas] = SXPParse(fullfile(folderName, 'Line_A5.s2p'));

SLineMeas = cat(4, SLine1Meas, SLine2Meas, SLine3Meas, SLine4Meas, SLine5Meas);

Port1Meas = struct('open', squeeze(SOpen(1, 1, :)), ...
    'short', squeeze(SShort(1, 1, :)), ...
    'match', squeeze(SMatch(1, 1, :)));

Port2Meas = struct('open', squeeze(SOpen(2, 2, :)), ...
    'short', squeeze(SShort(2, 2, :)), ...
    'match', squeeze(SMatch(2, 2, :)));


end