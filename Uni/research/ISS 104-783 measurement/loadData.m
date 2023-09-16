function [freq, Port1Meas, Port2Meas, SThruMeas, SLineMeas] = loadData(folderName)

[freq_DC_to_67_GHz, SOpen_DC_to_67_GHz] = SXPParse(fullfile(folderName, 'Open_DC_to_67_GHz.s2p'));
[freq_67_to_110_GHz, SOpen_67_to_110_GHz] = SXPParse(fullfile(folderName, 'Open_67_to_110_GHz.s2p'));
freq = cat(1, freq_DC_to_67_GHz, freq_67_to_110_GHz);
SOpen = cat(3, SOpen_DC_to_67_GHz, SOpen_67_to_110_GHz);

[~, SShort_DC_to_67_GHz] = SXPParse(fullfile(folderName, 'Short_A1_DC_to_67_GHz.s2p'));
[~, SShort_67_to_110_GHz] = SXPParse(fullfile(folderName, 'Short_A2_67_to_110_GHz.s2p'));
SShort = cat(3, SShort_DC_to_67_GHz, SShort_67_to_110_GHz);

[~, SMatch_DC_to_67_GHz] = SXPParse(fullfile(folderName, 'Match_A1_DC_to_67_GHz.s2p'));
[~, SMatch_67_to_110_GHz] = SXPParse(fullfile(folderName, 'Match_A2_67_to_110_GHz.s2p'));
SMatch = cat(3, SMatch_DC_to_67_GHz, SMatch_67_to_110_GHz);

[~, SThruMeas_DC_to_67_GHz] = SXPParse(fullfile(folderName, 'Thru_A1_DC_to_67_GHz.s2p'));
[~, SThruMeas_67_to_110_GHz] = SXPParse(fullfile(folderName, 'Thru_A2_67_to_110_GHz.s2p'));
SThruMeas = cat(3, SThruMeas_DC_to_67_GHz, SThruMeas_67_to_110_GHz);

[~, SLine1Meas_DC_to_67_GHz] = SXPParse(fullfile(folderName, 'Line_A1_DC_to_67_GHz.s2p'));
[~, SLine1Meas_67_to_110_GHz] = SXPParse(fullfile(folderName, 'Line_A1_67_to_110_GHz.s2p'));
SLine1Meas = cat(3, SLine1Meas_DC_to_67_GHz, SLine1Meas_67_to_110_GHz);

[~, SLine2Meas_DC_to_67_GHz] = SXPParse(fullfile(folderName, 'Line_A2_DC_to_67_GHz.s2p'));
[~, SLine2Meas_67_to_110_GHz] = SXPParse(fullfile(folderName, 'Line_A2_67_to_110_GHz.s2p'));
SLine2Meas = cat(3, SLine2Meas_DC_to_67_GHz, SLine2Meas_67_to_110_GHz);

[~, SLine3Meas_DC_to_67_GHz] = SXPParse(fullfile(folderName, 'Line_A3_DC_to_67_GHz.s2p'));
[~, SLine3Meas_67_to_110_GHz] = SXPParse(fullfile(folderName, 'Line_A3_67_to_110_GHz.s2p'));
SLine3Meas = cat(3, SLine3Meas_DC_to_67_GHz, SLine3Meas_67_to_110_GHz);

[~, SLine4Meas_DC_to_67_GHz] = SXPParse(fullfile(folderName, 'Line_A4_DC_to_67_GHz.s2p'));
[~, SLine4Meas_67_to_110_GHz] = SXPParse(fullfile(folderName, 'Line_A4_67_to_110_GHz.s2p'));
SLine4Meas = cat(3, SLine4Meas_DC_to_67_GHz, SLine4Meas_67_to_110_GHz);

[~, SLine5Meas_DC_to_67_GHz] = SXPParse(fullfile(folderName, 'Line_A5_DC_to_67_GHz.s2p'));
[~, SLine5Meas_67_to_110_GHz] = SXPParse(fullfile(folderName, 'Line_A5_67_to_110_GHz.s2p'));
SLine5Meas = cat(3, SLine5Meas_DC_to_67_GHz, SLine5Meas_67_to_110_GHz);

SLineMeas = cat(4, SLine1Meas, SLine2Meas, SLine3Meas, SLine4Meas, SLine5Meas);

Port1Meas = struct('open', squeeze(SOpen(1, 1, :)), ...
    'short', squeeze(SShort(1, 1, :)), ...
    'match', squeeze(SMatch(1, 1, :)));

Port2Meas = struct('open', squeeze(SOpen(2, 2, :)), ...
    'short', squeeze(SShort(2, 2, :)), ...
    'match', squeeze(SMatch(2, 2, :)));


end