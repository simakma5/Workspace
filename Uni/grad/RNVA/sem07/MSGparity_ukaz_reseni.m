A =[1 0 1 0 0 1; %D29*
    0 1 0 1 1 0; %D30*
    1 0 1 0 1 0; %d1
    1 1 0 1 0 0;
    1 1 1 0 1 1;
    0 1 1 1 0 0;
    1 0 1 1 1 1; %d5
    1 1 0 1 1 1;
    0 1 1 0 1 0;
    0 0 1 1 0 1;
    0 0 0 1 1 1;
    1 0 0 0 1 1; %d10
    1 1 0 0 0 1;
    1 1 1 0 0 0;
    1 1 1 1 0 1;
    1 1 1 1 1 0;
    0 1 1 1 1 1; %d15
    0 0 1 1 1 0;
    1 0 0 1 1 0;
    1 1 0 0 1 0;
    0 1 1 0 0 1;
    1 0 1 1 0 0; %d20
    0 1 0 1 1 0;
    0 0 1 0 1 1;
    1 0 0 1 0 1;
    0 1 0 0 1 1]; %d24

HT = [A;eye(6)];
preamble = [1 0 0 0 1 0 1 1];

% skupina 101 11:00
%       1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128	129	130	131	132	133	134	135	136	137	138	139	140	141	142	143	144	145	146	147	148	149	150	151	152	153	154	155	156	157	158	159	160	161	162	163	164	165	166	167	168	169	170	171	172	173	174	175	176	177	178	179	180	181	182	183	184	185	186	187	188	189	190	191	192	193	194	195	196	197	198	199	200	201
% Yp =  [0	0	0	0	1	0	1	1	0	1	1	0	0	0	0	0	1   0	0	0	1	0	1	1	0	0	0	0	0	1	0	1	0	1	0	1	1	1	0	0	1	0	0	0	0	1	1	0	0	0	1	1	1	1	0	1	0	0	1	0	1	0	1	1	0	1	1	0	0	1	0	0	1	1	0	0	0	1	0	1	0	1	0	0	1	1	0	1	0	0	0	1	0	0	0	0	0	0	0	1	0	1	0	0	1	1	1	0	1	0	0	1	1	1	0	1	1	1	0	0	1	1	1	0	0	0	0	1	1	1	1	0	1	1	0	1	1	1	0	0	0	1	1	0	0	0	0	1	0	1	1	0	1	1	0	1	1	0	0	0	0	1	0	0	0	1	1	0	0	0	0	1	1	0	1	0	0	1	1	0	0	0	1	1	0	1	1	0	1	1	0	0	1	0	0	0	1	1	0	1	0]
% i_preamb = 17 %zacatek preambule pro 101

% S3 =  0     0     1     1     0     1
% 10. bit (d8) 3. slova
% gps_WN = 343
% TOW = 57706
% date = 23-Mar-2006 00:10:22 (následující podrámec)


% skupina 102 12:45
%       1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128	129	130	131	132	133	134	135	136	137	138	139	140	141	142	143	144	145	146	147	148	149	150	151	152	153	154	155	156	157	158	159	160	161	162	163	164	165	166	167	168	169	170	171	172	173	174	175	176	177	178	179	180	181	182	183	184	185	186	187	188	189	190	191	192	193	194	195	196	197	198	199	200	201
Yp =  [0	0	1	1	0	1	1	1	0	0	0	0	0	0	0	0   1	1	1	1	1	1	0	1	1	0	1	0	0	0	1	0	0	0	1	0	1	1	0	0	0	0	0	1	0	1	0	1	0	1	1	1	0	0	1	0	0	0	0	1	1	0	0	1	1	0	0	1	1	1	0	1	0	1	0	0	1	1	0	1	1	0	1	1	0	1	0	1	0	0	0	1	0	1	0	1	0	1	1	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	0	0	0	0	0	1	0	1	1	1	0	1	0	0	1	1	0	1	1	0	1	1	0	1	1	0	1	1	1	1	0	1	0	0	1	1	0	1	1	1	1	0	1	0	1	1	0	1	0	1	0	0	1	1	1	0	1	1	0	0	0	0	0	0	1	0	1	1	1	1	1	1	1	0	0	1	0	0	0	0	1   1]
i_preamb = 31 %zacatek preambule pro 102

% S2 =   1     0     1     1     1     1
% 7. bit (d5) 2. slova
% gps_WN = 343
% TOW = 56406
% date = 22-Mar-2006 22:00:22 (následující podrámec)

preamb_corr = xcorr(Yp.*2-1,preamble.*2-1);
plot(preamb_corr(length(Yp):end))

Y1 = Yp(i_preamb-2:i_preamb+29);    %1. slovo
Y2 = Yp(i_preamb+28:i_preamb+59); %2. slovo
Y3 = Yp(i_preamb+58:i_preamb+89); %3. slovo

Y1(3:26) = mod(Y1(3:26)+Y1(2),2);
S1 = mod(Y1 * HT,2) % vypocet syndromu
Y2(3:26) = mod(Y2(3:26)+Y2(2),2);
S2 = mod(Y2 * HT,2)
Y3(3:26) = mod(Y3(3:26)+Y3(2),2);
S3 = mod(Y3 * HT,2)

E1 = sum((ones(32,1)*S1==HT)')==6; % predpoklad jednonasobne chyby => chybove slovo E
Y1 = xor(Y1, E1);
So1 = mod(Y1 * HT,2)

E2 = sum((ones(32,1)*S2==HT)')==6;
Y2 = xor(Y2, E2);
So2 = mod(Y2 * HT,2)

E3 = sum((ones(32,1)*S3==HT)')==6;
Y3 = xor(Y3, E3);
So3 = mod(Y3 * HT,2)

%  Prevod GPS_WN a TOW na datum a cas
refgpsDT = datenum(1980,1,6);         %  gps reference day
gps_WN = Y3(3:12)*2.^[9:-1:0]'
TOW = Y2(3:19)*2.^[16:-1:0]'
gpsroll = 1; %preteceni cisla tydne v noci 21/22 August 1999
leap_s = 14; %1.1.2006-31.12.2008 dle http://tf.nist.gov/pubs/bulletin/leapsecond.htm
gpsDT  = refgpsDT + (gps_WN + 1024*gpsroll)*7 + (TOW*6-leap_s)/(3600*24);
date = datestr(gpsDT,0)
