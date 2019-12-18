static double emu_tinker_bias_y[40][4] = {
    {0.22772415000000001,-0.37692982250000007,-0.0402944999999999,0.04022747500000001},
    {-0.28278585,0.19515017749999994,0.04740549999999999,-0.06182552499999994},
    {0.3212841500000001,-0.8808398225,-0.11989450000000001,0.06898547500000007},
    {0.17846414999999993,-0.6810798225000001,-0.09233449999999999,0.05497347500000005},
    {0.4098541499999999,-0.5539598225000001,-0.07733450000000008,0.0647434750000001},
    {-0.76878285,0.5381131774999999,0.08117550000000007,-0.06465652499999996},
    {0.40140415000000007,-0.5459038225,-0.08968449999999994,0.08774847500000005},
    {0.24942415000000007,0.026990177499999934,-0.019294499999999992,0.03807147500000008},
    {-0.39614284999999994,0.11911017749999997,0.04924550000000005,-0.06747352499999992},
    {0.3164341500000001,-0.16069982250000003,-0.03416449999999993,0.033380475000000076},
    {-0.2499258499999999,0.5025871774999999,0.06269550000000002,-0.06842252499999996},
    {-0.12223585000000003,-0.3533408225000001,-0.03049450000000009,0.017900475000000027},
    {0.2607141500000001,-0.11666582250000002,-0.012854499999999991,0.03635147500000002},
    {0.08741415000000008,0.26517217749999994,0.04785550000000005,-0.043800524999999924},
    {-0.76097985,0.24849317749999994,0.05062549999999999,-0.022498524999999936},
    {-0.06336585000000006,0.5376801775,0.07865549999999999,-0.07991652499999996},
    {-0.09263584999999996,0.46694977749999994,0.041795500000000096,-0.038999524999999924},
    {-0.06395584999999993,0.5896181775,0.06771550000000004,-0.06875352499999993},
    {0.37515414999999996,0.3734211774999999,0.02701550000000008,0.015070475000000028},
    {0.37989414999999993,0.42874727749999997,0.04378550000000003,-0.012818524999999914},
    {0.04260415000000006,-0.06311582250000008,0.01377549999999994,-0.017001524999999962},
    {0.15116415000000005,-0.18857682250000007,-0.031214499999999923,0.028068475000000093},
    {-0.46622184999999994,-0.00875982250000007,0.025565499999999908,-0.020490524999999926},
    {0.9793741499999999,-1.4424998225,-0.23205450000000005,0.1436824750000001},
    {0.5773441500000001,-0.40441882250000005,-0.06932450000000001,0.07315647500000011},
    {-0.0900458500000001,0.05627017749999996,-0.0025844999999999896,0.0005704750000000702},
    {-0.43227585,0.15328717749999995,0.048885499999999915,-0.061249524999999916},
    {0.06163414999999994,0.17322317749999994,0.002195500000000017,0.002554475000000056},
    {-1.11435885,0.6018431774999999,0.10177550000000002,-0.051648524999999945},
    {-0.58425185,0.33632737749999997,0.06437549999999992,-0.06131052499999995},
    {-0.05438585000000007,0.5742031774999999,0.06327550000000004,-0.05677052499999996},
    {-0.46905785,0.5551571774999999,0.08878549999999996,-0.09477352499999991},
    {0.48297415,0.06671617749999992,-0.005604500000000012,0.03804847500000008},
    {0.10267414999999991,-0.12002582250000005,-0.021424499999999957,0.03506147500000001},
    {-0.2600958499999999,0.08162417749999995,0.021935500000000108,-0.03365652499999994},
    {0.7054041500000001,-1.0345198225,-0.15660450000000004,0.13090447500000002},
    {-0.49813684999999996,0.16804117749999994,0.04431549999999995,-0.04419252499999993},
    {-0.10155585,-0.29502082250000006,-0.039974500000000024,0.02313947500000002},
    {0.2226541500000001,0.20897417749999994,-0.0015045000000000197,0.01563147500000006},
    {0.33760415,-0.04134282250000004,0.003785499999999997,0.021988475000000007}
};
static double emu_tinker_bias_ymean[4] = {
    1.30796585, -0.39295018, -1.1741355 ,  0.52458352
};

static double emu_tinker_bias_cosmopars[40][7] = {
    {0.0226832,0.11406,-0.816597,0.975589,3.09292,63.3657,2.91875},
    {0.0224507,0.117258,-1.13444,0.976523,3.14987,73.0997,3.17375},
    {0.0229982,0.108749,-0.684659,0.99745,3.09377,63.7085,3.25875},
    {0.0226759,0.112306,-0.743738,0.948081,3.0009,64.0419,3.55625},
    {0.0220846,0.106252,-0.766605,0.965093,3.11856,65.0491,2.66375},
    {0.0206646,0.129503,-1.32636,0.927846,3.02365,72.7542,2.96125},
    {0.0229415,0.111503,-0.710308,0.970642,3.01597,62.6987,2.70625},
    {0.0228137,0.119597,-0.866603,0.966298,3.16214,64.3716,3.93875},
    {0.0206769,0.123835,-1.16373,0.949131,3.14713,69.4,3.59875},
    {0.0213439,0.115813,-0.830933,0.947523,3.07177,62.3564,3.89625},
    {0.0218933,0.128951,-1.24115,0.961017,3.05046,72.086,4.23625},
    {0.0226212,0.108951,-0.860928,0.996032,3.15811,67.734,2.83375},
    {0.022506,0.116766,-0.879472,0.953951,3.04801,65.3799,2.87625},
    {0.0218955,0.117176,-1.1195,0.978815,3.06781,71.0834,3.00375},
    {0.0226491,0.127108,-1.11683,0.972365,3.09389,68.73,2.74875},
    {0.0214676,0.128502,-1.30291,0.933567,3.0943,74.1,3.72625},
    {0.0217741,0.120676,-1.1313,0.966216,3.01408,70.0743,3.76875},
    {0.0222812,0.119424,-1.24807,0.951972,3.03502,74.4391,3.21625},
    {0.0228741,0.11574,-1.03218,0.95331,3.0202,70.7516,4.27875},
    {0.0223826,0.113263,-1.09161,0.967291,3.09569,72.4296,3.68375},
    {0.0222537,0.122488,-0.989837,0.952888,3.11987,67.0552,3.38625},
    {0.0236151,0.117169,-0.865832,0.975785,3.13158,66.389,3.85375},
    {0.0215385,0.121009,-1.03161,0.958557,3.07232,68.0619,2.62125},
    {0.0227263,0.101218,-0.565849,0.974575,3.0192,62.0335,3.47125},
    {0.0224536,0.110277,-0.761456,0.958918,3.14394,63.0321,4.15125},
    {0.0209291,0.117057,-0.947878,0.934464,3.03679,65.7127,3.08875},
    {0.0223621,0.119181,-1.1251,0.94426,3.12789,71.7554,2.79125},
    {0.021358,0.113425,-0.964588,0.966388,3.01471,67.393,4.02375},
    {0.0217472,0.131777,-1.39992,0.95857,3.14656,74.7675,3.81125},
    {0.0223248,0.128915,-1.23607,0.940146,3.15905,71.4136,3.42875},
    {0.0219279,0.123914,-1.22368,0.955236,3.11828,73.4315,4.06625},
    {0.0212396,0.127579,-1.38165,0.956072,3.07641,73.76,3.34375},
    {0.0224969,0.112845,-0.925758,0.949484,3.04281,68.4034,3.98125},
    {0.0233732,0.11497,-0.875221,0.989173,3.14916,66.0523,3.64125},
    {0.0227557,0.122229,-1.03229,0.950028,3.10694,69.0704,3.13125},
    {0.0234108,0.107641,-0.61282,0.995568,3.14029,61.6947,3.04625},
    {0.0219989,0.121297,-1.10822,0.967361,3.17942,70.4113,3.30125},
    {0.0228745,0.109744,-0.848746,0.977623,3.07229,66.7262,3.51375},
    {0.0237124,0.114975,-0.955008,0.976591,3.05351,69.7468,4.10875},
    {0.0217407,0.12007,-0.941053,0.960207,3.0927,64.7043,4.19375}
};

static double emu_tinker_bias_variances[40][4] = {
    {0.0072293,0.00439722,2.81818e-06,4.74412e-05},
    {0.00699877,0.00269795,7.42337e-07,2.39677e-05},
    {0.0110443,0.00711122,4.65542e-06,0.000112801},
    {0.0106216,0.0102302,4.31172e-06,0.000103383},
    {0.00724683,0.00524066,4.31665e-06,6.43917e-05},
    {0.0106784,0.00234987,5.0869e-06,6.98768e-05},
    {0.00854694,0.00664262,5.57295e-06,8.14472e-05},
    {0.00686848,0.0047873,3.72477e-06,5.69853e-05},
    {0.00724118,0.00248583,8.77754e-07,2.31016e-05},
    {0.00877297,0.00692596,6.01501e-06,9.17487e-05},
    {0.00627794,0.00336004,9.26898e-07,2.86192e-05},
    {0.00788265,0.00320692,8.60393e-07,2.75182e-05},
    {0.00645664,0.00441144,3.48456e-06,5.00461e-05},
    {0.00516963,0.0030301,1.88609e-06,2.83136e-05},
    {0.00627674,0.002754,3.41576e-06,5.52287e-05},
    {0.00476999,0.00279074,1.30704e-06,2.45841e-05},
    {0.00625622,0.00428272,1.611e-06,3.75563e-05},
    {0.00486567,0.00359014,1.30515e-06,2.83193e-05},
    {0.00793602,0.00820381,8.56658e-06,0.000112397},
    {0.00547263,0.00462295,5.10187e-06,6.12737e-05},
    {0.00639258,0.00334604,1.60622e-06,3.16773e-05},
    {0.00728083,0.00523552,2.64198e-06,5.18353e-05},
    {0.0114682,0.00289903,2.18067e-06,3.56099e-05},
    {0.0289611,0.0255266,5.05678e-05,0.000513587},
    {0.0115307,0.00948345,1.29089e-05,0.000154464},
    {0.00693994,0.00459213,1.46837e-06,3.85683e-05},
    {0.00863687,0.00262847,1.12103e-06,2.78236e-05},
    {0.00777637,0.00787936,3.47484e-06,7.01753e-05},
    {0.00742029,0.00223942,1.12076e-05,0.000137364},
    {0.00756381,0.00224989,1.71963e-06,3.86059e-05},
    {0.00491035,0.00316571,1.4925e-06,3.04492e-05},
    {0.00734602,0.00208953,8.32388e-07,2.69424e-05},
    {0.00901434,0.00916022,1.04109e-05,0.000133636},
    {0.0066438,0.00432074,2.23461e-06,4.61122e-05},
    {0.00686723,0.00277957,7.85587e-07,2.64958e-05},
    {0.0119653,0.00774551,1.07802e-05,0.000130799},
    {0.00885545,0.00250722,1.62889e-06,3.42051e-05},
    {0.00961551,0.00626229,1.51812e-06,4.70177e-05},
    {0.00691408,0.00685252,4.84129e-06,7.6411e-05},
    {0.00646323,0.00500426,4.99086e-06,6.97103e-05}
};

static double emu_tinker_bias_rotation_matrix[4][4] = {
    {0.10646337980893085,0.021132636136643016,-0.97064255773295,-0.21464385751102177},
    {-0.978398102331227,-0.14405378008966646,-0.08332901175552901,-0.12264557715783554},
    {0.168930091936486,-0.9427314540214435,0.060196305393736965,-0.28124088332582375},
    {-0.053528741369879666,-0.3000988599285399,-0.21745276667457572,0.9272484253869187}
};


static double emu_tinker_bias_log_lambda[4][7] = {

    {-7.447915268014326,-7.464755353609271,-2.599556684860801,-5.8827162885183135,-3.4939080354809984,22.358830230679345,2.032157502031441},
    {-8.23942986061086,13.078292333315005,-1.6453454526931253,24.55006207730737,-2.208502981716542,6.529841736884682,2.702316841471984},
    {-2.3990180929022267,-4.984895914479709,0.6207751414355313,24.725128904515856,6.6630793421725425,6.460552261861607,5.096751130041617},
    {-4.200230133222608,-3.9710113322602116,2.102438559177407,23.3325814549079,5.148802333913723,10.75333618297953,9.909683276181012},

};

typedef struct {
    int tinker_bias_ncosmo;
    int tinker_bias_nparam;
    int tinker_bias_nsamp;
    int tinker_nparam;
} TinkerEmuParam;

TinkerEmuParam tinkerEmuParam = {
    .tinker_bias_ncosmo=7, //do not change
    .tinker_bias_nparam=4, // do not change
    .tinker_bias_nsamp=40, 
    .tinker_nparam = 6 // number of parameters of equation 7 in 1907.13167.; 
};

