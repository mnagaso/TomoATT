#ifndef ONE_D_MODELS_H
#define ONE_D_MODELS_H

#include <vector>
#include "config.h"


// ak135  http://rses.anu.edu.au/seismology/ak135/ak135t.html
// depth(km), P velocity(km/s)
inline const std::vector<std::vector<CUSTOMREAL>> model_1d_ak135 \
    {{     0.000,      5.8000},
     {    20.000,      5.8000},
     {    20.000,      6.5000},
     {    35.000,      6.5000},
     {    35.000,      8.0400},
     {    77.500,      8.0450},
     {   120.000,      8.0500},
     {   165.000,      8.1750},
     {   210.000,      8.3000},
     {   210.000,      8.3000},
     {   260.000,      8.4825},
     {   310.000,      8.6650},
     {   360.000,      8.8475},
     {   410.000,      9.0300},
     {   410.000,      9.3600},
     {   460.000,      9.5280},
     {   510.000,      9.6960},
     {   560.000,      9.8640},
     {   610.000,     10.0320},
     {   660.000,     10.2000},
     {   660.000,     10.7900},
     {   710.000,     10.9229},
     {   760.000,     11.0558},
     {   809.500,     11.1353},
     {   859.000,     11.2221},
     {   908.500,     11.3068},
     {   958.000,     11.3896},
     {  1007.500,     11.4705},
     {  1057.000,     11.5495},
     {  1106.500,     11.6269},
     {  1156.000,     11.7026},
     {  1205.500,     11.7766},
     {  1255.000,     11.8491},
     {  1304.500,     11.9200},
     {  1354.000,     11.9895},
     {  1403.500,     12.0577},
     {  1453.000,     12.1245},
     {  1502.500,     12.1912},
     {  1552.000,     12.2550},
     {  1601.500,     12.3185},
     {  1651.000,     12.3819},
     {  1700.500,     12.4426},
     {  1750.000,     12.5031},
     {  1799.500,     12.5631},
     {  1849.000,     12.6221},
     {  1898.500,     12.6804},
     {  1948.000,     12.7382},
     {  1997.500,     12.7956},
     {  2047.000,     12.8526},
     {  2096.500,     12.9096},
     {  2146.000,     12.9668},
     {  2195.500,     13.0222},
     {  2245.000,     13.0783},
     {  2294.500,     13.1336},
     {  2344.000,     13.1894},
     {  2393.500,     13.2465},
     {  2443.000,     13.3018},
     {  2492.500,     13.3585},
     {  2542.000,     13.4156},
     {  2591.500,     13.4741},
     {  2640.000,     13.5312},
     {  2690.000,     13.5900},
     {  2740.000,     13.6494},
     {  2740.000,     13.6494},
     {  2789.670,     13.6530},
     {  2839.330,     13.6566},
     {  2891.500,     13.6602},
     {  2891.500,      8.0000},
     {  2939.330,      8.0382},
     {  2989.660,      8.1283},
     {  3039.990,      8.2213},
     {  3090.320,      8.3122},
     {  3140.660,      8.4001},
     {  3190.990,      8.4861},
     {  3241.320,      8.5692},
     {  3291.650,      8.6496},
     {  3341.980,      8.7283},
     {  3392.310,      8.8036},
     {  3442.640,      8.8761},
     {  3492.970,      8.9461},
     {  3543.300,      9.0138},
     {  3593.640,      9.0792},
     {  3643.970,      9.1426},
     {  3694.300,      9.2042},
     {  3744.630,      9.2634},
     {  3794.960,      9.3205},
     {  3845.290,      9.3760},
     {  3895.620,      9.4297},
     {  3945.950,      9.4814},
     {  3996.280,      9.5306},
     {  4046.620,      9.5777},
     {  4096.950,      9.6232},
     {  4147.280,      9.6673},
     {  4197.610,      9.7100},
     {  4247.940,      9.7513},
     {  4298.270,      9.7914},
     {  4348.600,      9.8304},
     {  4398.930,      9.8682},
     {  4449.260,      9.9051},
     {  4499.600,      9.9410},
     {  4549.930,      9.9761},
     {  4600.260,     10.0103},
     {  4650.590,     10.0439},
     {  4700.920,     10.0768},
     {  4801.580,     10.1415},
     {  4851.910,     10.1739},
     {  4902.240,     10.2049},
     {  4952.580,     10.2329},
     {  5002.910,     10.2565},
     {  5053.240,     10.2745},
     {  5103.570,     10.2854},
     {  5153.500,     10.2890},
     {  5153.500,     11.0427},
     {  5204.610,     11.0585},
     {  5255.320,     11.0718},
     {  5306.040,     11.0850},
     {  5356.750,     11.0983},
     {  5407.460,     11.1166},
     {  5458.170,     11.1316},
     {  5508.890,     11.1457},
     {  5559.600,     11.1590},
     {  5610.310,     11.1715},
     {  5661.020,     11.1832},
     {  5711.740,     11.1941},
     {  5813.160,     11.2134},
     {  5863.870,     11.2219},
     {  5914.590,     11.2295},
     {  5965.300,     11.2364},
     {  6016.010,     11.2424},
     {  6066.720,     11.2477},
     {  6117.440,     11.2521},
     {  6168.150,     11.2557},
     {  6218.860,     11.2586},
     {  6269.570,     11.2606},
     {  6320.290,     11.2618},
     {  6371.000,     11.2622}};


// iasp91 model from http://ds.iris.edu/spud/earthmodel/9991809
inline const std::vector<std::vector<CUSTOMREAL>> model_1d_iasp91 \
    {{0.00,     5.8000},
     {1.00,     5.8000},
     {2.00,     5.8000},
     {3.00,     5.8000},
     {4.00,     5.8000},
     {5.00,     5.8000},
     {6.00,     5.8000},
     {7.00,     5.8000},
     {8.00,     5.8000},
     {9.00,     5.8000},
     {10.00,    5.8000},
     {11.00,    5.8000},
     {12.00,    5.8000},
     {13.00,    5.8000},
     {14.00,    5.8000},
     {15.00,    5.8000},
     {16.00,    5.8000},
     {17.00,    5.8000},
     {18.00,    5.8000},
     {19.00,    5.8000},
     {20.00,    5.8000},
     {20.00,    6.5000},
     {21.00,    6.5000},
     {22.00,    6.5000},
     {23.00,    6.5000},
     {24.00,    6.5000},
     {25.00,    6.5000},
     {26.00,    6.5000},
     {27.00,    6.5000},
     {28.00,    6.5000},
     {29.00,    6.5000},
     {30.00,    6.5000},
     {31.00,    6.5000},
     {32.00,    6.5000},
     {33.00,    6.5000},
     {34.00,    6.5000},
     {35.00,    6.5000},
     {35.00,    8.0400},
     {40.00,    8.0406},
     {45.00,    8.0412},
     {50.00,    8.0418},
     {60.00,    8.0429},
     {70.00,    8.0441},
     {80.00,    8.0453},
     {90.00,    8.0465},
     {100.00,   8.0476},
     {110.00,   8.0488},
     {120.00,   8.0500},
     {120.00,   8.0500},
     {130.00,   8.0778},
     {140.00,   8.1056},
     {150.00,   8.1333},
     {160.00,   8.1611},
     {170.00,   8.1889},
     {180.00,   8.2167},
     {190.00,   8.2444},
     {200.00,   8.2722},
     {210.00,   8.3000},
     {210.00,   8.3000},
     {220.00,   8.3365},
     {230.00,   8.3730},
     {240.00,   8.4095},
     {250.00,   8.4460},
     {260.00,   8.4825},
     {270.00,   8.5190},
     {280.00,   8.5555},
     {290.00,   8.5920},
     {300.00,   8.6285},
     {310.00,   8.6650},
     {320.00,   8.7015},
     {330.00,   8.7380},
     {340.00,   8.7745},
     {350.00,   8.8110},
     {360.00,   8.8475},
     {370.00,   8.8840},
     {380.00,   8.9205},
     {390.00,   8.9570},
     {400.00,   8.9935},
     {410.00,   9.0300},
     {410.00,   9.3600},
     {420.00,   9.3936},
     {430.00,   9.4272},
     {440.00,   9.4608},
     {450.00,   9.4944},
     {460.00,   9.5280},
     {470.00,   9.5616},
     {480.00,   9.5952},
     {490.00,   9.6288},
     {500.00,   9.6624},
     {510.00,   9.6960},
     {520.00,   9.7296},
     {530.00,   9.7632},
     {540.00,   9.7968},
     {550.00,   9.8304},
     {560.00,   9.8640},
     {570.00,   9.8976},
     {580.00,   9.9312},
     {590.00,   9.9648},
     {600.00,   9.9984},
     {610.00,   10.0320},
     {620.00,   10.0656},
     {630.00,   10.0992},
     {640.00,   10.1328},
     {650.00,   10.1664},
     {660.00,   10.2000},
     {660.00,   10.7900},
     {670.00,   10.8166},
     {680.00,   10.8432},
     {690.00,   10.8697},
     {700.00,   10.8963},
     {710.00,   10.9229},
     {720.00,   10.9495},
     {730.00,   10.9761},
     {740.00,   11.0026},
     {750.00,   11.0292},
     {760.00,   11.0558},
     {760.00,   11.0558},
     {770.00,   11.0738},
     {780.00,   11.0917},
     {790.00,   11.1095},
     {800.00,   11.1272},
     {900.00,   11.2997},
     {1000.00,  11.4640},
     {1100.00,  11.6208},
     {1200.00,  11.7707},
     {1300.00,  11.9142},
     {1400.00,  12.0521},
     {1500.00,  12.1849},
     {2000.00,  12.7944},
     {2500.00,  13.3697},
     {2700.00,  13.6076},
     {2740.00,  13.6564},
     {2740.00,  13.6564},
     {2750.00,  13.6587},
     {2800.00,  13.6703},
     {2850.00,  13.6818},
     {2889.00,  13.6908},
     {2889.00,  8.0088},
     {2900.00,  8.0280},
     {3000.00,  8.1995},
     {3100.00,  8.3642},
     {3200.00,  8.5222},
     {3300.00,  8.6735},
     {3400.00,  8.8180},
     {3500.00,  8.9558},
     {4000.00,  9.5437},
     {4500.00,  9.9633},
     {5153.90,  10.2578},
     {5153.90,  11.0914},
     {5500.00,  11.1644},
     {6000.00,  11.2270},
     {6371.00,  11.2409}};


// fortran test model
inline const std::vector<std::vector<CUSTOMREAL>> model_1d_prem \
    = {{    0.00,     5.80000},
       {   15.00,     5.80000},
       {   15.00,     6.80000},
       {   24.40,     6.80000},
       {   24.40,     8.11061},
       {   40.00,     8.10119},
       {   60.00,     8.08907},
       {   80.00,     8.07688},
       {  115.00,     8.05540},
       {  150.00,     8.03370},
       {  185.00,     8.01180},
       {  220.00,     7.98970},
       {  220.00,     8.55896},
       {  265.00,     8.64552},
       {  310.00,     8.73209},
       {  355.00,     8.81867},
       {  400.00,     8.90522},
       {  400.00,     9.13397},
       {  450.00,     9.38990},
       {  500.00,     9.64588},
       {  550.00,     9.90185},
       {  600.00,    10.15782},
       {  635.00,    10.21203},
       {  670.00,    10.26622},
       {  670.00,    10.75131},
       {  721.00,    10.91005},
       {  771.00,    11.06557},
       {  871.00,    11.24490},
       {  971.00,    11.41560},
       { 1071.00,    11.57828},
       { 1171.00,    11.73357},
       { 1271.00,    11.88209},
       { 1371.00,    12.02445},
       { 1471.00,    12.16126},
       { 1571.00,    12.29316},
       { 1671.00,    12.42075},
       { 1771.00,    12.54466},
       { 1871.00,    12.66550},
       { 1971.00,    12.78389},
       { 2071.00,    12.90045},
       { 2171.00,    13.01579},
       { 2271.00,    13.13055},
       { 2371.00,    13.24532},
       { 2471.00,    13.36074},
       { 2571.00,    13.47742},
       { 2671.00,    13.59597},
       { 2741.00,    13.68041},
       { 2771.00,    13.68753},
       { 2871.00,    13.71168},
       { 2891.00,    13.71660},
       { 2891.00,     8.06482},
       { 2971.00,     8.19939},
       { 3071.00,     8.36019},
       { 3171.00,     8.51298},
       { 3271.00,     8.65805},
       { 3371.00,     8.79573},
       { 3471.00,     8.92632},
       { 3571.00,     9.05015},
       { 3671.00,     9.16752},
       { 3771.00,     9.27867},
       { 3871.00,     9.38418},
       { 3971.00,     9.48409},
       { 4071.00,     9.57881},
       { 4171.00,     9.66865},
       { 4271.00,     9.75393},
       { 4371.00,     9.83496},
       { 4471.00,     9.91206},
       { 4571.00,     9.98554},
       { 4671.00,    10.05572},
       { 4771.00,    10.12291},
       { 4871.00,    10.18743},
       { 4971.00,    10.24959},
       { 5071.00,    10.30971},
       { 5149.50,    10.35568},
       { 5149.50,    11.02827},
       { 5171.00,    11.03643},
       { 5271.00,    11.07249},
       { 5371.00,    11.10542},
       { 5471.00,    11.13521},
       { 5571.00,    11.16186},
       { 5671.00,    11.18538},
       { 5771.00,    11.20576},
       { 5871.00,    11.22301},
       { 5971.00,    11.23712},
       { 6071.00,    11.24809},
       { 6171.00,    11.25593},
       { 6271.00,    11.26064},
       { 6371.00,    11.26220}};


#endif