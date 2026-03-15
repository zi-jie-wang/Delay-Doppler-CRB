f = 5.89e9;
B = 10e6;
M = 64;
Mcp = 0.25*M;
N = 50;
Q = M*N;
Ts = 1/B;
Tsym = (M+Mcp)*Ts;
Tcpi = N*Tsym;

P = 4;
tau = [0,8,4,6]*Tsym/M;
nu = [4.82,-3.23,1.38,-2.47]/(N*Tsym);
g = [-0.02-0.09j, 0.40 + 0.73j, 0.03 + 0.45j, 0.15-0.43j];

noise_var = 1;
pt = 10^(1);