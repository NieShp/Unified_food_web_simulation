function dy = SimpleFoodChain(B,par)
x=B(1);  % plant
y=B(2);  % herbivore
Q = (par.N - y*par.q)/x;
K = (par.N - y*par.q)/par.Q0;
dxdt = par.r * (1 - x/K) * x - par.a*x/(1+par.a*par.h*x) * y;
dydt = par.e * min(Q/par.q, 1) * par.a*x/(1+par.a*par.h*x) * y - par.m *y;
dy=[dxdt; dydt];


