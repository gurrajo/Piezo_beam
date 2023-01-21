function f = piezo_func_short(t,y)
%% Material constants and geometry info
n = 5;
m = 1;
V_tot = 16*60*0.26*10^-9;
V_piezo = 12*50*0.11*10^-9;

mass_p = 0.367*10^-3;
C = 5.8*10^-10;
piezo_L = 50*10^-3;
w_p = 12*10^-3;
h_p = 0.11*10^-3;
eps_s = C*h_p/(w_p*piezo_L);
rho_p = mass_p/V_tot;

rho = 1.4*10^3;
E = 3.55*10^9;
E_piezo = 4.2*10^9;
E_PET = 3.14*10^9;
E_p = (V_piezo*E_piezo + E_PET*(V_tot-V_piezo))/V_tot;

h_b = 0.5*10^-3;
w_b = 40*10^-3;
A = h_b*w_b;
h_p = 0.26*10^-3;
A_p = h_p*w_p;
%b_p = 8.87*10^8;
d = 13*10^-10;
b_p = E_p*d/(eps_s);

zeta = 4*10^-2;

J_p = A_p*(h_p/2 + h_b/2);

I = (w_b*h_b^3)/12;
I_p = ((w_p*h_p^3)/12 + A_p*(h_p/2 + h_b/2)^2);

%% Matrices and vectors generated from Piezo_beam.m and Piezo_beam_2.m
% M_p = M_p_new_cst_func(rho_p,A_p);
%K_p = K_p_new_cst_func(I,E,A,rho,rho_p,A_p);
% 
% M_b = M_b_new_cst_func(rho,A);
% K_b = K_b_new_cst_func(I,E,rho,A);
% 
% B_e = B_e_cst_func(eps_s);
% C_a = C_a_cst_func(eps_s,w_b,A_p);
% B_a = B_a_new_cst_func(b_p,J_p,w_b);

M_p = M_p_cst_func(rho_p,A_p);
K_p = K_p_cst_func(I_p,E_p);

M_b = M_b_cst_func(rho,A);
K_b = K_b_cst_func(I,E);

B_e = B_e_cst_func(eps_s);
C_a = C_a_cst_func(eps_s,w_b,A_p);
B_a = B_a_cst_func(b_p,J_p,w_b);


Xi = Xi_new_cst_func(E,I,rho,A,zeta);

K_a = K_b + K_p;
M_a = M_b + M_p;
Omega = inv(M_a)*K_a;

B_1 = inv(M_a)*B_a;
B_2 = inv(M_a)*B_e';

%% External load
% currently neglected
F = 0;
%% Define ODE function

f = zeros(n*2,1);

f(1:2:(n*2)-1) = y(2:2:n*2);
%y(end) = C_a*(B_a')*y(1:2:(n*2)-1); 
f(2:2:n*2) = B_2*F - B_1*C_a*(B_a')*y(1:2:(n*2)-1) - Xi*y(2:2:n*2) - Omega*y(1:2:(n*2)-1);

end