clear all; close all; clc
%% Define shape functions and get stiffness matrix components
% symbolic math variables
x = sym('x','real');
xi = sym('xi','real');
L = sym('l','real');
E_p = sym('E_p','real');
A_p = sym('A_p','real');
I_p = sym('I_p','real');
E = sym('E','real');
A = sym('A','real');
I = sym('I','real');
rho = sym('rho','real');
rho_p = sym('rho_p','real');
F_y = sym('F_x', 'real');
eps_S = sym('eps_S', 'real');
L_b = sym('L_b', 'real');
w_b = sym('w_b', 'real');
J_p = sym('J_p', 'real');
R = sym('R', 'real');
b_p = sym('b_p', 'real');
h_p = sym('h_p', 'real');

%shape functions
N1 = 1  -(3*((x/L)^2)) + 2*((x/L)^3); 
N2 = L*(x/L - 2*(x/L)^2 + (x/L)^3); 
N3 = 3*(x/L)^2 - 2*(x/L)^3; 
N4 = L*(-(x/L)^2 + (x/L)^3);

N = [N1 N2 N3 N4];


BQ = [1/(w_b*L)];
%differentiateLshape functions wrt isoparam. coordinates
dN1_dxi=gradient(N1,x);
dN2_dxi=gradient(N2,x);
dN3_dxi=gradient(N3,x);
dN4_dxi=gradient(N4,x);
dN1_dxi_2=gradient(dN1_dxi,x);
dN2_dxi_2=gradient(dN2_dxi,x);
dN3_dxi_2=gradient(dN3_dxi,x);
dN4_dxi_2=gradient(dN4_dxi,x);

%second derivatives
Be_2 = [dN1_dxi_2 dN2_dxi_2 dN3_dxi_2 dN4_dxi_2];

%element mass matrix
M_e = rho*A*int(N'*N,x,0,L) + rho_p*A_p*int(N'*N,x,0,L);

M_e_beam = rho*A*int(N'*N,x,0,L);

K_e_beam = E*I*int(Be_2'*Be_2,x,0,L);
%element external load vector
f_ext = int(N'*F_y,x, 0,L);

K_vv = (E_p*I_p*int(Be_2'*Be_2,0,L) + E*I*int(Be_2'*Be_2,0,L));
K_QQ = A_p/(eps_S*w_b^2*L);
K_Qv = -b_p*J_p*int(BQ'*Be_2,x,0,L);
K_vQ = -b_p*J_p*int(Be_2'*BQ,x,0,L);


%Generate matlab functions from symbolic math expressions
matlabFunction(M_e_beam,'File','M_e_beam_piezo_cst_func','Vars',{rho,A,L});
matlabFunction(K_e_beam,'File','K_e_beam_piezo_cst_func','Vars',{I,E,L});
matlabFunction(M_e,'File','Me_piezo_cst_func','Vars',{rho,A,L,rho_p,A_p});
matlabFunction(f_ext,'File','f_ext_piezo_cst_func','Vars',{F_y,L});
matlabFunction(K_vv,'File','K_vv_piezo_cst_func','Vars',{I,E,L,E_p,I_p});
matlabFunction(K_QQ,'File','K_QQ_piezo_cst_func','Vars',{A_p,w_b,eps_S,L});
matlabFunction(K_Qv,'File','K_Qv_piezo_cst_func','Vars',{b_p,w_b,J_p,L});
matlabFunction(K_vQ,'File','K_vQ_piezo_cst_func','Vars',{b_p,w_b,J_p,L});
