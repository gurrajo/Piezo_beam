clear all;close all; clc

n = 5; % number of modes
m = 1; % number of transducers

xi = sym('xi', [n,1], 'real');
x = sym('x', 'real');

E_p = sym('E_p',[m,1],  'real');
I_p = sym('I_p',[m,1], 'real');
rho_p = sym('rho_p',[m,1], 'real');
A_p = sym('A_p',[m,1], 'real');


E = sym('E', 'real');
I = sym('I', 'real');
rho = sym('rho', 'real');
A = sym('A', 'real');
zeta = sym('zeta','real');

x_trans = [0.02, 0.07];
L_b = 0.13;

C_n = [3.5160, 22.0345, 61.6972, 120.0902, 199.86];
a = sqrt(C_n)/L_b;
sig_n = [0.734096, 1.018466, 0.999225, 1.000033, 1];

omega_n = C_n*sqrt(E*I/(rho*A*L_b^4));


Y_k = (cosh(a*x) - cos(a*x))-(sig_n.*(sinh(a*x) - sin(a*x)));
for i = 1:n
    d_Y_k(i) = gradient(Y_k(i), x);
    d2_Y_k(i) = gradient(d_Y_k(i), x);
end

M_b = sym('M_b',[n, n], 'real');
M_b(:,:) = zeros(n, n);
K_b = sym('K_b',[n, n], 'real');
K_b(:,:) = zeros(n, n);
Xi = sym('Xi',[n, n], 'real');
Xi(:,:) = zeros(n, n);
for i = 1:n
    Xi(i,i) = zeta*omega_n(i);
    M_b(i,i) = rho*A*int(Y_k(i).*Y_k(i),x, 0, L_b);
    temp = int(d2_Y_k(i).*d2_Y_k(i),x, 0, L_b);
    K_b(i,i) = I*E*simplify(temp);
end

M_p = sym('M_p',[n, n], 'real');
M_p(:,:) = zeros(n, n);
K_p = sym('K_p',[n, n], 'real');
K_p(:,:) = zeros(n, n);
for i = 1:n
    for j = 1:n
        M_p(i,j) = rho_p*A_p*int(Y_k(i).*Y_k(j),x, x_trans(1,1), x_trans(1,2));
        temp = int(d2_Y_k(i).*d2_Y_k(j),x, x_trans(1,1), x_trans(1,2));
        K_p(i,j) = E_p*I_p*simplify(temp);
    end
end


matlabFunction(Xi,'File','Xi_cst_func','Vars',{E,I,rho,A,zeta});
matlabFunction(M_p,'File','M_p_cst_func','Vars',{rho_p,A_p});
matlabFunction(K_p,'File','K_p_cst_func','Vars',{I_p,E_p});

matlabFunction(M_b,'File','M_b_cst_func','Vars',{rho,A});
matlabFunction(K_b,'File','K_b_cst_func','Vars',{I,E});
