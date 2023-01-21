clear all;close all; clc

n = 5; % number of modes
m = 1; % number of transducers


xi = sym('xi', [n,1], 'real');
x = sym('x', 'real');

b_p = sym('b_p',[m,1],  'real');
J_p = sym('J_p',[m,1],  'real');
A_p = sym('A_p',[m,1], 'real');
eps_s = sym('eps_s','real');

L_b = 0.13;
w_b = sym('w_b', 'real');

x_trans = [0.02,0.07];

x_e = L_b/2;

C_n = [3.5160, 22.0345, 61.6972, 120.0902, 199.86];
a = sqrt(C_n)/L_b;
sig_n = [0.734096, 1.018466, 0.999225, 1.000033, 1];

Y_k = (cosh(a*x) - cos(a*x))-(sig_n.*(sinh(a*x) - sin(a*x)));

for i = 1:n
    d_Y_k(i) = gradient(Y_k(i), x);
    d2_Y_k(i) = gradient(d_Y_k(i), x);
end

B_e = (cosh(a*x_e) - cos(a*x_e))-(sig_n.*(sinh(a*x_e) - sin(a*x_e)));


C_a = sym('C_a',[m, m], 'real');
C_a(:,:) = zeros(m, m);
for i = 1:m
    C_a(i,i) = eps_s*w_b^2*(x_trans(i,2) - x_trans(i,1))/A_p(i);
end


B_a = sym('B_a',[n, m], 'real');
B_a(:,:) = zeros(n, m);
for i = 1:n
    for j = 1:m
        temp = int(d2_Y_k(i),x, x_trans(j,1), x_trans(j,2));
        B_a(i,j) = - b_p(j)*J_p(j)/(w_b*(x_trans(j,2) - x_trans(j,1)))*simplify(temp);
    end
end

matlabFunction(B_e,'File','B_e_cst_func','Vars', eps_s);
matlabFunction(C_a,'File','C_a_cst_func','Vars',{eps_s,w_b,A_p});
matlabFunction(B_a,'File','B_a_cst_func','Vars',{b_p,J_p,w_b});


