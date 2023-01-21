clear all; close all; clc
% uniformly distributed, oscillating load for cantilever beam
% load experiment data
data_exp1 = load('free_exp_1.mat');

%% calculate the damping ratio as a funtion of time
[maxs, pos] = findpeaks(data_exp1.disp_exp1,data_exp1.time_exp1, 'MinPeakDistance', 0.1);
maxs(1:4) = [];
pos(1:4) = [];
maxs(end) = [];
pos(end) = [];
start_ind = find(data_exp1.time_exp1 == pos(16));
data_exp1.start_point = start_ind;
for i = 1:length(maxs)-1
    delta = log(maxs(i)/maxs(i+1));
    zeta_arr(i) = delta/(sqrt(delta^2 + (2*pi)^2));
end

figure
hold on
plot(data_exp1.time_exp1,data_exp1.disp_exp1, 'r--', 'LineWidth', 1)
plot(pos,maxs, '*')
ylabel("y [mm]")
yyaxis right
ylabel("\zeta")
plot(pos(1:end-1), zeta_arr, 'k')
grid on
title("Experiment dampening analysis")
legend("Middle point movement","Local maximas", "\zeta Damping ratio")
xlabel("t [Seconds]")
xlim([0,5])

x1 = maxs(19);
x3 = maxs(20);
delta = log(x1/x3);
zeta = delta/(sqrt(delta^2 + (2*pi)^2));

data_exp1.time_exp1 = data_exp1.time_exp1 - data_exp1.time_exp1(data_exp1.start_point);

%% FEM stuff
g = 9.81;

rho = 1.4*10^3;
E = 3.55*10^9;
h_b = 0.5*10^-3;
w_b = 40*10^-3;
A = h_b*w_b;

I = (w_b*h_b^3)/12;
L_b = 0.13;  %meter

nel = 12;
x = linspace(0,L_b, nel+1);

%% Modal shape method
n = 5;
t_span = [0:0.0001:2];
y0 = zeros(n*2,1);
y0(1) = (data_exp1.initial_disp)/2;

[t_mod,y] = ode45('free_modal_func', t_span, y0, zeta);

xi = y(:,1:2:n*2-1);
C_n = [3.5160, 22.0345, 61.6972, 120.0902, 199.86];
a = sqrt(C_n)/L_b;
sig_n = [0.734096, 1.018466, 0.999225, 1.000033, 1];

for i = 1:n
    Y_k(i,:) = (cosh(a(i)*x) - cos(a(i)*x))-(sig_n(i)*(sinh(a(i)*x) - sin(a(i)*x)));
end
u_mod = xi*Y_k;



%% FEM part
end_t = 2;

t = linspace(0,end_t,2000000);
L = x(2);
Ke = K_e_beam_piezo_cst_func(I,E,L);
Me = M_e_beam_piezo_cst_func(rho,A,L);
fe = zeros(4,1);% without *sin(sig_load*t)
K = zeros(2*length(x),2*length(x));
M = zeros(2*length(x),2*length(x));
F = zeros(2*length(x), 1);



for i = 1:length(x)-1
    K(i*2-1:i*2+2,i*2-1:i*2+2) = K(i*2-1:i*2+2,i*2-1:i*2+2)+Ke;
    M(i*2-1:i*2+2,i*2-1:i*2+2) = M(i*2-1:i*2+2,i*2-1:i*2+2)+Me;
    F(i*2-1:i*2+2) = F(i*2-1:i*2+2)+fe;
end
omega = C_n*sqrt(E*I/(A*rho*L_b^4));
beta = (2*zeta*omega(2) - 2*zeta*omega(1))/(omega(2)^2 - omega(1)^2);
alpha = omega(1)*(2*zeta -beta*omega(1));
C = alpha*M + beta*K;


del_t = t(2)-t(1);
nnodes = length(x);
ndofs = nnodes*2;
ntime = length(t);
u_cent = zeros(ndofs, ntime);
dof_F=[1:ndofs];
dof_C = [1,2];
dof_F(dof_C)= [];
node_m = nel/2 + 1;
a_stat = zeros(ndofs,1);
initial_disp = data_exp1.disp_exp1(data_exp1.start_point)*10^-3;
a_stat(end-1) = 0.004184; % corresponds to initial displacement of middle node == maxs(18)
dof_F_stat = dof_F;
dof_F_stat(end-1) = [];
a_stat(dof_F_stat) = inv(K(dof_F_stat,dof_F_stat))*(-K(dof_F_stat,end-1)*a_stat(end-1));
u_cent = zeros(ndofs,ntime);
u_cent(:,1) = a_stat;

A = [zeros(ndofs,ndofs), eye(ndofs, ndofs);
     -inv(M)*K         , -inv(M)*C];

%u = zeros(ndofs, ntime);
u_prime = zeros(ndofs, ntime);
%u(:,1) = a_stat;
f_dofs = [3:ndofs];
nfdofs = length(f_dofs);
f_dofs_comp = [f_dofs, [ndofs+3:ndofs*2]];

c_dofs = [1,2];
a1 = 1/(del_t)^2;
a2 = 1/(2*del_t);
u_cent_0 = a_stat;
for i = 1:ntime-1
    if i == 1
        u_cent(f_dofs,i+1) = (inv(M(dof_F,dof_F)*a1+ C(dof_F,dof_F)*a2))*(F(dof_F) + (M(dof_F,dof_F)*2*a1 - K(dof_F,dof_F))*u_cent(dof_F,i) - (a1*M(dof_F,dof_F) - a2*C(dof_F, dof_F))*u_cent_0(dof_F));
    else
        u_cent(f_dofs,i+1) = (inv(M(dof_F,dof_F)*a1+ C(dof_F,dof_F)*a2))*(F(dof_F) + (M(dof_F,dof_F)*2*a1 - K(dof_F,dof_F))*u_cent(dof_F,i) - (a1*M(dof_F,dof_F) - a2*C(dof_F, dof_F))*u_cent(dof_F, i-1));
    end
%     Xi = vertcat(u(f_dofs,i), u_prime(f_dofs, i));
%     K1 = del_t*(A(f_dofs_comp,f_dofs_comp)*Xi);
%     K2 = del_t*((A(f_dofs_comp,f_dofs_comp)*(Xi + K1/2)));
%     K3 = del_t*((A(f_dofs_comp,f_dofs_comp)*(Xi + K2/2)));
%     K4 = del_t*((A(f_dofs_comp,f_dofs_comp)*(Xi + K3)));
%     Xi_1 = Xi + (1/6)*(K1 + 2*K2 + 2*K3 + K4);
%     u(f_dofs, i+1) = Xi_1(1:nfdofs);
%  %   u(1:2, i+1) = 0;
%     u_prime(f_dofs, i+1) = Xi_1(nfdofs+1:end);
% %    u_prime(1:2, i+1) = 0;
end


%% plots
omega_exp = 1/diff([pos(18), pos(19)]);
% [maxs_FEM, locs_FEM] = findpeaks(u_cent(node_m*2-1, :),t, 'MinPeakDistance', 0.1);
% omega_FEM = 1/diff([locs_FEM(2), locs_FEM(3)]);
% [maxs_FEM, locs_FEM] = findpeaks(u(node_m*2-1, :),t, 'MinPeakDistance', 0.1);
% omega_FEM = 1/diff([locs_FEM(2), locs_FEM(3)]);
figure
%plot(x, y(end,:), 'r', 'LineWidth', 2)
hold on 
%plot(x, u_cent(1:2:end,end-1), '--b', 'LineWidth', 1)

title("Full beam shape")
%legend("Analytical", "FEM")
ylabel("y [Meters]")
xlabel("x [Meters]")

figure
hold on
plot(data_exp1.time_exp1,data_exp1.disp_exp1, 'r', 'LineWidth', 1)
plot(t_mod, u_mod(:,node_m)*10^3,  'k--', 'LineWidth', 2)
plot(t, u_cent(node_m*2-1,:)*10^3,  'b', 'LineWidth', 1)
grid on
title("Middle point vertical displacement",'FontSize', 12)
legend("Experiment", "FEM", "Modal" ,'FontSize', 12)
xlim([0,1.25])
ylim([-1.5,1.5])
ylabel("y [mm]",'FontSize', 12)
xlabel("t [Seconds]",'FontSize', 12)