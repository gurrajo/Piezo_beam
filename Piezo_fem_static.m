function [a_stat, initial_disp] = Piezo_fem_static(disp_m,d)
g = 9.81;
m = 1; % number of transducers
V_tot = 16*60*0.26*10^-9;
V_piezo = 12*50*0.11*10^-9;

mass_p = 0.367*10^-3;
C = 5.8*10^-10;
piezo_L = 50*10^-3;
w_p = 12*10^-3;
h_p = 0.11*10^-3;
eps_S = C*h_p/(w_p*piezo_L);
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
b_p = E_p*d/eps_S;
R = 5*10^6;
zeta = 4*10^-2;

J_p = A_p*(h_p/2 + h_b/2);

I = (w_b*h_b^3)/12;
I_p = ((w_p*h_p^3)/12 + A_p*(h_p/2 + h_b/2)^2);
%% Load case (constant, evenly distributed, transverse load)
F_y = 0; %N/m

omega_load = 0;  %rad/sec



%% Cross section beam
L = 0.13;  %meter


%% FEM related input
nnodes = 10; %number of nodes
nel = nnodes-1; % number of elements
x = [0,0.01,0.02,0.07,0.08,0.09,0.1,0.11,0.12,0.13];
l_beam = 0.01;
l_p = 0.05;
dofs_per_node = 2; % Number of dofs per node
ndofs = nnodes*dofs_per_node; % number of degrees of freedom

%% ANCF related input

M_e_beam = M_e_beam_piezo_cst_func(rho,A,l_beam);
K_e_beam = K_e_beam_piezo_cst_func(I,E,l_beam);

M_vv = Me_piezo_cst_func(rho,A,l_p,rho_p,A_p);
f_ext= f_ext_piezo_cst_func(F_y,l_beam);
K_vv= K_vv_piezo_cst_func(I,E,l_p,E_p,I_p);
K_QQ= K_QQ_piezo_cst_func(A_p,w_p,eps_S,l_p);
K_Qv = K_Qv_piezo_cst_func(b_p,w_p,J_p,l_p);
K_vQ = K_vQ_piezo_cst_func(b_p,w_p,J_p,l_p);

F = zeros(ndofs,1);
K = zeros(ndofs,ndofs); 
M = zeros(ndofs,ndofs);

dof_F=[1:ndofs]; % list of free degrees of freedom

trans_nodes = [3,4];
trans_el = 1;
dof_C = [1,2];

dof_F(dof_C) = [];  %removing the prescribed dofs from dof_F

% assemble non time dependent stuff

trans_v_dofs = [trans_nodes(1)*2-1:trans_nodes(end)*2];


for el = 1:nel
    el_nodes = [el, el+1];
    % cordinate ordering, (x,dx,x,dx....Q,Q,Q)
    el_dofs = [el*2-1:el*2+2];
    
    if sum(ismember(el_nodes,trans_nodes)) == 2
        % if a transducer node is incuded in the element
        M_e = M_vv;
        K_e = K_vv;
        f_e = f_ext;
        
    else
        % If regular beam element
        M_e = M_e_beam;
        K_e = K_e_beam;
        f_e = f_ext;
    end
    K(el_dofs,el_dofs) = K(el_dofs,el_dofs)+K_e;
    M(el_dofs,el_dofs) = M(el_dofs,el_dofs)+M_e;

    F(el_dofs) = f_e;
end
% initialize for time dependent assembly
ntime = 500000; % time steps
tend = 3; % end time (sec)
t = linspace(0,tend,ntime+1); % time vector

a = zeros(ndofs,ntime); % displacement field, function of (degrees of freedom, time)
Q = zeros(1,ntime+1);
%% Initial displacement using beam normal function
C_n = [3.5160, 22.0345];
beta = sqrt(C_n(1))/L;
sig_n = 0.734096;

Y_k = (cosh(beta*x) - cos(beta*x))-(sig_n*(sinh(beta*x) - sin(beta*x)));

% remaining dofs have initial value 0, (y, dr2/dx)

%% Central differencing variables
del_t = t(2);
a1 = 1/del_t^2;
a2 = 1/(2*del_t);

C_Q = zeros(trans_el,trans_el);
K_QQ_comp = zeros(trans_el,trans_el);
K_Qv_comp = zeros(trans_el,(trans_el+1)*2);
for i = 1:trans_el
    el_nodes = [i*2-1:i*2+2];
    K_QQ_comp(i,i) = K_QQ;
    K_Qv_comp(i,el_nodes) = K_Qv_comp(i,el_nodes) + K_Qv;
    C_Q(i,i) = R/(del_t*l_beam);
end

%% Intial conditions
% vertical displacemnt accoring to first normal mode
eps = 1*10^-6;

a_stat = zeros(ndofs,1);
dof_F_stat = dof_F;
dof_F_stat(end-1) = [];
F_stat = zeros(ndofs,1);
node_m = 4;
a_stat(end-1) = disp_m*2;
a_stat(dof_F_stat) = inv(K(dof_F_stat,dof_F_stat))*(-K(dof_F_stat,end-1)*a_stat(end-1));
while abs(a_stat(node_m*2-1) - disp_m) > eps
    F_stat(trans_v_dofs) = - K_Qv_comp'*inv(K_QQ_comp)*(K_Qv_comp*a_stat(trans_v_dofs));
    a_stat_new = inv(K(dof_F_stat,dof_F_stat))*(F_stat(dof_F_stat)-K(dof_F_stat,end-1)*a_stat(end-1)); 
    a_stat(dof_F_stat) = 0.9*a_stat(dof_F_stat) + 0.1*a_stat_new;
    
    if abs(a_stat(node_m*2-1) - disp_m) > eps
        a_stat(end-1) = a_stat(end-1) + (disp_m - a_stat(node_m*2-1));
    end 
    
end

initial_disp = a_stat(end-1);
end