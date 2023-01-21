function [t,u,Q] = Piezo_fem_imp(a_stat,d)
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
a(:,1) = a_stat;
a(:,2) = a_stat;

omega = C_n*sqrt(E*I/(A*rho*L^4));
beta = (2*zeta*omega(2) - 2*zeta*omega(1))/(omega(2)^2 - omega(1)^2);
alpha = omega(1)*(2*zeta -beta*omega(1));
C = alpha*M + beta*K;
%% Setup and solve FE equations
% Assemble time dependent stiffness matrix 

for i = 1:ntime
    %converge = 0;
    %while converge == 0
    %Q_i = sum(Q(:,i));
    F = F*.0;
%    for el = 1: nel
%        el_nodes = [el,el+1];
%        el_dofs = [el*2-1:el*2+2];
%             if sum(ismember(el_nodes,trans_nodes)) == 2
%                 % if transducer
%                 m = (l_p*A_p*rho_p + l_p*A*rho);
%                 F(el_dofs) = F(el_dofs) + [0, -a(el_dofs(3),i)*g*m/2,0,a(el_dofs(3),i)*g*m/2]';
%             else
%                 m = l_beam*A*rho;
%                 F(el_dofs) = F(el_dofs) + [0, -a(el_dofs(3),i)*g*m/2,0,a(el_dofs(3),i)*g*m/2]';
%             end
%
%    end
    Q(i+1) = inv(K_QQ_comp)*(K_Qv_comp*a(trans_v_dofs,i));
    if i > 1
        F(trans_v_dofs) = F(trans_v_dofs) - K_Qv_comp'*Q(i+1);
        a(dof_F,i+1) = inv(M(dof_F,dof_F)*a1 + C(dof_F,dof_F)*a2)*(F(dof_F)  + (M(dof_F,dof_F)*2*a1 - K(dof_F,dof_F))*a(dof_F,i) - (a1*M(dof_F,dof_F) - a2*C(dof_F,dof_F))*a(dof_F, i-1));
    end
end

u = a(1:2:end,:);

end