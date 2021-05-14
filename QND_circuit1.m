% This code reads and plots the results of the circuit measuring concurrence (circuit 1) outputted from QND_Qiskit.py
clear all;

% Basic operations used in the circuit
CNOT=[1 0 0 0 ; 0 1 0 0 ; 0 0 0 1 ; 0 0 1 0];
CNOT_13=[1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0 ; 0 0 0 1 0 0 0 0 ; 0 0 0 0 0 1 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 0 1 ; 0 0 0 0 0 0 1 0];
unity=[1 0 ; 0 1];
Ry_1=[cos(pi/4) -sin(pi/4) ; sin(pi/4) cos(pi/4)];
Ry_2=[cos(-pi/4) -sin(-pi/4) ; sin(-pi/4) cos(-pi/4)];
Rx_1=[cos(pi/4) -i*sin(pi/4) ; -i*sin(pi/4) cos(pi/4)];
Rx_2=[cos(-pi/4) -i*sin(-pi/4) ; -i*sin(-pi/4) cos(-pi/4)];

% Projectors
P0_1q=[1 0 ; 0 0];
P1_1q=[0 0 ; 0 1];

% Projectors on qubits C, D
P0=kron(unity,kron(unity,P0_1q));
P1=kron(unity,kron(unity,P1_1q));

% Unitary operators of the QND circuit
U_C=(kron(Rx_2,kron(Rx_2,unity)))*kron(unity,CNOT)*CNOT_13*(kron(Rx_1,kron(Rx_1,unity)));
U_V=(kron(Ry_2,kron(Ry_2,kron(unity,unity))))*(kron(unity,CNOT_13))*(kron(CNOT_13,unity))*(kron(unity,kron(unity,CNOT)))*(kron(Ry_1,kron(Ry_1,kron(unity,unity))));
U_P=(kron(unity,CNOT_13))*(kron(CNOT_13,unity))*(kron(unity,kron(unity,CNOT)));

% Load initial tomography
tomo_i=load('tomo_initial.txt');

rho_i_exp=read_rho(tomo_i);
phi_i_tomo=tomo_i(:,1);
theta_i_tomo=tomo_i(:,3);
    
QND = load('ND_initial_meas_C.txt');
final1=load('ND_final_meas_C1.txt');
final2=load('ND_final_meas_C2.txt');
final_nops=load('ND_final_meas_Cnops.txt');

C_i_QND=QND(:,4);

for i=1:length(tomo_i(:,1))
   C_i_tomo(i)=concurrence(rho_i_exp(:,:,i)); 
end

phi_f_exp_1=final2(:,1);
theta_f_exp_1=final2(:,3);
phi_f_exp_nops=final_nops(:,1);
theta_f_exp_nops=final_nops(:,3);
rho_f_th_0=zeros(4,4,length(tomo_i(:,1)));
rho_f_th_1=zeros(4,4,length(tomo_i(:,1)));
rho_f_th_nops=zeros(4,4,length(tomo_i(:,1)));
rho_f_th_8=zeros(8,8,length(tomo_i(:,1)));

rho_i_th=zeros(4,4,length(tomo_i(:,1)));
% Compute initial analytical state
for i=1:length(tomo_i(:,1))
    Ry=[cos(phi_i_tomo(i)/2) -sin(phi_i_tomo(i)/2) ; sin(phi_i_tomo(i)/2) cos(phi_i_tomo(i)/2)];
    CRy=[1 0 0 0 ; 0 1 0 0 ; 0 0 cos(theta_i_tomo(i)/2) -sin(theta_i_tomo(i)/2); 0 0 sin(theta_i_tomo(i)/2) cos(theta_i_tomo(i)/2)];
    U_prep=CRy*kron(Ry,unity);
    rho_i_th(:,:,i)=(U_prep*[1;0;0;0])*(U_prep*[1;0;0;0])';
    rho_i_th_8(:,:,i)=kron(rho_i_th(:,:,i),P0_1q);
end

% Compute final analytical state
U=U_C;
for i=1:length(tomo_i(:,1))
    rho_f_th_8(:,:,i)=U*rho_i_th_8(:,:,i)*U';

    % Postselection
    rho_f_th_0(:,:,i)=(1/trace(rho_f_th_8(:,:,i)*P0)).*PartialTrace(P0*rho_f_th_8(:,:,i)*P0,[3],[2,2,2]);
    rho_f_th_1(:,:,i)=(1/trace(rho_f_th_8(:,:,i)*P1)).*PartialTrace(P1*rho_f_th_8(:,:,i)*P1,[3],[2,2,2]);
    rho_f_th_nops(:,:,i)=PartialTrace(rho_f_th_8(:,:,i),[3],[2,2,2]);
end

phi_i_QND=QND(:,1);
theta_i_QND=QND(:,3);
rho_f_exp_1=read_rho(final2);
rho_f_exp_nops=read_rho(final_nops);

for i=1:length(final2(:,1))
   C_f_exp_1(i)=concurrence(rho_f_exp_1(:,:,i));
end

for i=1:length(final_nops(:,1))
   C_f_exp_nops(i)=concurrence(rho_f_exp_nops(:,:,i)); 
end

idx = find(theta_f_exp_1>=pi-0.05 & theta_f_exp_1<=pi+0.05 & phi_f_exp_1>=0 & phi_f_exp_1<=pi);
for i=1:length(idx)
    new_theta_f_exp_1(i)=theta_f_exp_1(idx(i));
    new_phi_f_exp_1(i)=phi_f_exp_1(idx(i));
    new_C_f_exp_1(i)=C_f_exp_1(idx(i));

    idx2=find(theta_i_tomo==new_theta_f_exp_1(i) & phi_i_tomo==new_phi_f_exp_1(i));
    Fidelity_f_1(i)=Fidelity(rho_f_th_1(:,:,idx2),rho_f_exp_1(:,:,idx(i)));
end

idx = find(theta_f_exp_nops>=pi-0.05 & theta_f_exp_nops<=pi+0.05 & phi_f_exp_nops>=0 & phi_f_exp_nops<=pi);
for i=1:length(idx)
    new_theta_f_exp_nops(i)=theta_f_exp_nops(idx(i));
    new_phi_f_exp_nops(i)=phi_f_exp_nops(idx(i));
    new_C_f_exp_nops(i)=C_f_exp_nops(idx(i));

    idx2=find(theta_i_tomo==new_theta_f_exp_nops(i) & phi_i_tomo==new_phi_f_exp_nops(i));
    Fidelity_f_nops(i)=Fidelity(rho_f_th_nops(:,:,idx2),rho_f_exp_nops(:,:,idx(i)));
end

idx = find(theta_i_QND>=pi-0.05 & theta_i_QND<=pi+0.05 & phi_i_QND>=0 & phi_i_QND<=pi);
for i=1:length(idx)
    new_theta_i_QND(i)=theta_i_QND(idx(i));
    new_phi_i_QND(i)=phi_i_QND(idx(i));
    new_C_i_QND(i)=C_i_QND(idx(i));
end

figure
ax = gca;
hold on
plot(new_phi_f_exp_1,Fidelity_f_1,'cx')
hold on
plot(new_phi_f_exp_nops,Fidelity_f_nops,'k+')
xlabel('\phi')
ylabel('Fidelity')
ylim([0 1])
grid on
set(gca, 'FontName', 'LM Roman 12', 'FontSize', 16)
title('Fidelity, measurement C, \theta=\pi')
legend('final state, postselected 1','final state, no postselection')
saveas(gca,'FC_i_1D.png')
saveas(gca,'FC_i_1D.fig')

figure
ax = gca;
plot(new_phi_i_QND,new_C_i_QND,'b.');
hold on
plot(new_phi_f_exp_nops,new_C_f_exp_nops,'k+');
hold on
plot(new_phi_f_exp_1,new_C_f_exp_1,'cx');
xlabel('\phi')
ylabel('C')
ylim([0 1])
grid on
set(gca, 'FontName', 'LM Roman 12', 'FontSize', 16)
title('C, \theta=\pi')
legend('QND','final state, no postselection','final state, postselected 1')
saveas(gca,'C_i_1D.png')
saveas(gca,'C_i_1D.fig')