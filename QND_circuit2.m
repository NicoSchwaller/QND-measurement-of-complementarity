% This code reads and plots the results of the general circuit (circuit 2) outputted from QND_Qiskit.py
clear all;

% Choose the quantity to read
VPC=3;

% Basic operations used in the circuit
CNOT=[1 0 0 0 ; 0 1 0 0 ; 0 0 0 1 ; 0 0 1 0];
CNOT_13=[1 0 0 0 0 0 0 0 ; 0 1 0 0 0 0 0 0 ; 0 0 1 0 0 0 0 0 ; 0 0 0 1 0 0 0 0 ; 0 0 0 0 0 1 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 0 1 ; 0 0 0 0 0 0 1 0];
unity=[1 0 ; 0 1];
Ry_1=[cos(pi/4) -sin(pi/4) ; sin(pi/4) cos(pi/4)];
Ry_2=[cos(-pi/4) -sin(-pi/4) ; sin(-pi/4) cos(-pi/4)];
Rx_1=[cos(pi/4) -i*sin(pi/4) ; -i*sin(pi/4) cos(pi/4)];
Rx_2=[cos(-pi/4) -i*sin(-pi/4) ; -i*sin(-pi/4) cos(-pi/4)];

% Projectors
P0=[1 0 ; 0 0];
P1=[0 0 ; 0 1];

% Projectors on qubits C, D
P00=kron(unity,kron(unity,kron(P0,P0)));
P01=kron(unity,kron(unity,kron(P0,P1)));
P10=kron(unity,kron(unity,kron(P1,P0)));
P11=kron(unity,kron(unity,kron(P1,P1)));

% Unitary operators of the QND circuit
U_C=(kron(Rx_2,kron(Rx_2,kron(unity,unity))))*(kron(unity,CNOT_13))*(kron(CNOT_13,unity))*(kron(unity,kron(unity,CNOT)))*(kron(Rx_1,kron(Rx_1,kron(Ry_1,unity))));
U_V=(kron(Ry_2,kron(Ry_2,kron(unity,unity))))*(kron(unity,CNOT_13))*(kron(CNOT_13,unity))*(kron(unity,kron(unity,CNOT)))*(kron(Ry_1,kron(Ry_1,kron(unity,unity))));
U_P=(kron(unity,CNOT_13))*(kron(CNOT_13,unity))*(kron(unity,kron(unity,CNOT)));

% Load initial tomography
tomo_i=load('tomo_initial.txt');
rho_i_th=zeros(4,4,length(tomo_i(:,1)));
rho_i_exp=read_rho(tomo_i);
phi_i_tomo=tomo_i(:,1);
theta_i_tomo=tomo_i(:,3);

% Compute initial analytical state
for i=1:length(tomo_i(:,1))
    Ry=[cos(phi_i_tomo(i)/2) -sin(phi_i_tomo(i)/2) ; sin(phi_i_tomo(i)/2) cos(phi_i_tomo(i)/2)];
    CRy=[1 0 0 0 ; 0 1 0 0 ; 0 0 cos(theta_i_tomo(i)/2) -sin(theta_i_tomo(i)/2); 0 0 sin(theta_i_tomo(i)/2) cos(theta_i_tomo(i)/2)];
    U_prep=CRy*kron(Ry,unity);
    rho_i_th(:,:,i)=(U_prep*[1;0;0;0])*(U_prep*[1;0;0;0])';
    rho_i_th_16(:,:,i)=kron(rho_i_th(:,:,i),kron(P0, P0));
end

% Initialize density matrices
rho_f_th_00=zeros(4,4,length(tomo_i(:,1)));
rho_f_th_01=zeros(4,4,length(tomo_i(:,1)));
rho_f_th_10=zeros(4,4,length(tomo_i(:,1)));
rho_f_th_11=zeros(4,4,length(tomo_i(:,1)));
rho_f_th_nops=zeros(4,4,length(tomo_i(:,1)));
rho_f_th_16=zeros(16,16,length(tomo_i(:,1)));

%% Concurrence

if VPC==3
    % Compute final analytical state
    U=U_C;
    for i=1:length(tomo_i(:,1))

        rho_f_th_16(:,:,i)=U*rho_i_th_16(:,:,i)*U';

        % Postselection
        rho_f_th_00(:,:,i)=(1/trace(rho_f_th_16(:,:,i)*P00)).*PartialTrace(P00*rho_f_th_16(:,:,i)*P00,[3,4],[2,2,2,2]);
        rho_f_th_01(:,:,i)=(1/trace(rho_f_th_16(:,:,i)*P01)).*PartialTrace(P01*rho_f_th_16(:,:,i)*P01,[3,4],[2,2,2,2]);
        rho_f_th_10(:,:,i)=(1/trace(rho_f_th_16(:,:,i)*P10)).*PartialTrace(P10*rho_f_th_16(:,:,i)*P10,[3,4],[2,2,2,2]);
        rho_f_th_11(:,:,i)=(1/trace(rho_f_th_16(:,:,i)*P11)).*PartialTrace(P11*rho_f_th_16(:,:,i)*P11,[3,4],[2,2,2,2]);
        rho_f_th_nops(:,:,i)=PartialTrace(rho_f_th_16(:,:,i),[3,4],[2,2,2,2]);
    end

    QND = load('ND_initial_meas_C.txt');
    final1=load('ND_final_meas_C1.txt');
    final2=load('ND_final_meas_C2.txt');
    final3=load('ND_final_meas_C3.txt');
    final4=load('ND_final_meas_C4.txt');
    final_nops=load('ND_final_meas_Cnops.txt');
    
    C_i_QND=QND(:,4);
    depth=QND(:,5);
    u1=QND(:,6);
    u2=QND(:,7);
    u3=QND(:,8);
    cx=QND(:,9);
    depth_raw=QND(:,10);
    u1_raw=QND(:,11);
    u2_raw=QND(:,12);
    u3_raw=QND(:,13);
    cx_raw=QND(:,14);
            
    for i=1:length(tomo_i(:,1))
       C_i_tomo(i)=concurrence(rho_i_exp(:,:,i)); 
    end

    phi_f_exp_00=final1(:,1);
    theta_f_exp_00=final1(:,3);
    phi_f_exp_01=final2(:,1);
    theta_f_exp_01=final2(:,3);
    phi_f_exp_10=final3(:,1);
    theta_f_exp_10=final3(:,3);
    phi_f_exp_11=final4(:,1);
    theta_f_exp_11=final4(:,3);
    phi_f_exp_nops=final_nops(:,1);
    theta_f_exp_nops=final_nops(:,3);
    phi_i_QND=QND(:,1);
    theta_i_QND=QND(:,3);
    rho_f_exp_00=read_rho(final1);
    rho_f_exp_01=read_rho(final2);
    rho_f_exp_10=read_rho(final3);
    rho_f_exp_11=read_rho(final4);
    rho_f_exp_nops=read_rho(final_nops);

    for i=1:length(final1(:,1))
       C_f_exp_00(i)=concurrence(rho_f_exp_00(:,:,i));
    end
    for i=1:length(final2(:,1))
       C_f_exp_01(i)=concurrence(rho_f_exp_01(:,:,i));
    end
    
    for i=1:length(final3(:,1))
       C_f_exp_10(i)=concurrence(rho_f_exp_10(:,:,i));
    end
    
    for i=1:length(final4(:,1))
       C_f_exp_11(i)=concurrence(rho_f_exp_11(:,:,i));
    end
       
    for i=1:length(final_nops(:,1))
       C_f_exp_nops(i)=concurrence(rho_f_exp_nops(:,:,i)); 
    end
    
    idx = find(theta_f_exp_00>=pi-0.05 & theta_f_exp_00<=pi+0.05 & phi_f_exp_00>=0 & phi_f_exp_00<=pi);
    for i=1:length(idx)
        new_theta_f_exp_00(i)=theta_f_exp_00(idx(i));
        new_phi_f_exp_00(i)=phi_f_exp_00(idx(i));
        new_C_f_exp_00(i)=C_f_exp_00(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_00(i) & phi_i_tomo==new_phi_f_exp_00(i));
        Fidelity_f_00(i)=Fidelity(rho_f_th_00(:,:,idx2),rho_f_exp_00(:,:,idx(i)));
    end

    idx = find(theta_f_exp_01>=pi-0.05 & theta_f_exp_01<=pi+0.05 & phi_f_exp_01>=0 & phi_f_exp_01<=pi);
    for i=1:length(idx)
        new_theta_f_exp_01(i)=theta_f_exp_01(idx(i));
        new_phi_f_exp_01(i)=phi_f_exp_01(idx(i));
        new_C_f_exp_01(i)=C_f_exp_01(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_01(i) & phi_i_tomo==new_phi_f_exp_01(i));
        Fidelity_f_01(i)=Fidelity(rho_f_th_01(:,:,idx2),rho_f_exp_01(:,:,idx(i)));
    end
    
    idx = find(theta_f_exp_10>=pi-0.05 & theta_f_exp_10<=pi+0.05 & phi_f_exp_10>=0 & phi_f_exp_10<=pi);
    for i=1:length(idx)
        new_theta_f_exp_10(i)=theta_f_exp_10(idx(i));
        new_phi_f_exp_10(i)=phi_f_exp_10(idx(i));
        new_C_f_exp_10(i)=C_f_exp_10(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_10(i) & phi_i_tomo==new_phi_f_exp_10(i));
        Fidelity_f_10(i)=Fidelity(rho_f_th_10(:,:,idx2),rho_f_exp_10(:,:,idx(i)));
    end
    
    idx = find(theta_f_exp_11>=pi-0.05 & theta_f_exp_11<=pi+0.05 & phi_f_exp_11>=0 & phi_f_exp_11<=pi);
    for i=1:length(idx)
        new_theta_f_exp_11(i)=theta_f_exp_11(idx(i));
        new_phi_f_exp_11(i)=phi_f_exp_11(idx(i));
        new_C_f_exp_11(i)=C_f_exp_11(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_11(i) & phi_i_tomo==new_phi_f_exp_11(i));
        Fidelity_f_11(i)=Fidelity(rho_f_th_11(:,:,idx2),rho_f_exp_11(:,:,idx(i)));
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
        new_u1(i)=u1(idx(i));
        new_u2(i)=u2(idx(i));
        new_u3(i)=u3(idx(i));
        new_cx(i)=cx(idx(i));
        new_depth(i)=depth(idx(i));
        new_u1_raw(i)=u1_raw(idx(i));
        new_u2_raw(i)=u2_raw(idx(i));
        new_u3_raw(i)=u3_raw(idx(i));
        new_cx_raw(i)=cx_raw(idx(i));
        new_depth_raw(i)=depth_raw(idx(i));
    end
    
    idx = find(theta_i_tomo>=pi-0.05 & theta_i_tomo<=pi+0.05 & phi_i_tomo>=0 & phi_i_tomo<=pi);
    for i=1:length(idx)
        new_theta_i_tomo(i)=theta_i_tomo(idx(i));
        new_phi_i_tomo(i)=phi_i_tomo(idx(i));
        new_C_i_tomo(i)=C_i_tomo(idx(i));
        Fidelity_i_tomo(i)=Fidelity(rho_i_th(:,:,idx(i)),rho_i_exp(:,:,idx(i)));
    end

    figure
    ax = gca;
    plot(new_phi_f_exp_00,Fidelity_f_00,'rx')
    hold on
    plot(new_phi_f_exp_01,Fidelity_f_01,'cx')
    hold on
    plot(new_phi_f_exp_10,Fidelity_f_10,'c^')
    hold on
    plot(new_phi_f_exp_11,Fidelity_f_11,'ro')
    hold on
    plot(new_phi_f_exp_nops,Fidelity_f_nops,'k+')
    hold on
    plot(new_phi_i_tomo,Fidelity_i_tomo,'bo')
    xlabel('\phi')
    ylabel('Fidelity')
    ylim([0 1]); grid on;
    set(gca, 'FontName', 'LM Roman 12', 'FontSize', 16)
    title('Fidelity, measurement C, \theta=\pi')
    legend('final state, 00 postselected','final state, 01 postselected','final state, 10 postselected','final state, 11 postselected','final state, no postselection','initial state')
    saveas(gca,'FC_i_1D.png')
    saveas(gca,'FC_i_1D.fig')

    figure
    ax = gca;
    plot(new_phi_i_tomo,new_C_i_tomo,'bo');
    hold on
    plot(new_phi_i_QND,new_C_i_QND,'b.');
    hold on
    plot(new_phi_f_exp_nops,new_C_f_exp_nops,'k+');
    hold on
    plot(new_phi_f_exp_00,new_C_f_exp_00,'rx');
    hold on
    plot(new_phi_f_exp_01,new_C_f_exp_01,'cx');
    hold on
    plot(new_phi_f_exp_10,new_C_f_exp_10,'c^');
    hold on
    plot(new_phi_f_exp_11,new_C_f_exp_11,'ro');
    xlabel('\phi')
    ylabel('C')
    ylim([0 1]); grid on;
    set(gca, 'FontName', 'LM Roman 12', 'FontSize', 16)
    title('C, \theta=\pi')
    legend('initial tomography','QND','final state, no postselection','final state, 00 postselected','final state, 01 postselected','final state, 10 postselected','final state, 11 postselected')
    saveas(gca,'C_i_1D.png')
    saveas(gca,'C_i_1D.fig')
end

%% Visibility

if VPC==1
    
    % Compute final analytical state
    U=U_V;
    for i=1:length(tomo_i(:,1))

            rho_f_th_16(:,:,i)=U*rho_i_th_16(:,:,i)*U';

            % Postselection
            rho_f_th_00(:,:,i)=(1/trace(rho_f_th_16(:,:,i)*P00)).*PartialTrace(P00*rho_f_th_16(:,:,i)*P00,[3,4],[2,2,2,2]);
            rho_f_th_01(:,:,i)=(1/trace(rho_f_th_16(:,:,i)*P01)).*PartialTrace(P01*rho_f_th_16(:,:,i)*P01,[3,4],[2,2,2,2]);
            rho_f_th_10(:,:,i)=(1/trace(rho_f_th_16(:,:,i)*P10)).*PartialTrace(P10*rho_f_th_16(:,:,i)*P10,[3,4],[2,2,2,2]);
            rho_f_th_11(:,:,i)=(1/trace(rho_f_th_16(:,:,i)*P11)).*PartialTrace(P11*rho_f_th_16(:,:,i)*P11,[3,4],[2,2,2,2]);
            rho_f_th_nops(:,:,i)=PartialTrace(rho_f_th_16(:,:,i),[3,4],[2,2,2,2]);
    end
    
    QND = load('ND_initial_meas_V.txt');
    final1=load('ND_final_meas_V1.txt');
    final2=load('ND_final_meas_V2.txt');
    final3=load('ND_final_meas_V3.txt');
    final4=load('ND_final_meas_V4.txt');
    final_nops=load('ND_final_meas_Vnops.txt');
    
    rho_f_exp_00=read_rho(final1);
    rho_f_exp_01=read_rho(final2);
    rho_f_exp_10=read_rho(final3);
    rho_f_exp_11=read_rho(final4);
    rho_f_exp_nops=read_rho(final_nops);
    
    VA_i_QND=QND(:,4);
    VB_i_QND=QND(:,5);
    
    for i=1:length(tomo_i(:,1))
       VA_i_tomo(i)=visibility('A',rho_i_exp(:,:,i));
       VB_i_tomo(i)=visibility('B',rho_i_exp(:,:,i));
    end

    phi_f_exp_00=final1(:,1);
    theta_f_exp_00=final1(:,3);
    try
        phi_f_exp_01=final2(:,1);
        theta_f_exp_01=final2(:,3);
    catch
        phi_f_exp_01=[];
        theta_f_exp_01=[];
    end
    try
        phi_f_exp_10=final3(:,1);
        theta_f_exp_10=final3(:,3);
    catch
        phi_f_exp_10=[];
        theta_f_exp_10=[];
    end
        try
        phi_f_exp_11=final4(:,1);
        theta_f_exp_11=final4(:,3);
    catch
        phi_f_exp_11=[];
        theta_f_exp_11=[];
    end
    phi_f_exp_nops=final_nops(:,1);
    theta_f_exp_nops=final_nops(:,3);
    phi_i_tomo=tomo_i(:,1);
    theta_i_tomo=tomo_i(:,3);
    phi_i_QND=QND(:,1);
    theta_i_QND=QND(:,3);

    for i=1:length(final1(:,1))
       VA_f_exp_00(i)=visibility('A',rho_f_exp_00(:,:,i));
       VB_f_exp_00(i)=visibility('B',rho_f_exp_00(:,:,i));
    end
    try
    for i=1:length(final2(:,1))
       VA_f_exp_01(i)=visibility('A',rho_f_exp_01(:,:,i));
       VB_f_exp_01(i)=visibility('B',rho_f_exp_01(:,:,i));
    end
    catch
       VA_f_exp_01=[];
       VB_f_exp_01=[];
    end
    try
    for i=1:length(final3(:,1))
       VA_f_exp_10(i)=visibility('A',rho_f_exp_10(:,:,i));
       VB_f_exp_10(i)=visibility('B',rho_f_exp_10(:,:,i));
    end
    catch
       VA_f_exp_10=[];
       VB_f_exp_10=[];
    end
    
    try
    for i=1:length(final4(:,1))
       VA_f_exp_11(i)=visibility('A',rho_f_exp_11(:,:,i));
       VB_f_exp_11(i)=visibility('B',rho_f_exp_11(:,:,i));
    end
    catch
       VA_f_exp_11=[];
       VB_f_exp_11=[];
    end
       
    for i=1:length(final_nops(:,1))
       VA_f_exp_nops(i)=visibility('A',rho_f_exp_nops(:,:,i)); 
       VB_f_exp_nops(i)=visibility('B',rho_f_exp_nops(:,:,i));
    end
    
    val_theta=0;
    idx = find(theta_f_exp_00>=val_theta-0.05 & theta_f_exp_00<=val_theta+0.05 & phi_f_exp_00>=0 & phi_f_exp_00<=pi);
    for i=1:length(idx)
        new_theta_f_exp_00(i)=theta_f_exp_00(idx(i));
        new_phi_f_exp_00(i)=phi_f_exp_00(idx(i));
        new_VA_f_exp_00(i)=VA_f_exp_00(idx(i));
        new_VB_f_exp_00(i)=VB_f_exp_00(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_00(i) & phi_i_tomo==new_phi_f_exp_00(i));
        try
            Fidelity_f_00(i)=Fidelity(rho_f_th_00(:,:,idx2),rho_f_exp_00(:,:,idx(i)));
        catch
            Fidelity_f_01(i)=NaN;
        end
    end
    
    idx = find(theta_f_exp_01>=val_theta-0.05 & theta_f_exp_01<=val_theta+0.05 & phi_f_exp_01>=0 & phi_f_exp_01<=pi);
    for i=1:length(idx)
        new_theta_f_exp_01(i)=theta_f_exp_01(idx(i));
        new_phi_f_exp_01(i)=phi_f_exp_01(idx(i));
        new_VA_f_exp_01(i)=VA_f_exp_01(idx(i));
        new_VB_f_exp_01(i)=VB_f_exp_01(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_01(i) & phi_i_tomo==new_phi_f_exp_01(i));
        try
            Fidelity_f_01(i)=Fidelity(rho_f_th_01(:,:,idx2),rho_f_exp_01(:,:,idx(i)));
        catch
            Fidelity_f_01(i)=NaN;
        end
    end
    
    idx = find(theta_f_exp_10>=val_theta-0.05 & theta_f_exp_10<=val_theta+0.05 & phi_f_exp_10>=0 & phi_f_exp_10<=pi);
    for i=1:length(idx)
        new_theta_f_exp_10(i)=theta_f_exp_10(idx(i));
        new_phi_f_exp_10(i)=phi_f_exp_10(idx(i));
        new_VA_f_exp_10(i)=VA_f_exp_10(idx(i));
        new_VB_f_exp_10(i)=VB_f_exp_10(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_10(i) & phi_i_tomo==new_phi_f_exp_10(i));
        try
            Fidelity_f_10(i)=Fidelity(rho_f_th_10(:,:,idx2),rho_f_exp_10(:,:,idx(i)));
        catch
            Fidelity_f_10(i)=NaN;
        end
    end
    
    idx = find(theta_f_exp_11>=val_theta-0.05 & theta_f_exp_11<=val_theta+0.05 & phi_f_exp_11>=0 & phi_f_exp_11<=pi);
    for i=1:length(idx)
        new_theta_f_exp_11(i)=theta_f_exp_11(idx(i));
        new_phi_f_exp_11(i)=phi_f_exp_11(idx(i));
        new_VA_f_exp_11(i)=VA_f_exp_11(idx(i));
        new_VB_f_exp_11(i)=VB_f_exp_11(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_11(i) & phi_i_tomo==new_phi_f_exp_11(i));
        try
            Fidelity_f_11(i)=Fidelity(rho_f_th_11(:,:,idx2),rho_f_exp_11(:,:,idx(i)));
        catch
            Fidelity_f_11(i)=NaN;
        end
    end
    
    idx = find(theta_f_exp_nops>=val_theta-0.05 & theta_f_exp_nops<=val_theta+0.05 & phi_f_exp_nops>=0 & phi_f_exp_nops<=pi);
    for i=1:length(idx)
        new_theta_f_exp_nops(i)=theta_f_exp_nops(idx(i));
        new_phi_f_exp_nops(i)=phi_f_exp_nops(idx(i));
        new_VA_f_exp_nops(i)=VA_f_exp_nops(idx(i));
        new_VB_f_exp_nops(i)=VB_f_exp_nops(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_nops(i) & phi_i_tomo==new_phi_f_exp_nops(i));
        Fidelity_f_nops(i)=Fidelity(rho_f_th_nops(:,:,idx2),rho_f_exp_nops(:,:,idx(i)));
    end
    
    idx = find(theta_i_QND>=val_theta-0.05 & theta_i_QND<=val_theta+0.05 & phi_i_QND>=0 & phi_i_QND<=pi);
    for i=1:length(idx)
        new_theta_i_QND(i)=theta_i_QND(idx(i));
        new_phi_i_QND(i)=phi_i_QND(idx(i));
        new_VA_i_QND(i)=VA_i_QND(idx(i));
        new_VB_i_QND(i)=VB_i_QND(idx(i));
    end
    
    idx = find(theta_i_tomo>=val_theta-0.05 & theta_i_tomo<=val_theta+0.05 & phi_i_tomo>=0 & phi_i_tomo<=pi);
    for i=1:length(idx)
        new_theta_i_tomo(i)=theta_i_tomo(idx(i));
        new_phi_i_tomo(i)=phi_i_tomo(idx(i));
        new_VA_i_tomo(i)=VA_i_tomo(idx(i));
        new_VB_i_tomo(i)=VB_i_tomo(idx(i));

        Fidelity_i_tomo(i)=Fidelity(rho_i_th(:,:,idx(i)),rho_i_exp(:,:,idx(i)));
    end

    figure
    ax = gca;
    plot(new_phi_f_exp_00,Fidelity_f_00,'rx')
    hold on
    plot(new_phi_f_exp_01,Fidelity_f_01,'cx')
    hold on
    plot(new_phi_f_exp_10,Fidelity_f_10,'c^')
    hold on
    plot(new_phi_f_exp_11,Fidelity_f_11,'ro')
    hold on
    plot(new_phi_f_exp_nops,Fidelity_f_nops,'k+')
    hold on
    plot(new_phi_i_tomo,Fidelity_i_tomo,'bo')
    xlabel('\phi')
    ylabel('Fidelity')
    ylim([0 1]); grid on;
    set(gca, 'FontName', 'LM Roman 12', 'FontSize', 16)
    title('Fidelity, measurement V, \theta=\0')
    legend('final state, 00 postselected','final state, 01 postselected','final state, 10 postselected','final state, 11 postselected','final state, no postselection','initial state')
    saveas(gca,'FVA_i_1D.png')
    saveas(gca,'FVA_i_1D.fig')
 
    Fidelity_f_nops=[];
    Fidelity_f_00=[];
    Fidelity_f_01=[];
    Fidelity_f_10=[];
    Fidelity_f_11=[];
    Fidelity_i_tomo=[];

    figure
    ax = gca;
    plot(new_phi_i_tomo,new_VA_i_tomo,'bo');
    hold on
    plot(new_phi_i_QND,new_VA_i_QND,'b.');
    hold on
    plot(new_phi_f_exp_nops,new_VA_f_exp_nops,'k+');
    hold on
    plot(new_phi_f_exp_00,new_VA_f_exp_00,'rx');
    hold on
    plot(new_phi_f_exp_01,new_VA_f_exp_01,'gv');
    hold on
    plot(new_phi_f_exp_10,new_VA_f_exp_10,'g^');
    hold on
    plot(new_phi_f_exp_11,new_VA_f_exp_11,'ro');
    xlabel('\phi')
    ylabel('V_A')
    ylim([0 1]); grid on;
    set(gca, 'FontName', 'LM Roman 12', 'FontSize', 16)
    title('V_A, \theta=0')
    legend('initial tomography','QND','final state, no postselection','final state, 00 postselected','final state, 01 postselected','final state, 10 postselected','final state, 11 postselected')
    saveas(gca,'VA_i_1D.png')
    saveas(gca,'VA_i_1D.fig')

    clear new_theta_f_exp_00
    clear new_theta_f_exp_01
    clear new_theta_f_exp_10
    clear new_theta_f_exp_11
    clear new_theta_f_exp_nops
    clear new_phi_f_exp_00
    clear new_phi_f_exp_01
    clear new_phi_f_exp_10
    clear new_phi_f_exp_11
    clear new_phi_f_exp_nops
    clear new_VA_f_exp_00
    clear new_VB_f_exp_00
    clear new_VA_f_exp_01
    clear new_VB_f_exp_01
    clear new_VA_f_exp_10
    clear new_VB_f_exp_10
    clear new_VA_f_exp_11
    clear new_VB_f_exp_11
    clear new_theta_i_QND
    clear new_phi_i_QND
    clear new_VA_i_QND
    clear new_VB_i_QND
    clear new_theta_i_tomo
    clear new_phi_i_tomo
    clear new_VA_i_tomo
    clear new_VB_i_tomo

    val_theta=3*pi/2;
    idx = find(theta_f_exp_00>=val_theta-0.05 & theta_f_exp_00<=val_theta+0.05 & phi_f_exp_00>=0 & phi_f_exp_00<=2*pi);
    for i=1:length(idx)
        new_theta_f_exp_00(i)=theta_f_exp_00(idx(i));
        new_phi_f_exp_00(i)=phi_f_exp_00(idx(i));
        new_VA_f_exp_00(i)=VA_f_exp_00(idx(i));
        new_VB_f_exp_00(i)=VB_f_exp_00(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_00(i) & phi_i_tomo==new_phi_f_exp_00(i));
        try
            Fidelity_f_00(i)=Fidelity(rho_f_th_00(:,:,idx2),rho_f_exp_00(:,:,idx(i)));
        catch
            Fidelity_f_00(i)=NaN;
        end
    end
    
    idx = find(theta_f_exp_01>=val_theta-0.05 & theta_f_exp_01<=val_theta+0.05 & phi_f_exp_01>=0 & phi_f_exp_01<=2*pi);
    for i=1:length(idx)
        new_theta_f_exp_01(i)=theta_f_exp_01(idx(i));
        new_phi_f_exp_01(i)=phi_f_exp_01(idx(i));
        new_VA_f_exp_01(i)=VA_f_exp_01(idx(i));
        new_VB_f_exp_01(i)=VB_f_exp_01(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_01(i) & phi_i_tomo==new_phi_f_exp_01(i));
        try
            Fidelity_f_01(i)=Fidelity(rho_f_th_01(:,:,idx2),rho_f_exp_01(:,:,idx(i)));
        catch
            Fidelity_f_01(i)=NaN;
        end
    end
    
    idx = find(theta_f_exp_10>=val_theta-0.05 & theta_f_exp_10<=val_theta+0.05 & phi_f_exp_10>=0 & phi_f_exp_10<=2*pi);
    for i=1:length(idx)
        new_theta_f_exp_10(i)=theta_f_exp_10(idx(i));
        new_phi_f_exp_10(i)=phi_f_exp_10(idx(i));
        new_VA_f_exp_10(i)=VA_f_exp_10(idx(i));
        new_VB_f_exp_10(i)=VB_f_exp_10(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_10(i) & phi_i_tomo==new_phi_f_exp_10(i));
        try
            Fidelity_f_10(i)=Fidelity(rho_f_th_10(:,:,idx2),rho_f_exp_10(:,:,idx(i)));
        catch
            Fidelity_f_10(i)=NaN;
        end
    end
    
    idx = find(theta_f_exp_11>=val_theta-0.05 & theta_f_exp_11<=val_theta+0.05 & phi_f_exp_11>=0 & phi_f_exp_11<=2*pi);
    for i=1:length(idx)
        new_theta_f_exp_11(i)=theta_f_exp_11(idx(i));
        new_phi_f_exp_11(i)=phi_f_exp_11(idx(i));
        new_VA_f_exp_11(i)=VA_f_exp_11(idx(i));
        new_VB_f_exp_11(i)=VB_f_exp_11(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_11(i) & phi_i_tomo==new_phi_f_exp_11(i));
        try
            Fidelity_f_11(i)=Fidelity(rho_f_th_11(:,:,idx2),rho_f_exp_11(:,:,idx(i)));
        catch
            Fidelity_f_11(i)=NaN;
        end
    end
    
    idx = find(theta_f_exp_nops>=val_theta-0.05 & theta_f_exp_nops<=val_theta+0.05 & phi_f_exp_nops>=0 & phi_f_exp_nops<=2*pi);
    for i=1:length(idx)
        new_theta_f_exp_nops(i)=theta_f_exp_nops(idx(i));
        new_phi_f_exp_nops(i)=phi_f_exp_nops(idx(i));
        new_VA_f_exp_nops(i)=VA_f_exp_nops(idx(i));
        new_VB_f_exp_nops(i)=VB_f_exp_nops(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_nops(i) & phi_i_tomo==new_phi_f_exp_nops(i));
        try
            Fidelity_f_nops(i)=Fidelity(rho_f_th_nops(:,:,idx2),rho_f_exp_nops(:,:,idx(i)));
        catch
            Fidelity_f_nops(i)=NaN;
        end
    end
    
    idx = find(theta_i_QND>=val_theta-0.05 & theta_i_QND<=val_theta+0.05 & phi_i_QND>=0 & phi_i_QND<=2*pi);
    for i=1:length(idx)
        new_theta_i_QND(i)=theta_i_QND(idx(i));
        new_phi_i_QND(i)=phi_i_QND(idx(i));
        new_VA_i_QND(i)=VA_i_QND(idx(i));
        new_VB_i_QND(i)=VB_i_QND(idx(i));
    end
    
    idx = find(theta_i_tomo>=val_theta-0.05 & theta_i_tomo<=val_theta+0.05 & phi_i_tomo>=0 & phi_i_tomo<=2*pi);
    for i=1:length(idx)
        new_theta_i_tomo(i)=theta_i_tomo(idx(i));
        new_phi_i_tomo(i)=phi_i_tomo(idx(i));
        new_VA_i_tomo(i)=VA_i_tomo(idx(i));
        new_VB_i_tomo(i)=VB_i_tomo(idx(i));
        Fidelity_i_tomo(i)=Fidelity(rho_i_th(:,:,idx(i)),rho_i_exp(:,:,idx(i)));
    end

    
    figure
    ax = gca;
    plot(new_phi_f_exp_00,Fidelity_f_00,'rx')
    hold on
    plot(new_phi_f_exp_01,Fidelity_f_01,'cx')
    hold on
    plot(new_phi_f_exp_10,Fidelity_f_10,'c^')
    hold on
    plot(new_phi_f_exp_11,Fidelity_f_11,'ro')
    hold on
    plot(new_phi_f_exp_nops,Fidelity_f_nops,'k+')
    hold on
    plot(new_phi_i_tomo,Fidelity_i_tomo,'bo')
    xlabel('\phi')
    ylabel('Fidelity')
    ylim([0 1]); grid on;
    set(gca, 'FontName', 'LM Roman 12', 'FontSize', 16)
    title('Fidelity, measurement V, \theta=3\pi/2')
    legend('final state, 00 postselected','final state, 01 postselected','final state, 10 postselected','final state, 11 postselected','final state, no postselection','initial state')
    saveas(gca,'FVB_i_1D.png')
    saveas(gca,'FVB_i_1D.fig')
    
    figure
    ax = gca;
    plot(new_phi_i_tomo,new_VB_i_tomo,'bo');
    hold on
    plot(new_phi_i_QND,new_VB_i_QND,'b.');
    hold on
    plot(new_phi_f_exp_nops,new_VB_f_exp_nops,'k+');
    hold on
    plot(new_phi_f_exp_00,new_VB_f_exp_00,'rx');
    hold on
    plot(new_phi_f_exp_01,new_VB_f_exp_01,'gv');
    hold on
    plot(new_phi_f_exp_10,new_VB_f_exp_10,'g^');
    hold on
    plot(new_phi_f_exp_11,new_VB_f_exp_11,'ro');
    xlabel('\phi')
    ylabel('V_B')
    ylim([0 1]); grid on;
    set(gca, 'FontName', 'LM Roman 12', 'FontSize', 16)
    title('V_B, \theta=3\pi/2')
    legend('initial tomography','QND','final state, no postselection','final state, 00 postselected','final state, 01 postselected','final state, 10 postselected','final state, 11 postselected')
    saveas(gca,'VB_i_1D.png')
    saveas(gca,'VB_i_1D.fig')
end

%% Predictability

if VPC==2
    % Compute final analytical state
    U=U_P;
    for i=1:length(tomo_i(:,1))

            rho_f_th_16(:,:,i)=U*rho_i_th_16(:,:,i)*U';

            % Postselection
            rho_f_th_00(:,:,i)=(1/trace(rho_f_th_16(:,:,i)*P00)).*PartialTrace(P00*rho_f_th_16(:,:,i)*P00,[3,4],[2,2,2,2]);
            rho_f_th_01(:,:,i)=(1/trace(rho_f_th_16(:,:,i)*P01)).*PartialTrace(P01*rho_f_th_16(:,:,i)*P01,[3,4],[2,2,2,2]);
            rho_f_th_10(:,:,i)=(1/trace(rho_f_th_16(:,:,i)*P10)).*PartialTrace(P10*rho_f_th_16(:,:,i)*P10,[3,4],[2,2,2,2]);
            rho_f_th_11(:,:,i)=(1/trace(rho_f_th_16(:,:,i)*P11)).*PartialTrace(P11*rho_f_th_16(:,:,i)*P11,[3,4],[2,2,2,2]);
            rho_f_th_nops(:,:,i)=PartialTrace(rho_f_th_16(:,:,i),[3,4],[2,2,2,2]);
    end
    
    QND = load('ND_initial_meas_P.txt');
    final1=load('ND_final_meas_P1.txt');
    final4=load('ND_final_meas_P4.txt');
    final_nops=load('ND_final_meas_Pnops.txt');
    
    PA_i_QND=QND(:,4);
    PB_i_QND=QND(:,5);
    
    for i=1:length(tomo_i(:,1))
       PA_i_tomo(i)=predictability('A',rho_i_exp(:,:,i));
       PB_i_tomo(i)=predictability('B',rho_i_exp(:,:,i));
    end

    try
        phi_f_exp_00=final1(:,1);
        theta_f_exp_00=final1(:,3);
    catch
        phi_f_exp_00=[];
        theta_f_exp_00=[];
    end
    
    try
        phi_f_exp_01=final2(:,1);
        theta_f_exp_01=final2(:,3);
    catch
        phi_f_exp_01=[];
        theta_f_exp_01=[];
    end
    
    rho_f_exp_00=read_rho(final1);
    rho_f_exp_11=read_rho(final4);
    rho_f_exp_nops=read_rho(final_nops);
    try
        phi_f_exp_10=final3(:,1);
        theta_f_exp_10=final3(:,3);
    catch
        phi_f_exp_10=[];
        theta_f_exp_10=[];
    end
    try
        phi_f_exp_11=final4(:,1);
        theta_f_exp_11=final4(:,3);
    catch
        phi_f_exp_11=[];
        theta_f_exp_11=[];
    end
    phi_f_exp_nops=final_nops(:,1);
    theta_f_exp_nops=final_nops(:,3);
    phi_i_tomo=tomo_i(:,1);
    theta_i_tomo=tomo_i(:,3);
    phi_i_QND=QND(:,1);
    theta_i_QND=QND(:,3);

    for i=1:length(final1(:,1))
       PA_f_exp_00(i)=predictability('A',rho_f_exp_00(:,:,i));
       PB_f_exp_00(i)=predictability('B',rho_f_exp_00(:,:,i));
    end
    try
        for i=1:length(final2(:,1))
           PA_f_exp_01(i)=predictability('A',rho_f_exp_01(:,:,i));
           PB_f_exp_01(i)=predictability('B',rho_f_exp_01(:,:,i));
        end
    catch
            PA_f_exp_01(i)=NaN;
           PB_f_exp_01(i)=NaN;
    end
    try
        for i=1:length(final3(:,1))
           PA_f_exp_10(i)=predictability('A',rho_f_exp_10(:,:,i));
           PB_f_exp_10(i)=predictability('B',rho_f_exp_10(:,:,i));
        end
    catch
       PA_f_exp_10(i)=NaN;
       PB_f_exp_10(i)=NaN;
    end
    
    for i=1:length(final4(:,1))
       PA_f_exp_11(i)=predictability('A',rho_f_exp_11(:,:,i));
       PB_f_exp_11(i)=predictability('B',rho_f_exp_11(:,:,i));
    end
       
    for i=1:length(final_nops(:,1))
       PA_f_exp_nops(i)=predictability('A',rho_f_exp_nops(:,:,i)); 
       PB_f_exp_nops(i)=predictability('B',rho_f_exp_nops(:,:,i));
    end
    
    val_theta=pi;
    idx = find(theta_f_exp_00>=val_theta-0.05 & theta_f_exp_00<=val_theta+0.05 & phi_f_exp_00>=pi/2 & phi_f_exp_00<=3*pi/2);
    for i=1:length(idx)
        new_theta_f_exp_00(i)=theta_f_exp_00(idx(i));
        new_phi_f_exp_00(i)=phi_f_exp_00(idx(i));
        new_PA_f_exp_00(i)=PA_f_exp_00(idx(i));
        new_PB_f_exp_00(i)=PB_f_exp_00(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_00(i) & phi_i_tomo==new_phi_f_exp_00(i));
        try
        Fidelity_f_00(i)=Fidelity(rho_f_th_00(:,:,idx2),rho_f_exp_00(:,:,idx2(i)));
        catch
            Fidelity_f_00(i)=NaN;
        end
    end
    
    idx = find(theta_f_exp_01>=val_theta-0.05 & theta_f_exp_01<=val_theta+0.05 & phi_f_exp_01>=pi/2 & phi_f_exp_01<=3*pi/2);
    for i=1:length(idx)
        new_theta_f_exp_01(i)=theta_f_exp_01(idx(i));
        new_phi_f_exp_01(i)=phi_f_exp_01(idx(i));
        new_PA_f_exp_01(i)=PA_f_exp_01(idx(i));
        new_PB_f_exp_01(i)=PB_f_exp_01(idx(i));
        idx2=find(theta_i_tomo==new_theta_f_exp_01(i) & phi_i_tomo==new_phi_f_exp_01(i));
        try
            Fidelity_f_01(i)=Fidelity(rho_f_th_01(:,:,idx2),rho_f_exp_01(:,:,idx(i)));
        catch
            Fidelity_f_01(i)=NaN;
        end
    end
    
    idx = find(theta_f_exp_10>=val_theta-0.05 & theta_f_exp_10<=val_theta+0.05 & phi_f_exp_10>=pi/2 & phi_f_exp_10<=3*pi/2);
    for i=1:length(idx)
        new_theta_f_exp_10(i)=theta_f_exp_10(idx(i));
        new_phi_f_exp_10(i)=phi_f_exp_10(idx(i));
        new_PA_f_exp_10(i)=PA_f_exp_10(idx(i));
        new_PB_f_exp_10(i)=PB_f_exp_10(idx(i));
        
        idx2=find(theta_i_tomo==new_theta_f_exp_10(i) & phi_i_tomo==new_phi_f_exp_10(i));
        try
            Fidelity_f_10(i)=Fidelity(rho_f_th_10(:,:,idx2),rho_f_exp_10(:,:,idx(i)));
        catch
            Fidelity_f_10(i)=NaN;
        end
    end
    
    idx = find(theta_f_exp_11>=val_theta-0.05 & theta_f_exp_11<=val_theta+0.05 & phi_f_exp_11>=pi/2 & phi_f_exp_11<=3*pi/2);
    for i=1:length(idx)
        new_theta_f_exp_11(i)=theta_f_exp_11(idx(i));
        new_phi_f_exp_11(i)=phi_f_exp_11(idx(i));
        new_PA_f_exp_11(i)=PA_f_exp_11(idx(i));
        new_PB_f_exp_11(i)=PB_f_exp_11(idx(i));
        
        idx2=find(theta_i_tomo==new_theta_f_exp_11(i) & phi_i_tomo==new_phi_f_exp_11(i));
        try
            Fidelity_f_11(i)=Fidelity(rho_f_th_11(:,:,idx2),rho_f_exp_11(:,:,idx(i)));
        catch
            Fidelity_f_11(i)=NaN;
        end
    end
    
    idx = find(theta_f_exp_nops>=val_theta-0.05 & theta_f_exp_nops<=val_theta+0.05 & phi_f_exp_nops>=pi/2 & phi_f_exp_nops<=3*pi/2);
    for i=1:length(idx)
        new_theta_f_exp_nops(i)=theta_f_exp_nops(idx(i));
        new_phi_f_exp_nops(i)=phi_f_exp_nops(idx(i));
        new_PA_f_exp_nops(i)=PA_f_exp_nops(idx(i));
        new_PB_f_exp_nops(i)=PB_f_exp_nops(idx(i));

        idx2=find(theta_i_tomo==new_theta_f_exp_nops(i) & phi_i_tomo==new_phi_f_exp_nops(i));
        try
            Fidelity_f_nops(i)=Fidelity(rho_f_th_nops(:,:,idx2),rho_f_exp_nops(:,:,idx(i)));
        catch
            Fidelity_f_nops(i)=NaN;
        end
    end
    
    idx = find(theta_i_QND>=val_theta-0.05 & theta_i_QND<=val_theta+0.05 & phi_i_QND>=pi/2 & phi_i_QND<=3*pi/2);
    for i=1:length(idx)
    new_theta_i_QND(i)=theta_i_QND(idx(i));
    new_phi_i_QND(i)=phi_i_QND(idx(i));
    new_PA_i_QND(i)=PA_i_QND(idx(i));
    new_PB_i_QND(i)=PB_i_QND(idx(i));
    end
    
    idx = find(theta_i_tomo>=val_theta-0.05 & theta_i_tomo<=val_theta+0.05 & phi_i_tomo>=pi/2 & phi_i_tomo<=3*pi/2);
    for i=1:length(idx)
    new_theta_i_tomo(i)=theta_i_tomo(idx(i));
    new_phi_i_tomo(i)=phi_i_tomo(idx(i));
    new_PA_i_tomo(i)=PA_i_tomo(idx(i));
    new_PB_i_tomo(i)=PB_i_tomo(idx(i));
    Fidelity_i_tomo(i)=Fidelity(rho_i_th(:,:,idx(i)),rho_i_exp(:,:,idx(i)));
    end

    figure
    ax = gca;
    plot(new_phi_i_tomo,new_PA_i_tomo,'b*');
    hold on
    plot(new_phi_i_QND,new_PA_i_QND,'b.');
    hold on
    plot(new_phi_f_exp_nops,new_PA_f_exp_nops,'k+');
    hold on
    plot(new_phi_f_exp_00,new_PA_f_exp_00,'rx');
    hold on
    plot(new_phi_f_exp_11,new_PA_f_exp_11,'ro');
    xlabel('\phi')
    ylabel('P_A')
    ylim([0 1]); grid on;
    set(gca, 'FontName', 'LM Roman 12', 'FontSize', 16)
    title('P_A, \theta=\pi')
    legend('initial tomography','QND','final state, no postselection','final state, 00 postselected','final state, 11 postselected')
    saveas(gca,'PA_i_1D.png')
    saveas(gca,'PA_i_1D.fig')

    figure
    ax = gca;
    plot(new_phi_i_tomo,new_PB_i_tomo,'b*');
    hold on
    plot(new_phi_i_QND,new_PB_i_QND,'b.');
    hold on
    plot(new_phi_f_exp_nops,new_PB_f_exp_nops,'k+');
    hold on
    plot(new_phi_f_exp_00,new_PB_f_exp_00,'rx');
    hold on
    plot(new_phi_f_exp_11,new_PB_f_exp_11,'ro');
    xlabel('\phi')
    ylabel('P_B')
    ylim([0 1]); grid on;
    set(gca, 'FontName', 'LM Roman 12', 'FontSize', 16)
    title('P_B, \theta=\pi')
    legend('initial tomography','QND','final state, no postselection','final state, 00 postselected','final state, 11 postselected')
    saveas(gca,'PB_i_1D.png')
    saveas(gca,'PB_i_1D.fig')

    figure
    ax = gca;
    plot(new_phi_f_exp_00,Fidelity_f_00,'rx')
    hold on
    plot(new_phi_f_exp_11,Fidelity_f_11,'ro')
    hold on
    plot(new_phi_f_exp_nops,Fidelity_f_nops,'k+')
    hold on
    plot(new_phi_i_tomo,Fidelity_i_tomo,'bo')
    xlabel('\phi')
    ylabel('Fidelity')
    ylim([0 1]); grid on;
    set(gca, 'FontName', 'LM Roman 12', 'FontSize', 16)
    title('Fidelity, measurement P, \theta=\pi')
    legend('final state, 00 postselected','final state, 11 postselected','final state, no postselection','initial state')
    saveas(gca,'FP_i_1D.png')
    saveas(gca,'FP_i_1D.fig')
end