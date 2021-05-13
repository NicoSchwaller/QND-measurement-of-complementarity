# Needed for functions
import numpy as np
from copy import deepcopy
from collections import Counter

# Transpiling stuff
from qiskit.compiler import transpile
from qiskit.providers.aer.noise import NoiseModel

# Import Qiskit classes
import qiskit
from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister, IBMQ, execute

# Tomography functions
from qiskit.ignis.verification.tomography import state_tomography_circuits, StateTomographyFitter

# Import measurement calibration functions
from qiskit.ignis.mitigation.measurement import (complete_meas_cal, CompleteMeasFitter)

# In[1]: Setting account and experimental parameters

IBMQ.enable_account('token')
provider=IBMQ.get_provider(hub='ibm-q', group='open', project='main')

# Parameters of the experiment ---------------------------------------------------------------------------------------

# Choose backend
backend=provider.get_backend('ibmqx2')

# Number of shots per circuit
shots_tomo_f=5000
shots_tomo_i=5000
shots_QND=5000

# Minimal number of outcomes above which tomography is performed
mc=50

# Increment of angle phi
delta_phi=np.pi/32

optimization=0
noise_simulation=1

circuit_types=['simple','complete']


# In[2]: Models and functions

# Generate the calibration circuits
qr = qiskit.QuantumRegister(4)
meas_calibs, state_labels = complete_meas_cal(qubit_list=[0,1,2,3], qr=qr, circlabel='mcal')

# Execute the calibration circuits without noise
simulator = qiskit.Aer.get_backend('qasm_simulator')
job = qiskit.execute(meas_calibs, backend=simulator, shots=8192)
cal_results = job.result()

# The calibration matrix without noise is the identity matrix
meas_fitter = CompleteMeasFitter(cal_results, state_labels, circlabel='mcal')

#Select backend to compute noise model
device = backend


#Create noise model
properties = device.properties()
coupling_map = device.configuration().coupling_map
noise_model = NoiseModel.from_backend(backend)
basis_gates = noise_model.basis_gates

# Execute the calibration circuits using the noise model
noisy_job = qiskit.execute(meas_calibs, backend=simulator, shots=18192, noise_model=noise_model)
results = noisy_job.result()

# Calculate the calibration matrix
meas_fitter = CompleteMeasFitter(results, state_labels, circlabel='mcal')
meas_filter = meas_fitter.filter

# What is the measurement fidelity?
print("Average Measurement Fidelity: %f" % meas_fitter.readout_fidelity())

def partial_tracey(rho, qubit_2_keep):
    """ Calculate the partial trace for qubit system
    Parameters
    ----------
    rho: np.ndarray
        Density matrix
    qubit_2_keep: list
        Index of qubit to be kept after taking the trace
    Returns
    -------
    rho_res: np.ndarray
        Density matrix after taking partial trace
    """
    num_qubit = int(np.log2(rho.shape[0]))
    qubit_axis = [(i, num_qubit + i) for i in range(num_qubit)
                  if i not in qubit_2_keep]
    minus_factor = [(i, 2 * i) for i in range(len(qubit_axis))]
    minus_qubit_axis = [(q[0] - m[0], q[1] - m[1])
                        for q, m in zip(qubit_axis, minus_factor)]
    rho_res = np.reshape(rho, [2, 2] * num_qubit)
    qubit_left = num_qubit - len(qubit_axis)
    for i, j in minus_qubit_axis:
        rho_res = np.trace(rho_res, axis1=i, axis2=j)
    if qubit_left > 1:
        rho_res = np.reshape(rho_res, [2 ** qubit_left] * 2)
    
    return rho_res

def output_qubits(circuit):
    indices = slice(0, 2)
    return circuit.qubits[indices]

def tomography_processing(job, circuit):
    raw_results = job.result()
    new_results = deepcopy(raw_results)

    # Aggregate the results to erase information about all classical registers
    # except the one used for tomography
    for result_i, result in enumerate(new_results.results):
        # The tomography creg is the first one in the list
        creg_label, n_tomography_bits = result.header.creg_sizes[0]

        # Modify the result to forget all other cregs
        result.header.creg_sizes = [[creg_label, n_tomography_bits]]
        result.header.clbit_labels =             result.header.clbit_labels[:n_tomography_bits]
        result.header.memory_slots = n_tomography_bits

        # Aggregate counts
        old_counts = raw_results.get_counts(result_i)
        new_counts = Counter()

        for outcome, count in old_counts.items():
            new_counts[outcome.split()[-1]] += count

        new_results.results[result_i].data.counts = result.data.counts.from_dict(dict(new_counts))

    # An empty circuit with the same number of quantum registers but no
    # classical registers
    outputs = output_qubits(circuit)
    fake_circuit = QuantumCircuit(*circuit.qregs)
    fake_tomographs = state_tomography_circuits(fake_circuit, outputs)

    return StateTomographyFitter(new_results, fake_tomographs).fit()

def line(rho):
    rho_line=str(rho[0,0]).strip('(').strip(')')+' '+str(rho[0,1]).strip('(').strip(')')+' '+str(rho[0,2]).strip('(').strip(')')+' '+str(rho[0,3]).strip('(').strip(')')+' '+str(rho[1,0]).strip('(').strip(')')+' '+str(rho[1,1]).strip('(').strip(')')+' '+str(rho[1,2]).strip('(').strip(')')+' '+str(rho[1,3]).strip('(').strip(')')+' '+str(rho[2,0]).strip('(').strip(')')+' '+str(rho[2,1]).strip('(').strip(')')+' '+str(rho[2,2]).strip('(').strip(')')+' '+str(rho[2,3]).strip('(').strip(')')+' '+str(rho[3,0]).strip('(').strip(')')+' '+str(rho[3,1]).strip('(').strip(')')+' '+str(rho[3,2]).strip('(').strip(')')+' '+str(rho[3,3]).strip('(').strip(')')
    return rho_line


# In[3]: General circuit performing the QND measurement of V, P, and C

# Output state tomography
# if block_QND_tomo is set to 0, the tomography of the output state is performed after the QND circuit
# if block_QND_tomo is set to 1, only the QND measurement of the state is performed
block_QND_tomo=0

# Path of folder in which to save results
path="C:/Users/Nico/Documents/QND_complementarity/circuit2"

# Iteration: tomography, QND measurement of V, P and C
for p in range(1,5):
  if p==1:
      nondemolition=0
      print('initial tomography, compute VPC')
      VPC=''
  if p==2:
      VPC=1
      nondemolition=1
      print('QND measurement of V')
  if p==3:
      VPC=2
      nondemolition=1
      print('QND measurement of P')
  if p==4:
      VPC=3
      nondemolition=1
      print('QND measurement of C')

  qr = QuantumRegister(4)
  qr2 = QuantumRegister(3)
  cr = ClassicalRegister(3)

  # -----------------------------------------------------
  # Choosing the states to analyze (range of angles and increments)
  if VPC==1:
      delta_theta=3*np.pi/2
      theta_start=0
      theta_end=3*np.pi/2
      phi_start=0
      phi_end=2*np.pi
  if VPC==2:
      delta_theta=np.pi
      theta_start=np.pi
      theta_end=np.pi
      phi_start=np.pi/2
      phi_end=3*np.pi/2
  if VPC==3:
      delta_theta=np.pi
      theta_start=np.pi
      theta_end=np.pi
      phi_start=0*np.pi/32
      phi_end=np.pi
  if p==1:
      theta_start=0*np.pi
      theta_end=3*np.pi/2
      delta_theta=1*np.pi/2
      phi_start=0*np.pi
      phi_end=2*np.pi
  delta_var_lambda=np.pi/16
  var_lambda_start=0
  var_lambda_end=0

  if nondemolition==0:
    vec_C_AB_initial_tomo=[]
    vec_VA_initial_tomo=[]
    vec_PA_initial_tomo=[]
    vec_VB_initial_tomo=[]
    vec_PB_initial_tomo=[]
    logi=open(path+"/tomo_initial.txt","w")
  elif nondemolition==1: 
    if VPC==3:
      logND=open(path+"/ND_initial_meas_C.txt","w")
      logfnops=open(path+"/ND_final_meas_Cnops.txt","w")
      logf1=open(path+"/ND_final_meas_C1.txt","w")
      logf2=open(path+"/ND_final_meas_C2.txt","w")
      logf3=open(path+"/ND_final_meas_C3.txt","w")
      logf4=open(path+"/ND_final_meas_C4.txt","w")
    if VPC==1:
      logND=open(path+"/ND_initial_meas_V.txt","w")
      logfnops=open(path+"/ND_final_meas_Vnops.txt","w")
      logf1=open(path+"/ND_final_meas_V1.txt","w")
      logf2=open(path+"/ND_final_meas_V2.txt","w")
      logf3=open(path+"/ND_final_meas_V3.txt","w")
      logf4=open(path+"/ND_final_meas_V4.txt","w")
    if VPC==2:
      logND=open(path+"/ND_initial_meas_P.txt","w")
      logfnops=open(path+"/ND_final_meas_Pnops.txt","w")
      logf1=open(path+"/ND_final_meas_P1.txt","w")
      logf2=open(path+"/ND_final_meas_P2.txt","w")
      logf3=open(path+"/ND_final_meas_P3.txt","w")
      logf4=open(path+"/ND_final_meas_P4.txt","w")

  # Scan lambda
  for i in range(1,round((var_lambda_end-var_lambda_start)/delta_var_lambda)+2):
    var_lambda=var_lambda_start+i*delta_var_lambda
    print('i=',i,'/',round((var_lambda_end-var_lambda_start)/delta_var_lambda)+1,'---------------------')
    
    # Scan theta
    for j in range(1,round((theta_end-theta_start)/delta_theta)+2):
        print('j-1=',j-1,'/',round((theta_end-theta_start)/delta_theta)+1)
        theta=theta_start+(j-1)*delta_theta

        # Scan phi
        for k in range(1,round((phi_end-phi_start)/delta_phi)+2):
          phi=phi_start+(k-1)*delta_phi
          print('k=',k)
          circuit = QuantumCircuit(qr)

          # Input state -------------------------------------------------------------------------------------------------------------
          circuit.ry(phi,qr[0])
          circuit.cu3(theta,0,0,qr[0],qr[1])
          # -------------------------------------------------------------------------------------------------------------------------------

          if nondemolition==1:

            # Nondemolition circuit to measure V, P or C
            if VPC==1:
              circuit.ry(np.pi/2,qr[0])
              circuit.ry(np.pi/2,qr[1])
            if VPC==3:
              circuit.rx(np.pi/2,qr[0])
              circuit.rx(np.pi/2,qr[1])
              circuit.ry(np.pi/2,qr[2])
              circuit.cx(qr[2],qr[3])
            circuit.cx(qr[0],qr[2])
            circuit.cx(qr[1],qr[3])
            if VPC==1:
              circuit.ry(-np.pi/2,qr[0])
              circuit.ry(-np.pi/2,qr[1])
            if VPC==3:
              circuit.rx(-np.pi/2,qr[0])
              circuit.rx(-np.pi/2,qr[1])
          
          if nondemolition==1:
            circuit_t=deepcopy(circuit)

            # ---------------
            # QND measurement 
            # ---------------

            crc = ClassicalRegister(2)
            circuit.add_register(crc)
            circuit.measure(qr[2],crc[0])
            circuit.measure(qr[3],crc[1])
          
            # transpiled circuit
            transp_circ = transpile(circuit, backend=backend, seed_transpiler=11)

            # optimized circuit
            if optimization==1:
                circuit = transpile(circuit, backend=backend, seed_transpiler=11, optimization_level=3)
            elif optimization==0:
                circuit = transpile(circuit, backend=backend, seed_transpiler=11)
                
            # noise simulation
            if noise_simulation==1:
                job_counts = execute(circuit, simulator, noise_model=noise_model,coupling_map=coupling_map,basis_gates=basis_gates, shots=shots_QND)
            elif noise_simulation==0:
                job_counts = execute(circuit, backend, shots=shots_QND)
            result_counts=job_counts.result()
            counts=result_counts.get_counts()

            # Compute probabilities from measurement outcomes
            p_00=0
            p_01=0
            p_10=0
            p_11=0
            for state in counts:
                  if state=='00':
                      p_00=p_00+counts[state]/shots_QND
                  if state=='01':
                      p_10=p_10+counts[state]/shots_QND
                  if state=='10':
                      p_01=p_01+counts[state]/shots_QND
                  if state=='11':
                      p_11=p_11+counts[state]/shots_QND
                        
            # QND computation from the outcomes' statistics and registering characteristics of the circuits (depth, number of gates, ...)
            if VPC==3:
                  C_AB_initial=abs(p_10+p_01-p_00-p_11)

                  # Characteristics of the circuits
                  used_gates = transp_circ.count_ops()
                  tgates={'u1':0,'u2':0,'u3':0,'cx':0}
                  for gate in ('u1','u2','u3','cx'):
                      if gate in used_gates:
                          tgates[gate]=used_gates[gate]
                  tdepth = transp_circ.depth()
                  tfilename='transp_circ_'+str(round(phi,1))+'_'+str(round(theta,1))+'_C.png'

                  used_gates = circuit.count_ops() 
                  gates={'u1':0,'u2':0,'u3':0,'cx':0}
                  for gate in ('u1','u2','u3','cx'):
                      if gate in used_gates:
                          gates[gate]=used_gates[gate]
                  depth = circuit.depth()
                  filename ='circuit_'+str(round(phi,1))+'_'+str(round(theta,1))+'_C.png'

                  # String to be added in output file
                  string_i=str(phi) + ' ' +str(var_lambda) + ' ' +str(theta)+ ' '+str(C_AB_initial)+' '+str(depth)+' '+str(gates['u1'])+' '+str(gates['u2'])+' '+str(gates['u3'])+' '+str(gates['cx'])+' '+str(tdepth)+' '+str(tgates['u1'])+' '+str(tgates['u2'])+' '+str(tgates['u3'])+' '+str(tgates['cx'])+'\n'

            if VPC==2:
                  PA_initial=abs(p_00+p_01-p_10-p_11)
                  PB_initial=abs(p_00+p_10-p_01-p_11)

                  # Characteristics of the circuits
                  used_gates = transp_circ.count_ops()
                  tgates={'u1':0,'u2':0,'u3':0,'cx':0}
                  for gate in ('u1','u2','u3','cx'):
                      if gate in used_gates:
                          tgates[gate]=used_gates[gate]
                  tdepth = transp_circ.depth()
                  tfilename='transp_circ_'+str(round(phi,1))+'_'+str(round(theta,1))+'_P.png'

                  used_gates = circuit.count_ops() 
                  gates={'u1':0,'u2':0,'u3':0,'cx':0}
                  for gate in ('u1','u2','u3','cx'):
                      if gate in used_gates:
                          gates[gate]=used_gates[gate]
                  depth = circuit.depth()
                  filename ='circuit_'+str(round(phi,1))+'_'+str(round(theta,1))+'_P.png'

                  # String to be added in output file
                  string_i=str(phi) + ' ' +str(var_lambda) + ' ' +str(theta)+ ' '+str(PA_initial)+ ' '+str(PB_initial)+' '+str(depth)+' '+str(gates['u1'])+' '+str(gates['u2'])+' '+str(gates['u3'])+' '+str(gates['cx'])+' '+str(tdepth)+' '+str(tgates['u1'])+' '+str(tgates['u2'])+' '+str(tgates['u3'])+' '+str(tgates['cx'])+'\n'
            if VPC==1:
                  VA_initial=abs(p_00+p_01-p_10-p_11)
                  VB_initial=abs(p_00+p_10-p_01-p_11)

                  # Characteristics of the circuits
                  used_gates = transp_circ.count_ops()
                  tgates={'u1':0,'u2':0,'u3':0,'cx':0}
                  for gate in ('u1','u2','u3','cx'):
                      if gate in used_gates:
                          tgates[gate]=used_gates[gate]
                  tdepth = transp_circ.depth()
                  tfilename='transp_circ_'+str(round(phi,1))+'_'+str(round(theta,1))+'_V.png'

                  used_gates = circuit.count_ops() 
                  gates={'u1':0,'u2':0,'u3':0,'cx':0}
                  for gate in ('u1','u2','u3','cx'):
                      if gate in used_gates:
                          gates[gate]=used_gates[gate]
                  depth = circuit.depth()
                  filename =path+'circuit_'+str(round(phi,1))+'_'+str(round(theta,1))+'_V.png'

                  # String to be added in output file
                  string_i=str(phi) + ' ' +str(var_lambda) + ' ' +str(theta)+ ' '+str(VA_initial)+ ' '+str(VB_initial)+' '+str(depth)+' '+str(gates['u1'])+' '+str(gates['u2'])+' '+str(gates['u3'])+' '+str(gates['cx'])+' '+str(tdepth)+' '+str(tgates['u1'])+' '+str(tgates['u2'])+' '+str(tgates['u3'])+' '+str(tgates['cx'])+'\n'
            logND.write(string_i)

            # --------------------------------
            # QND measurement and tomography of the outgoing state 
            # --------------------------------

            if block_QND_tomo==0:

                  string_VP_fnops=''
                  string_VP_f1=''
                  string_VP_f2=''
                  string_VP_f3=''
                  string_VP_f4=''
                    
                  qst_VDC = state_tomography_circuits(circuit_t,[qr[0], qr[1]])
                  qst_VDC_no_anc = deepcopy(qst_VDC)

                  ca = ClassicalRegister(2)
                  for qst_VDC_circ in qst_VDC:
                    qst_VDC_circ.add_register(ca)
                    qst_VDC_circ.measure(qr[2],ca[0])
                    qst_VDC_circ.measure(qr[3],ca[1])

                  if noise_simulation==1:
                    job = execute(qst_VDC, simulator, noise_model=noise_model,coupling_map=coupling_map,basis_gates=basis_gates, shots=shots_tomo_f)
                    job_nops=execute(qst_VDC_no_anc, simulator, noise_model=noise_model,coupling_map=coupling_map,basis_gates=basis_gates, shots=shots_tomo_f)
                  elif noise_simulation==0:
                    job = execute(qst_VDC, backend, shots=shots_tomo_f)
                    job_nops = execute(qst_VDC_no_anc, backend, shots=shots_tomo_f)
                    
                  raw_results=job.result()
                  resultnops = job_nops.result()
                  new_result1 = deepcopy(raw_results)
                  new_result2 = deepcopy(raw_results)
                  new_result3 = deepcopy(raw_results)
                  new_result4 = deepcopy(raw_results)

                  # Retrieve outcomes of the ancillae qubits
                  for resultidx, _ in enumerate(raw_results.results):
                    old_counts = raw_results.get_counts(resultidx)
                    new_counts1 = {}
                    new_counts2 = {}
                    new_counts3 = {}
                    new_counts4 = {}
                    
                    new_result1.results[resultidx].header.creg_sizes = [new_result1.results[resultidx].header.creg_sizes[0]]
                    new_result1.results[resultidx].header.clbit_labels = new_result1.results[resultidx].header.clbit_labels[0:-1]
                    new_result1.results[resultidx].header.memory_slots = 2
                    new_result2.results[resultidx].header.creg_sizes = [new_result2.results[resultidx].header.creg_sizes[0]]
                    new_result2.results[resultidx].header.clbit_labels = new_result2.results[resultidx].header.clbit_labels[0:-1]
                    new_result2.results[resultidx].header.memory_slots = 2
                    new_result3.results[resultidx].header.creg_sizes = [new_result3.results[resultidx].header.creg_sizes[0]]
                    new_result3.results[resultidx].header.clbit_labels = new_result3.results[resultidx].header.clbit_labels[0:-1]
                    new_result3.results[resultidx].header.memory_slots = 2
                    new_result4.results[resultidx].header.creg_sizes = [new_result4.results[resultidx].header.creg_sizes[0]]
                    new_result4.results[resultidx].header.clbit_labels = new_result4.results[resultidx].header.clbit_labels[0:-1]
                    new_result4.results[resultidx].header.memory_slots = 2

                    for reg_key in old_counts:
                        reg_bits = reg_key.split(' ')
                        if VPC==1 or VPC==2 or VPC==3:
                          if reg_bits[0]=='00':
                            new_counts1[reg_bits[1]]=old_counts[reg_key]
                          if reg_bits[0]=='10':
                            new_counts2[reg_bits[1]]=old_counts[reg_key]
                          if reg_bits[0]=='01':
                            new_counts3[reg_bits[1]]=old_counts[reg_key]
                          if reg_bits[0]=='11':
                            new_counts4[reg_bits[1]]=old_counts[reg_key]
                            
                    new_result1.results[resultidx].data.counts = new_result1.results[resultidx].data.from_dict(new_counts1).to_dict()
                    new_result2.results[resultidx].data.counts = new_result2.results[resultidx].data.from_dict(new_counts2).to_dict()
                    new_result3.results[resultidx].data.counts = new_result3.results[resultidx].data.from_dict(new_counts3).to_dict()
                    new_result4.results[resultidx].data.counts = new_result4.results[resultidx].data.from_dict(new_counts4).to_dict()

                  tomo_VDCnops = StateTomographyFitter(resultnops, qst_VDC_no_anc)  
                  tomo_VDC1 = StateTomographyFitter(new_result1, qst_VDC_no_anc)
                  tomo_VDC2 = StateTomographyFitter(new_result2, qst_VDC_no_anc)
                  tomo_VDC3 = StateTomographyFitter(new_result3, qst_VDC_no_anc)
                  tomo_VDC4 = StateTomographyFitter(new_result4, qst_VDC_no_anc)

                  rhonops = tomo_VDCnops.fit()
                  rho_ABnops=partial_tracey(rhonops,[1,0])
                  string_VP_fnops=string_VP_f1+str(phi) + ' ' +str(var_lambda) + ' ' +str(theta)+' '+line(rho_ABnops)+'\n'

                  if len(new_counts1) != 0 and not(len(new_counts1)==1 and new_counts1.get(list(new_counts1.keys())[0]) < mc) and not(len(new_counts1)==2 and new_counts1.get(list(new_counts1.keys())[0])+new_counts1.get(list(new_counts1.keys())[1]) < mc) and not(len(new_counts1)==3 and new_counts1.get(list(new_counts1.keys())[0])+new_counts1.get(list(new_counts1.keys())[1])+new_counts1.get(list(new_counts1.keys())[2]) < mc) and not(len(new_counts1)==4 and new_counts1.get(list(new_counts1.keys())[0])+new_counts1.get(list(new_counts1.keys())[1])+new_counts1.get(list(new_counts1.keys())[2])+new_counts1.get(list(new_counts1.keys())[3]) < mc):
                    #print('counts1=',new_counts1)
                    rho1 = tomo_VDC1.fit()
                    rho_AB1=partial_tracey(rho1,[1,0])
                    string_VP_f1=string_VP_f1+str(phi) + ' ' +str(var_lambda) + ' ' +str(theta)+' '+line(rho_AB1)+'\n'

                  if len(new_counts2) != 0 and not(len(new_counts2)==1 and new_counts2.get(list(new_counts2.keys())[0]) < mc) and not(len(new_counts2)==2 and new_counts2.get(list(new_counts2.keys())[0])+new_counts2.get(list(new_counts2.keys())[1]) < mc) and not(len(new_counts2)==3 and new_counts2.get(list(new_counts2.keys())[0])+new_counts2.get(list(new_counts2.keys())[1])+new_counts2.get(list(new_counts2.keys())[2]) < mc) and not(len(new_counts2)==4 and new_counts2.get(list(new_counts2.keys())[0])+new_counts2.get(list(new_counts2.keys())[1])+new_counts2.get(list(new_counts2.keys())[2])+new_counts2.get(list(new_counts2.keys())[3]) < mc):
                    rho2 = tomo_VDC2.fit()
                    rho_AB2=partial_tracey(rho2,[1,0])
                    string_VP_f2=string_VP_f2+str(phi) + ' ' +str(var_lambda) + ' ' +str(theta)+' '+line(rho_AB2)+'\n'

                  if len(new_counts3) != 0 and not(len(new_counts3)==1 and new_counts3.get(list(new_counts3.keys())[0]) < mc) and not(len(new_counts3)==2 and new_counts3.get(list(new_counts3.keys())[0])+new_counts3.get(list(new_counts3.keys())[1]) < mc) and not(len(new_counts3)==3 and new_counts3.get(list(new_counts3.keys())[0])+new_counts3.get(list(new_counts3.keys())[1])+new_counts3.get(list(new_counts3.keys())[2]) < mc) and not(len(new_counts3)==4 and new_counts3.get(list(new_counts3.keys())[0])+new_counts3.get(list(new_counts3.keys())[1])+new_counts3.get(list(new_counts3.keys())[2])+new_counts3.get(list(new_counts3.keys())[3]) < mc):
                    rho3 = tomo_VDC3.fit()
                    rho_AB3=partial_tracey(rho3,[1,0])
                    string_VP_f3=string_VP_f3+str(phi) + ' ' +str(var_lambda) + ' ' +str(theta)+' '+line(rho_AB3)+'\n'

                  if len(new_counts4) != 0 and not(len(new_counts4)==1 and new_counts4.get(list(new_counts4.keys())[0]) < mc) and not(len(new_counts4)==2 and new_counts4.get(list(new_counts4.keys())[0])+new_counts4.get(list(new_counts4.keys())[1]) < mc) and not(len(new_counts4)==3 and new_counts4.get(list(new_counts4.keys())[0])+new_counts4.get(list(new_counts4.keys())[1])+new_counts4.get(list(new_counts4.keys())[2]) < mc) and not(len(new_counts4)==4 and new_counts4.get(list(new_counts4.keys())[0])+new_counts4.get(list(new_counts4.keys())[1])+new_counts4.get(list(new_counts4.keys())[2])+new_counts4.get(list(new_counts4.keys())[3]) < mc):
                    rho4 = tomo_VDC4.fit()
                    rho_AB4=partial_tracey(rho4,[1,0])
                    string_VP_f4=string_VP_f4+str(phi) + ' ' +str(var_lambda) + ' ' +str(theta)+' '+line(rho_AB4)+'\n'

                  logfnops.write(string_VP_fnops)
                  logf1.write(string_VP_f1)
                  logf2.write(string_VP_f2)
                  logf3.write(string_VP_f3)
                  logf4.write(string_VP_f4)

          # -----------------------------
          # Tomography of the input state
          # -----------------------------

          if nondemolition==0:
            qst_VDC = state_tomography_circuits(circuit,[qr[0],qr[1]])
            if noise_simulation==1:
                job = execute(qst_VDC, simulator, noise_model=noise_model,coupling_map=coupling_map,basis_gates=basis_gates, shots=shots_tomo_i)
            elif noise_simulation==0:
                job = execute(qst_VDC, backend, shots=shots_tomo_i)
            tomo_VDC = StateTomographyFitter(job.result(), qst_VDC)
            rho = tomo_VDC.fit()
            string_it=str(phi) + ' ' +str(var_lambda) + ' ' +str(theta)+' '+line(rho)+'\n'
            logi.write(string_it)

  if nondemolition==1:
    logND.close()
    logfnops.close()
    logf1.close()
    logf2.close()
    logf3.close()
    logf4.close()

  if nondemolition==0:
    logi.close()


# In[4]: Simple circuit performing the QND measurement of C

# Output state tomography
# if block_QND_tomo is set to 0, the tomography of the output state is performed after the QND circuit
# if block_QND_tomo is set to 1, only the QND measurement of the state is performed
block_QND_tomo=0

# Path of folder in which to save results
path="C:/Users/Nico/Documents/QND_complementarity/circuit1"

# Measurement of concurrence
print('QND measurement of C')

qr = QuantumRegister(3)
cr = ClassicalRegister(3)

delta_theta=np.pi
theta_depart=np.pi
theta_fin=np.pi
phi_depart=0*np.pi
phi_fin=np.pi
delta_var_lambda=np.pi/16
var_lambda_depart=0
var_lambda_fin=0

logND=open(path+"/ND_initial_meas_C.txt","w")
logfnops=open(path+"/ND_final_meas_Cnops.txt","w")
logf1=open(path+"/ND_final_meas_C1.txt","w")
logf2=open(path+"/ND_final_meas_C2.txt","w")


for i in range(1,round((var_lambda_fin-var_lambda_depart)/delta_var_lambda)+2):
    var_lambda=var_lambda_depart+i*delta_var_lambda
    print('i=',i,'/',round((var_lambda_fin-var_lambda_depart)/delta_var_lambda)+1,'---------------------')
    
    for j in range(1,round((theta_fin-theta_depart)/delta_theta)+2):
        print('j-1=',j-1,'/',round((theta_fin-theta_depart)/delta_theta)+1)
        theta=theta_depart+(j-1)*delta_theta
    
        for k in range(1,round((phi_fin-phi_depart)/delta_phi)+2):
          phi=phi_depart+(k-1)*delta_phi
          print('k=',k)
    
          circuit = QuantumCircuit(qr)
    
          # Input state -------------------------------------------------------------------------------------------------------------
    
          circuit.ry(phi,qr[0])
          circuit.cu3(theta,0,0,qr[0],qr[1])
    
          # -------------------------------------------------------------------------------------------------------------------------------
    
          circuit.barrier()
    
          # Nondemolition circuit to measure VPC
          circuit.rx(np.pi/2,qr[0])
          circuit.rx(np.pi/2,qr[1])
          circuit.cx(qr[0],qr[2])
          circuit.cx(qr[1],qr[2])
          circuit.rx(-np.pi/2,qr[0])
          circuit.rx(-np.pi/2,qr[1])
    
          # -------------------------------------------------------------------------------------------------------
          # Nondemolition measurement of the input state, and tomography of the outgoing state
          # -------------------------------------------------------------------------------------------------------
          circuit_t=deepcopy(circuit)
    
          # ---------------
          # QND measurement 
          # ---------------
    
          crc = ClassicalRegister(1)
          circuit.add_register(crc)
          circuit.measure(qr[2],crc[0])
    
          # Transpiled circuit
          transp_circ = transpile(circuit, backend=backend, seed_transpiler=11)
    
          # Optimized circuit
          if optimization==1:
                circuit = transpile(circuit, backend=backend, seed_transpiler=11, optimization_level=3)
          elif optimization==0:
                circuit = transpile(circuit, backend=backend, seed_transpiler=11)
                
          # Noise simulation  
          if noise_simulation==1:
                job_counts = execute(circuit, simulator, noise_model=noise_model,coupling_map=coupling_map,basis_gates=basis_gates, shots=shots_QND)       
          elif noise_simulation==0:
                job_counts = execute(circuit, backend, shots=shots_QND)
    
          result_counts=job_counts.result()
          counts=result_counts.get_counts()
            
          # Convert outcomes into probabilities
          p_0=0
          p_1=0
          for state in counts:
                  if state=='0':
                      p_0=p_0+counts[state]/shots_QND
                  if state=='1':
                      p_1=p_1+counts[state]/shots_QND
    
          # Compute C from probabilities              
          C_AB_initial=abs(p_0-p_1)
    
          # Characteristics of the circuit
          used_gates = transp_circ.count_ops()
          tgates={'u1':0,'u2':0,'u3':0,'cx':0}
          for gate in ('u1','u2','u3','cx'):
              if gate in used_gates:
                  tgates[gate]=used_gates[gate]
          tdepth = transp_circ.depth()
          tfilename='transp_circ_'+str(round(phi,1))+'_'+str(round(theta,1))+'_C.png'
    
          used_gates = circuit.count_ops() 
          gates={'u1':0,'u2':0,'u3':0,'cx':0}
          for gate in ('u1','u2','u3','cx'):
              if gate in used_gates:
                  gates[gate]=used_gates[gate]
          depth = circuit.depth()
          filename ='circuit_'+str(round(phi,1))+'_'+str(round(theta,1))+'_C.png'
    
          # String to write in the output file
          string_i=str(phi) + ' ' +str(var_lambda) + ' ' +str(theta)+ ' '+str(C_AB_initial)+' '+str(depth)+' '+str(gates['u1'])+' '+str(gates['u2'])+' '+str(gates['u3'])+' '+str(gates['cx'])+' '+str(tdepth)+' '+str(tgates['u1'])+' '+str(tgates['u2'])+' '+str(tgates['u3'])+' '+str(tgates['cx'])+'\n'
    
          logND.write(string_i)
    
          # --------------------------------
          # QND measurement and tomography of the outgoing state 
          # --------------------------------
    
          string_VP_fnops=''
          string_VP_f1=''
          string_VP_f2=''
        
          qst_VDC = state_tomography_circuits(circuit_t,[qr[0], qr[1]])
          qst_VDC_no_anc = deepcopy(qst_VDC)
    
          ca = ClassicalRegister(1)
          for qst_VDC_circ in qst_VDC:
            qst_VDC_circ.add_register(ca)
            qst_VDC_circ.measure(qr[2],ca[0])
    
          if noise_simulation==1:
            job = execute(qst_VDC, simulator, noise_model=noise_model,coupling_map=coupling_map,basis_gates=basis_gates, shots=shots_tomo_f)
            job_nops = execute(qst_VDC_no_anc, simulator, noise_model=noise_model,coupling_map=coupling_map,basis_gates=basis_gates, shots=shots_tomo_f)
          elif noise_simulation==0:
            job = execute(qst_VDC, backend, shots=shots_tomo_f)
            job_nops = execute(qst_VDC_no_anc, backend, shots=shots_tomo_f)
    
          raw_results=job.result()
          results_nops=job_nops.result()
    
          new_result1 = deepcopy(raw_results)
          new_result2 = deepcopy(raw_results)
    
          # Retrieve outcomes of the ancilla qubit
          for resultidx, _ in enumerate(raw_results.results):
            old_counts = raw_results.get_counts(resultidx)
    
            new_counts1 = {}
            new_counts2 = {}
            new_countsnops= {}
    
            new_result1.results[resultidx].header.creg_sizes = [new_result1.results[resultidx].header.creg_sizes[0]]
            new_result1.results[resultidx].header.clbit_labels = new_result1.results[resultidx].header.clbit_labels[0:-1]
            new_result1.results[resultidx].header.memory_slots = 2
    
            new_result2.results[resultidx].header.creg_sizes = [new_result2.results[resultidx].header.creg_sizes[0]]
            new_result2.results[resultidx].header.clbit_labels = new_result2.results[resultidx].header.clbit_labels[0:-1]
            new_result2.results[resultidx].header.memory_slots = 2
    
            for reg_key in old_counts:
                reg_bits = reg_key.split(' ')
    
                if VPC==1 or VPC==2 or VPC==3:
                  if reg_bits[0]=='0' or reg_bits[0]=='1':
                    new_countsnops[reg_bits[1]]=old_counts[reg_key]
                  if reg_bits[0]=='0':
                    new_counts1[reg_bits[1]]=old_counts[reg_key]
                  if reg_bits[0]=='1':
                    new_counts2[reg_bits[1]]=old_counts[reg_key]
    
    
            new_result1.results[resultidx].data.counts = new_result1.results[resultidx].data.from_dict(new_counts1).to_dict()
            new_result2.results[resultidx].data.counts = new_result2.results[resultidx].data.from_dict(new_counts1).to_dict()
    
          tomo_VDCnops = StateTomographyFitter(results_nops, qst_VDC_no_anc)  
          tomo_VDC1 = StateTomographyFitter(new_result1, qst_VDC_no_anc)
          tomo_VDC2 = StateTomographyFitter(new_result2, qst_VDC_no_anc)
    
    
          rhonops = tomo_VDCnops.fit()
          rho_ABnops=partial_tracey(rhonops,[1,0])
          string_VP_fnops=string_VP_f1+str(phi) + ' ' +str(var_lambda) + ' ' +str(theta)+' '+line(rho_ABnops)+'\n'
    
          if len(new_counts1) != 0 and not(len(new_counts1)==1 and new_counts1.get(list(new_counts1.keys())[0]) < mc) and not(len(new_counts1)==2 and new_counts1.get(list(new_counts1.keys())[0])+new_counts1.get(list(new_counts1.keys())[1]) < mc) and not(len(new_counts1)==3 and new_counts1.get(list(new_counts1.keys())[0])+new_counts1.get(list(new_counts1.keys())[1])+new_counts1.get(list(new_counts1.keys())[2]) < mc) and not(len(new_counts1)==4 and new_counts1.get(list(new_counts1.keys())[0])+new_counts1.get(list(new_counts1.keys())[1])+new_counts1.get(list(new_counts1.keys())[2])+new_counts1.get(list(new_counts1.keys())[3]) < mc):
            rho1 = tomo_VDC1.fit()
            rho_AB1=partial_tracey(rho1,[1,0])
            string_VP_f1=string_VP_f1+str(phi) + ' ' +str(var_lambda) + ' ' +str(theta)+' '+line(rho_AB1)+'\n'
    
          if len(new_counts2) != 0 and not(len(new_counts2)==1 and new_counts2.get(list(new_counts2.keys())[0]) < mc) and not(len(new_counts2)==2 and new_counts2.get(list(new_counts2.keys())[0])+new_counts2.get(list(new_counts2.keys())[1]) < mc) and not(len(new_counts2)==3 and new_counts2.get(list(new_counts2.keys())[0])+new_counts2.get(list(new_counts2.keys())[1])+new_counts2.get(list(new_counts2.keys())[2]) < mc) and not(len(new_counts2)==4 and new_counts2.get(list(new_counts2.keys())[0])+new_counts2.get(list(new_counts2.keys())[1])+new_counts2.get(list(new_counts2.keys())[2])+new_counts2.get(list(new_counts2.keys())[3]) < mc):
            rho2 = tomo_VDC2.fit()
            rho_AB2=partial_tracey(rho2,[1,0])
            string_VP_f2=string_VP_f2+str(phi) + ' ' +str(var_lambda) + ' ' +str(theta)+' '+line(rho_AB2)+'\n'
    
          logfnops.write(string_VP_fnops)
          logf1.write(string_VP_f1)
          logf2.write(string_VP_f2)

logND.close()
logfnops.close()
logf1.close()
logf2.close()
