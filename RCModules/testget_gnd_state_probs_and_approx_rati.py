#testget_gnd_state_probs_and_approx_ratio
from pyaqc.RCModules.aqc_rigetti import get_gnd_state_probs_and_approx_ratio
#ONE AMBULANCE adjacency table
from pyaqc.RCModules.aqc_rigetti import get_gnd_state_probs_and_approx_ratio

from pyaqc.RCModules.PlatformConversion import print_QUBOdetails, CreateAmbulanceAdjacency,distance #   qc_all_env and rigoct38_env
import numpy as np 
import unicodedata as ud
betagreek =  ud.lookup('GREEK SMALL LETTER BETA')
gammagreek =  ud.lookup('GREEK SMALL LETTER gamma')
thetagreek =  ud.lookup('GREEK SMALL LETTER THETA')
################  Set up type of problem to solve ################
Ansatz_type={}
Ansatz_type['n_ambulance'] = 1
Ansatz_type['n_destinations'] = 5
Ansatz_type['n_qubits'] = Ansatz_type['n_destinations']
Ansatz_type['HammingWeightOfConstraint']  =  Ansatz_type['n_qubits'] - 1

gridWidth = Ansatz_type['n_qubits'] 
Ansatz_type['ConstraintMultiplier'] = 14#200             # I have used 14 for 5 qubit with constraints and traditional mixer, 140 for q11, tim 200
# 1500 for q17
# 800 for q8
Ansatz_type['use_XYMixer_constraints'] =  0   # remove the constraints when using an XY mixer
select_qubo_model = 1   # we nearly always use this
################ Create 1) 'Adjacency' table of the problem, and 2) 'Adjacency_constraint' of the problems constraints 3) Adjacency_feasible ################
#'Adjacency' table of the whole problem

Ansatz_type['Adjacency'] = CreateAmbulanceAdjacency(   
    gridWidth,
    n_qubits = Ansatz_type['n_qubits'], 
    ConstraintMultiplier = Ansatz_type['ConstraintMultiplier'], 
    Adddistance = 1,
    remove_constraint= Ansatz_type['use_XYMixer_constraints'] , 
    HammingWeightOfConstraint= Ansatz_type['HammingWeightOfConstraint'] ,
    qubo_model=select_qubo_model)
#Adjacency of just the problems constraints

Ansatz_type['Adjacency_constraint'] = CreateAmbulanceAdjacency(   
    gridWidth,
    n_qubits = Ansatz_type['n_qubits'], 
    ConstraintMultiplier = Ansatz_type['ConstraintMultiplier'], 
    Adddistance = 0,
    remove_constraint= Ansatz_type['use_XYMixer_constraints'] , 
    HammingWeightOfConstraint= Ansatz_type['HammingWeightOfConstraint'] ,
    qubo_model=select_qubo_model)
# Adjacency that defines whether a given state is 'feasible'
Ansatz_type['Adjacency_feasible'] = CreateAmbulanceAdjacency(   
    gridWidth,
    n_qubits = Ansatz_type['n_qubits'], 
    ConstraintMultiplier = Ansatz_type['ConstraintMultiplier'], 
    Adddistance = 0,
    remove_constraint= 0, 
    HammingWeightOfConstraint= Ansatz_type['HammingWeightOfConstraint'] ,
    qubo_model=select_qubo_model)
Adjacency = Ansatz_type['Adjacency']
n_qubits = Ansatz_type['n_qubits']
################ ... and print ################
if select_qubo_model:print('Qubo ', end='')
else: print('Ising ', end='')
print_QUBOdetails(Ansatz_type['Adjacency'] ,Ansatz_type['n_qubits'], 'Ambulance problem')
############################# convert Adjacency to Pauli sum of observeable edges and nodes #############################
from pyaqc.RCModules.aqc_rigetti import create_ring_mixer, Adjacency_qubo_to_Regetti,Adjacency_Ising_to_Regetti,Dicke_state
ListPauli_termsMy,SumPauli_termsMy  = Adjacency_qubo_to_Regetti( Adjacency, n_qubits)

print(SumPauli_termsMy)

from pyaqc.RCModules.aqc_rigetti import get_gnd_state_probs_and_approx_ratio,ansatz_init_XYMixer,ansatz_prog_init,ansatz_EV_init
from pyaqc.RCModules.aqc_rigetti import create_ring_mixer,create_local_ring_mixer,Dicke_state , Dicke_state_local
from pyquil.paulis import *  #  exponential_map, exponentiate_commuting_pauli_sum
from pyquil.api import *
from pyquil.gates import *

def get_mixer(n_qubits,remove_constraint,n_ambulance ):
    """
    Return Rigetti observables that can be exponentiated, and a description 'Ansatz_type'
    type: Pauli Observables
    """
    MyMixerHam = []
    MyMixerHam.append(0)                # Perhaps required by QAOA()
    if not remove_constraint:           # when the extra hard constraints are in place use an X mixer, otherwise an XY mixer
        for n in range(n_qubits):       
            MyMixerHam[0] += 1*sX(n) 
        Ansatz_type['Mixer'] ='TRADITIONAL /  Farhi  mixer'
    else:
        XYmixer = create_ring_mixer(n_qubits)
        MyMixerHam[0] = XYmixer         #XY
        Ansatz_type['Mixer'] ='XY Ring mixer'
        if n_ambulance==2:
            MyMixerHam, prog_init_three_XYMixers = three_XYMixers_Dicke_states(n_qubits, n_destinations) 
            Ansatz_type['prog_init'] = 'Three Dicke states'
            Ansatz_type['Mixer'] ='Three XY Ring mixers' 
    
    if n_ambulance==2:
        labels= ['Destination XYmixer','Start0 XYmixer','Start1 XYmixer']
        for n_mixer, mixer in enumerate(MyMixerHam):
            print(labels[n_mixer])
            print(mixer)
    else: 
        print(Ansatz_type['Mixer'])
        print(MyMixerHam[0])        
        
    return MyMixerHam, Ansatz_type
def get_prog_init(n_qubits,remove_constraint,n_ambulance ):
    """
    Return the Quantum circuit of the initial state (eg Hamamard for Xmixer)
    """
    prog_init = Program()

    if not remove_constraint:
        #A) 
        print('Hadamard initial state.   ', n_qubits, ' = n_qubits')
        for n in range(n_qubits):
            prog_init += H(n)
    else:
        if n_ambulance==2:
            #B)
            #MyMixerHam, prog_init_three_XYMixers = three_XYMixers_Dicke_states(n_qubits, n_destinations)  
            prog_init = prog_init_three_XYMixers       # this is just a reminder that the calculation was made earlier in the 'Mixer used' cell
            print('The initial state is made of three separate Dicke_state s')
        else:
            #C)
            print('Dicke_state initial state',HammingWeightOfConstraint, ' = HW')
            prog_init = Dicke_state(n_qubits, HammingWeightOfConstraint)

    return prog_init
# Build an ansatz based on building blocks that act on just two qubit

def ansatz_init(n_qubits,depth_form=1):
    
    theta_eg=[] 
    last_param_of_qubit = [0 for x in range(n_qubits)] 
    def ansatz_prog(theta, create_theta=0):
        """
        Build an ansatz based on building blocks that act on just two qubit, this structure results in smallish causal cones
        When create_theta=1
            Return  an example of theta for the parameterised circuit, with the last rotation on a qubit being set to pi/2, all others to zero.
        Else
            Return The Rigetti program of the parameterised circuit, 2 qubit blocks each with a CNOT gate followed by RY rotations on each. 
        """
        d=0 # dth element of parameters_ansatz 
         
        prog = Program()
        #first column of RY 
        if 0:
            for q in range(n_qubits):
                
                if create_theta: theta_eg.append(0)
                else:prog += RY(theta[d],q)
                d+=1
          
        def block_2_qubits(q,d,theta):
            'return program of CNOT(q,q+1) RY(theta[d],q) RY(theta[d+1],q+1) '
            prog = Program()
            #CNOT gate after first column
            prog +=CNOT(q,q+1)
            if create_theta: theta_eg.append(0)
            else:prog += RY(theta[d],q)
            last_param_of_qubit[q] =d
            d+=1
            q+=1
            if create_theta: theta_eg.append(0)
            else:prog += RY(theta[d],q)
            last_param_of_qubit[q] =d
            d+=1
            return prog,d
        for D in range (depth_form):
            #First column of 2 qubit blocks on CNOT(q_even,1+q_even)
            for q in range(n_qubits):
                if q%2==0 and q+1 < n_qubits:
                    prog_block,d = block_2_qubits(q,d,theta)
                    prog +=prog_block
            #second column of 2 qubit blocks on CNOT(q_odd,1+q_odd)
            for q in range(n_qubits):
                if q%2==1 and q+1 < n_qubits:
                    prog_block,d = block_2_qubits(q,d,theta)
                    prog +=prog_block
        
        if create_theta:
            for last_d in last_param_of_qubit:
                theta_eg[last_d] = np.pi/2
            return theta_eg
        if not create_theta:return prog
    # Create an intial theta
    ansatz_prog(theta_eg, create_theta=1)
    
    return ansatz_prog, theta_eg
from pyaqc.RCModules.aqc_rigetti import get_gnd_state_probs_and_approx_ratio_VQE_init,min_energy,get_gnd_state_prob

# Test either a QAOA, or a VQE, ansatz approx_ratio calculation
Use_VQE =1
if Use_VQE:
    state_feasible = '11011'
    ansatz_prog, theta_eg = ansatz_init(Ansatz_type['n_qubits'], depth_form=1)
    get_gnd_state_probs_and_approx_ratio = get_gnd_state_probs_and_approx_ratio_VQE_init(ansatz_prog, prt_details=1, state_feasible=state_feasible,**Ansatz_type)
    print(theta_eg,'\n',get_gnd_state_probs_and_approx_ratio(theta_eg))
    if 0:print(ansatz_prog(theta_eg))
    if 1:
    
        E_gnd_state = min_energy( Ansatz_type['Adjacency'],Ansatz_type['n_qubits'])

        Prob = get_gnd_state_prob(theta_eg,ansatz_prog,SumPauli_termsMy,Ansatz_type['n_qubits'], E_gnd_state,prt_details=False)
        print(E_gnd_state, ' = E_gnd_state, with prob of ',Prob )

else:
    MyMixerHam, Ansatz_type = get_mixer(n_qubits,Ansatz_type['remove_constraint'],Ansatz_type['n_ambulance'] )
    prog_init = get_prog_init(n_qubits,Ansatz_type['remove_constraint'],Ansatz_type['n_ambulance'] )

    # code plan: find bg for a q5 1 A problem
    q5p1TradXmHqq = [6.898, -0.022] # nm -252.8 #11.1% Energy_max =-20 , Energy_min = -40 when constraint are complied with and have zero energy contribution to the cost function
    opt_betagamma = q5p1TradXmHqq = [6.898, -0.022] # nm -252.8 #11.1%
    p = len(opt_betagamma)//2
    state_feasible = '11011'

    ansatz_prog = ansatz_prog_init(prog_init,   ListPauli_termsMy, n_qubits, Ansatz_type['n_destinations'],   p,MyMixerHam)
if 0:
    prob_gnd_state, approx_ratio,prob_feasible_state = get_gnd_state_probs_and_approx_ratio(opt_betagamma,ansatz_prog,SumPauli_termsMy,n_qubits,Adjacency_constraint=Ansatz_type['Adjacency_constraint'], state_feasible=state_feasible, prt_details=1)
    print(prob_gnd_state, approx_ratio,prob_feasible_state)
from pyaqc.RCModules.aqc_rigetti import  get_gnd_state_prob, min_energy
