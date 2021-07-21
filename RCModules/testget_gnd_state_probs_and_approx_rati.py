#testget_gnd_state_probs_and_approx_ratio
from pyaqc.RCModules.aqc_rigetti import get_gnd_state_probs_and_approx_ratio
#ONE AMBULANCE adjacency table
from pyaqc.RCModules.PlatformConversion import print_QUBOdetails, CreateAmbulanceAdjacency,distance #   qc_all_env and rigoct38_env
import numpy as np 
import unicodedata as ud
betagreek =  ud.lookup('GREEK SMALL LETTER BETA')
gammagreek =  ud.lookup('GREEK SMALL LETTER gamma')

################  Set up type of problem to solve ################
n_destinations = 5
n_qubits = n_destinations
HammingWeightOfConstraint =  n_qubits - 1

gridWidth = n_qubits
ConstraintMultiplier = 14              # I have used 14 for 5 qubit with constraints and traditional mixer, 140 for q11, tim 200
remove_constraint=  0   # remove the hard constraints when using an XY mixer
select_qubo_model = 1   # we nearly always use this
################ Create 1) adjacency table of the problem, and 2) of the hard constraints ################
Adjacency_constraint = CreateAmbulanceAdjacency(   
    gridWidth,
    n_qubits = n_qubits, 
    ConstraintMultiplier = ConstraintMultiplier, 
    Adddistance = 0,
    remove_constraint= remove_constraint, 
    HammingWeightOfConstraint= HammingWeightOfConstraint,
    qubo_model=select_qubo_model)
Adjacency = CreateAmbulanceAdjacency(   
    gridWidth,
    n_qubits = n_qubits, 
    ConstraintMultiplier = ConstraintMultiplier, 
    Adddistance = 1,
    remove_constraint= remove_constraint, 
    HammingWeightOfConstraint= HammingWeightOfConstraint,
    qubo_model=select_qubo_model)
################ ... and print ################
if select_qubo_model:print('Qubo ', end='')
else: print('Ising ', end='')
print_QUBOdetails(Adjacency,n_qubits, 'Ambulance problem')
n_ambulance=1
Ansatz_type = {'n_qubits':n_qubits, 'p':0,'Mixer':0,'prog_init':0, 'ConstraintMultiplier':0, 'remove_constraint':0, 'Optimizer_maxiter': 0}
Ansatz_type['ConstraintMultiplier'] = ConstraintMultiplier
Ansatz_type['remove_constraint'] = remove_constraint
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

MyMixerHam, Ansatz_type = get_mixer(n_qubits,remove_constraint,n_ambulance )
prog_init = get_prog_init(n_qubits,remove_constraint,n_ambulance )

# code plan: find bg for a q5 1 A problem
q5p1TradXmHqq = [6.898, -0.022] # nm -252.8 #11.1%
opt_betagamma = q5p1TradXmHqq = [6.898, -0.022] # nm -252.8 #11.1%
p = len(opt_betagamma)//2
state_feasible = '11011'
ansatz_prog = ansatz_prog_init(prog_init,   ListPauli_termsMy, n_qubits, n_destinations,   p,MyMixerHam)
prob_gnd_state, approx_ratio,prob_feasible_state = get_gnd_state_probs_and_approx_ratio(opt_betagamma,ansatz_prog,SumPauli_termsMy,n_qubits,Adjacency_constraint=Adjacency_constraint, state_feasible=state_feasible, prt_details=1)
print(prob_gnd_state, approx_ratio,prob_feasible_state)