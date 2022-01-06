# aqc_rigetti.py
# aqc functions that rely on the rigetti pyquil funtions
from pyquil import *
from pyquil.gates import *
import numpy as np
from pyquil.api import WavefunctionSimulator
from pyquil.paulis import ID, sZ, sX, sY , exponential_map, exponentiate_commuting_pauli_sum
from pyaqc.RCModules.QAOARC import decimal_state_to_binary_reversed, Energy_of_binary_state
print('rig_dec_2021_env')
"""
from 2021_env and rig_dec_2021_env
# Function alternative to Rigetti's QAOA()
def ansatz_init_XYMixer( BetaGamma, ListPauli_termsMy, n_qubits, n_destinations,p, XY_or_Xmixer_list):
def ansatz_EV_init(prog_init,  pauli_cost_terms_list, n_qubits, n_destinations,   p,MyMixerHam):
def ansatz_prog_init(prog_init,   pauli_cost_terms_list, n_qubits, n_destinations,   p,MyMixerHam):
def get_gnd_state_probs_and_approx_ratio(opt_betagamma,ansatz_prog,SumPauli_termsMy,n_qubits, Adjacency_constraint=None, state_feasible=None): BEING DEPRECATED
def min_energy( Adjacency, n_qubits, state_feasible=None):

def Dicke_state_local( 
p: 'Type:Program()',
qudit_start:'FIRST qubit in ring mixer(Type: int)'=0, 
qudit_end: 'LAST qubit in ring mixer (Type: int)'=2,  HWeight: 'int'=1):
def Dicke_state(n_qubits: int , HWeight: int):
def fixed_hw_states(init_prog0, qudit_start, qudit_end, HW, init_state='000111')
def three_XYMixers_Dicke_states(n_qubits, n_destinations,init_state=None)

def create_local_ring_mixer( qudit_start:'FIRST qubit in ring mixer'=0, qudit_end:'LAST qubit in ring mixer'=2):
def create_ring_mixer(n_qubits):

def Adjacency_Ising_to_Regetti( Adjacency: 'Adjacency table in form: qubo[(r,c)]', n_qubits):
def Adjacency_qubo_to_Regetti( Adjacency: 'Adjacency  table in form: qubo[(r,c)]', n_qubits):

def get_prog_init( **Ansatz_type )



def get_mixer( print_details=True,**Ansatz_type )
def get_mixer(n_qubits,use_XY_mixer,n_ambulance )
def get_prog_init(n_qubits,use_XY_mixer,n_ambulance )
get_QAOA_circuit(n_qubits, n_destinations,n_ambulance,p,use_XYMixer_constraints,HammingWeightOfConstraint,Adjacency,Adjacency_constraint,Adjacency_feasible, prt_details=1, state_feasible=None,Ansatz_type=None)
def get_QAOA_circuit(p, prt_details=1, state_feasible=None,**Ansatz_type)
def get_gnd_state_probs_and_approx_ratio_simple_init(p, prt_details=1, state_feasible=None,**Ansatz_type)

get_gnd_state_probs_and_approx_ratio_simple_init(n_qubits, n_destinations,p,Adjacency_constraint,Adjacency_feasible, prt_details=1, state_feasible=None):

def get_approx_ratio_init(ansatz_prog,prog_init,ListPauli_termsMy,SumPauli_termsMy,n_qubits,p,Adjacency_constraint,Adjacency_feasible,  prt_details=1,state_feasible=None )

def get_gnd_state_probs_and_approx_ratio(opt_betagamma,  prt_details=1) 
"""
def ansatz_init_XYMixer( BetaGamma, pauli_cost_terms_list, n_qubits, n_destinations,p, XY_or_Xmixer_list):
    """
    param  initial_angle  eg # uses betagammon = [beta0,beta1, beta2, ,gamma0, gamma1,gamma2]. Note there is ONE beta  for each gamma. 
        Type: list float
    'pauli_cost_terms_list'     A list of the terms of each edge in the Adjacency. This must not be a sum of all the terms
        Type:list of pauli observables
    param n_destinations  the number of ambulance destinations 
        Type:int
    param XY_or_Xmixer_list A Pauli list for a ringmixer or an Xmixer
    Returns prog , p repeats of ,  exponentials of the cost paulisums in  XYmixer_list, then in SumPauli_termsMy (cost terms).
        Type: Program() 
    """
    #Check the correct number of angles have been supplied (1 beta and 1 gamma per p, )
    beta_list = BetaGamma[0:p  ]  
    if 2* p != len( BetaGamma):
        raise AssertionError( 'Number of angles supplied:',len( BetaGamma), 'do not match the 2*p', 2* p)
    prog = Program() 
    gamma_list = BetaGamma[p:  ]
    if 0: print(gamma_list, ' = gamma_list')
    SumPauli_all_terms = 0
    for pauli_sum in pauli_cost_terms_list:     
        SumPauli_all_terms += pauli_sum
    for pn in range(p):        
        for term in SumPauli_all_terms:                       #SumPauli_termsMy
            prog +=exponential_map(term)(gamma_list[pn])     #gamma_i rotations one for each edge in the ring mixer per p
        for n, XYmixer in enumerate( XY_or_Xmixer_list) :
             for termXY in XYmixer:
                    prog +=exponential_map(termXY)(beta_list[pn])
        
    return prog
def ansatz_EV_init(prog_init,  pauli_cost_terms_list, n_qubits, n_destinations,   p, MyMixerHam):
    """
    Returns (type: python function) that takes 'initial_angle' as a parameter and returns the expectation of the ansatz at that angle
    param 
        'pauli_cost_terms_list'     A list of the terms of each edge in the Adjacency. This must not be a sum of all the terms
        'prog_init' The initial program, eg a series of Hadamard gates, or a fixed Hamming weight program, or a Dickie state.
    """
    SumPauli_all_terms = 0
    for pauli_sum in pauli_cost_terms_list:     
        SumPauli_all_terms += pauli_sum
    def ansatz_EV( initial_angle ):
        #prog = prog_init + ansatz_init_XYMixer(initial_angle,pauli_cost_terms_list, n_qubits, n_destinations, p,MyMixerHam )
        prog = prog_init + ansatz_init_XYMixer(initial_angle,SumPauli_all_terms, n_qubits, n_destinations, p,MyMixerHam )
        EV = WavefunctionSimulator().expectation(prep_prog=prog, pauli_terms=SumPauli_all_terms)
        #EV = WavefunctionSimulator().expectation(prep_prog=Known_state_sol, pauli_terms=pauli_cost_terms)
        #print(prog)
        return np.real(EV  )        #convert the potentially complex number to its real component

    return ansatz_EV
def ansatz_prog_init(prog_init, pauli_cost_terms_list, n_qubits, n_destinations,   p,MyMixerHam):
    """
    Returns the name of a function ansatz_prog. When called with an angle, for example; ansatz_prog( initial_angle ), 
    type: Rigetti Program
    It  represents the QAOA Ansatz ( cost Phase, and mixer), excluding the initial state.
    param MyMixerHam is a list of Pauli operators, which can be X or XY operators
    """
    def ansatz_prog( initial_angle ):
        return prog_init + ansatz_init_XYMixer(initial_angle,pauli_cost_terms_list, n_qubits, n_destinations, p ,MyMixerHam)

    return ansatz_prog

def get_gnd_state_prob(opt_Angle,ansatz_prog,SumPauli_termsMy,n_qubits, E_gnd_state,prt_details=False):
    prob_gnd_state = 0
    # EV of state '00...00'
    
    solution_ansatz = ansatz_prog(opt_Angle)        #   dict {state_i:prob_i}
    
    probsAbsolute  = WavefunctionSimulator().wavefunction(solution_ansatz).get_outcome_probs()
    opt_probs_absolute = [x for x in probsAbsolute.values()]
    
    for n in range(2**n_qubits):        
        if opt_probs_absolute[n] > 0.00001 or 1:    #|qn...q0>
            Expectation_of_state = Energy_of_binary_state( n, SumPauli_termsMy,n_qubits)    
            if Expectation_of_state <= E_gnd_state:   
                prob_gnd_state += opt_probs_absolute[n]
                EV_min = Expectation_of_state
    return prob_gnd_state

def get_gnd_state_probs_and_approx_ratio(opt_betagamma,ansatz_prog,SumPauli_termsMy,n_qubits,Adjacency_constraint=None, 
state_feasible=None, prt_details=False):
    """
    BEING DEPRECATED 
        instead use get_gnd_state_prob(theta,ansatz_prog,SumPauli_termsMy,Ansatz_type['n_qubits'], E_gnd_state,prt_details=False)
    return  A) the cumulative probability of the lowest energy state sampled, 
            B) 'approx_ratio' from the feasible states, (defined as states that minimize 'Adjacency_constraint') the average energy calculated without any constraint energy / ground state energy.
                    Or -1 if 'Adjacency_constraint' is not supplied
                    As defined by Wang in arXiv:1904.09314v2 [quant-ph] 21 May 2020, page 3 col 2
                    When the ground state without constraints is positive the reciprocal of 'approx_ratio' is calculated
            C) The probability that one of the suggested states is a feasible state(ie complies with all constraint )
    type: float,float

    param
        ansatz_prog
            type:function that returns a Rigetti.Program()
    kwargs
        
        'Adjacency_constraint'  The adjacency table which contain only the constraint energies , and no Soft constraints (eg distance)
            type: dict example {(0, 0): 5,  (1, 1): -1,  (2, 2): -1}
        This is required to calculate the approximation ratio
        'state_feasible' a string of a state that complies with all the constraints in Adjacency_constraint. For example;  '11011'.capitalize. 
            type: string of a binary state
        This speeds up the approximation ratio calculation
        
    """
    prob_gnd_state = 0
    prob_feasible_state = 0
    EV_feasible = 0 
    approx_ratio = -1
    seen_first_compliant_state = 0
    
    EV_min  = Expectation_of_state = Energy_of_binary_state( 0, SumPauli_termsMy,n_qubits)
    
    # EV of state '00...00'
    #EV_min = Energy_of_binary_state( state_str, SumPauli_termsMy,n_qubits)
    solution_ansatz = ansatz_prog(opt_betagamma)        #   dict {state_i:prob_i}
    n_feasible=0
    
    if not Adjacency_constraint==None: 
        ' get  Minimum energy of the Adjacency by considering all 2^n_qubits states. OR if state_feasible is supplied the adjacenys energy in that state  type:float'
        min_energy_constraint = min_energy( Adjacency_constraint, n_qubits, state_feasible=state_feasible.replace(' ',''))
        ListPauli_terms_constr,SumPauli_termsMy_constr  = Adjacency_qubo_to_Regetti( Adjacency_constraint, n_qubits)
        
    probsAbsolute  = WavefunctionSimulator().wavefunction(solution_ansatz).get_outcome_probs()
    opt_probs_absolute = [x for x in probsAbsolute.values()]
    for n in range(2**n_qubits):
        
        if opt_probs_absolute[n] > 0.00001 or 1:    #|qn...q0>
            Expectation_of_state = Energy_of_binary_state( n, SumPauli_termsMy,n_qubits)    
            # calculation of energy, using n as the state, NOT recalculating the state by building its quantum circuit eg Program(X(0),I(1))
            #Expectation_of_state = Energy_of_binary_state( n, SumPauli_termsMy,n_qubits)

            if Expectation_of_state <= EV_min:   
                
                psi_opt, state_string = decimal_state_to_binary_reversed(n,n_qubits)
                if Expectation_of_state < EV_min:
                    prob_gnd_state = 0
                prob_gnd_state += opt_probs_absolute[n]
                EV_min = Expectation_of_state
            if not Adjacency_constraint==None:
                # If the Adjacency_constraint has an energy that is the minimium possible, the constraints have been met..
                #...then the energy calculated without any constraint energy is averaged
                if Energy_of_binary_state( n, SumPauli_termsMy_constr, n_qubits) == min_energy_constraint:
                    if seen_first_compliant_state==0:
                        Energy_max = EV_min - min_energy_constraint 
                        #print(Energy_max , ' = Energy_max, first estimate')

                    seen_first_compliant_state = 1
                    EV_feasible += opt_probs_absolute[n] * (Expectation_of_state - min_energy_constraint)
                    prob_feasible_state +=opt_probs_absolute[n]
                    n_feasible +=1
                    #print(opt_probs_absolute[n] , (Expectation_of_state - min_energy_constraint),EV_feasible)
                    if (EV_feasible/prob_feasible_state < 2 or (Expectation_of_state - min_energy_constraint) <2) and 0:
                        print(EV_feasible/prob_feasible_state, 'cum wt avg')
                        print(n,str(bin(n)[2:]) ,'is feasible',min_energy_constraint,(Expectation_of_state - min_energy_constraint))
                    # find max energy when constraint energy is at a minimum
                    Energy_max = max(Energy_max,Expectation_of_state-min_energy_constraint)
                    

    if not prob_feasible_state==0:
        if EV_min-min_energy_constraint-Energy_max ==0: approx_ratio = 0
        else:approx_ratio = (EV_feasible/prob_feasible_state-Energy_max)/(EV_min-min_energy_constraint-Energy_max) 
        #print(EV_feasible/prob_feasible_state-Energy_max,'=(EV_feasible/prob_feasible_state-Energy_max)',EV_min-min_energy_constraint-Energy_max,'= EV_min-min_energy_constraint-Energy_max')
    if prt_details:
        print(EV_feasible/prob_feasible_state, ' = EV_feasible/prob_feasible_state,', EV_min, ' = EV_min,', min_energy_constraint, ' = min_energy_constraint.',prob_feasible_state, ' =prob_feasible_state,',n_feasible, ' = n_feasible,')
        print(Energy_max , ' = Energy_max, with constraint Energy=0')
    return prob_gnd_state, approx_ratio,prob_feasible_state

def min_energy( Adjacency, n_qubits, state_feasible=None):
    """
    Return  Minimum energy of the Adjacency by considering all 2^n_qubits states. OR if state_feasible is supplied the adjacency's energy in that state (less run time).
    type:float
    """

    ListPauli_terms_constr,SumPauli_termsMy_constr  = Adjacency_qubo_to_Regetti( Adjacency, n_qubits)
    min_energy_constraint = Energy_of_binary_state( 0, SumPauli_termsMy_constr,n_qubits)
    if state_feasible==None:
        for n in range(2**n_qubits):
            min_energy_constraint = min(min_energy_constraint,Energy_of_binary_state( n, SumPauli_termsMy_constr,n_qubits))
    else:
        if not len(state_feasible)== n_qubits:
            raise AssertionError('In get_gnd_state_probs_and_approx_ratio() call to min_energy(), len(state_feasible)',len(state_feasible), 'does not equal', n_qubits, 'n_qubits')
        feasible_int = int(state_feasible,2)        #convert binary string in base 2 to integer base 10

        min_energy_constraint = min(min_energy_constraint,Energy_of_binary_state( feasible_int, SumPauli_termsMy_constr,n_qubits))
    return min_energy_constraint

def Dicke_state_local( 
p: 'Type:Program()',
qudit_start:'FIRST qubit in ring mixer(Type: int)'=0, 
qudit_end: 'LAST qubit in ring mixer (Type: int)'=2,  HWeight: 'int'=1):
    """
    Returns type: Rigetti program. 
    
    Returns a program that is an even superposition of all states, in a subsystem of from qudit_start to qudit_end ,   of Hamming Weight = HWeight
    
    As suggested by H. Wang, S. Ashhab, and F. Nori, Physical Review A 79, 042335 (2009).
    """
    
    """ A chapter is defined as rotations of the same qubit ('rotated_qubit')  and NOT gates all targeting  another qubit ('noted_qubit')
The rotations are controlled by
    noted_qubit and rotated_qubit+1
Number of rotations per chapter
    noted_qubit-2       r in range(1,noted_qubit-1):
Rotation angle for rth rotation in that chapter
    2*np.arccos(np.sqrt(r/(noted_qubit + 1)

"""
    # check Program argument has qubits that are to be modified
    if min( p.get_qubits() ) > qudit_start or ( max( p.get_qubits() ) < qudit_end):
        raise AssertionError('Qubits specified are not in the Program p')
    for n in range(qudit_start,qudit_end+1):
        ##Initial state
        if n < qudit_end +1  - HWeight: 
            p += I(n)
        else:
            p += X(n)
        #   'noted_qubit' refers to q(noted_qubit + qudit_start ) and 
        #   'rotated_qubit' refers to q( rotated_qubit+ qudit_start )
    for noted_qubit in range(qudit_end - qudit_start, 0,-1):
        #Max Number of rotations in first/any chapter
        #print('new chapter')
        n_rotations = qudit_end - qudit_start
        if n_rotations < 0: n_rotations = 0
        test=1        
        for r in range(1,n_rotations+1):
            rotated_qubit = noted_qubit - r
            if rotated_qubit >= 0:
                #print('rotated_qubit=', rotated_qubit, 'noted_qubit=',noted_qubit, 'rth rotation= ',r)
                #p +=CNOT(rotated_qubit,noted_qubit )        #CNOT(cntrl, target)
                p +=CNOT(rotated_qubit+ qudit_start,noted_qubit+ qudit_start )        #CNOT(cntrl, target)
                if noted_qubit == rotated_qubit+1:
                    #print('single control')
                    #p += RY(2*np.arccos(np.sqrt(r/(noted_qubit + 1)) ), rotated_qubit).controlled(noted_qubit)
                    p += RY(2*np.arccos(np.sqrt(r/(noted_qubit + 1)) ), rotated_qubit + qudit_start).controlled(noted_qubit+ qudit_start)
                else:
                    #p += RY(2*np.arccos(np.sqrt(r/(noted_qubit + 1)) ), rotated_qubit).controlled(noted_qubit).controlled(rotated_qubit+1)
                    p += RY(2*np.arccos(np.sqrt(r/(noted_qubit + 1)) ), rotated_qubit+ qudit_start).controlled(noted_qubit+ qudit_start).controlled(rotated_qubit+1+ qudit_start)
                #p +=CNOT(rotated_qubit,noted_qubit)
                p +=CNOT(rotated_qubit+ qudit_start,noted_qubit+ qudit_start)
    return p

"""
Each Chapter selects a different qubit to be the target of the CNOT gates, and has one less rotations that the previous chapter
                Chapter1 (qudit_end  == noted_qubit)        n_rotations = qudit_end - qudit_start                                                               
0
qudit_start                                                                                                         .    C_RY( )   .          
            rotated_qubit                                                                  .    C_RY( )   .                 .        
            in_between_qubit                                .       C_RY( )   .                     .    
            in_between_qubit    .       C_RY( )    .                . 
qudit_end   noted_qubit        (Cnot)    .     (Cnot)      (Cnot)   .       (Cnot)       (Cnot)     .   (Cnot)    (Cnot)     .   (Cnot)     
qubit_n

0            Chapter2 (qudit_end - 1  == noted_qubit)       n_rotations -= 1
qudit_start rotated_qubit                                                                   .    C_RY( )   .
            in_between_qubit                                .       C_RY( )   .                     .    
            in_between_qubit    .       C_RY( )    .                . 
            noted_qubit         (Cnot)    .     (Cnot)      (Cnot)  .       (Cnot)       (Cnot)     .   (Cnot)
qudit_end
qubit_n
"""
def Dicke_state(n_qubits: int , HWeight: int):
    """
    Returns type: Rigetti program. 
    
    Returns a program that is an even superposition of all states, in a system of n_qubits,   of Hamming Weight = HWeight
    
    As suggested by H. Wang, S. Ashhab, and F. Nori, Physical Review A 79, 042335 (2009).
    """
    p = Program()
    """ A chapter is defined as rotations of the same qubit ('rotated_qubit')  and NOT gates targeting the same qubit ('noted_qubit')
The rotations are controlled by
    noted_qubit and rotated_qubit+1
Number of rotations per chapter
    noted_qubit-2       r in range(1,noted_qubit-1):
Rotation angle for rth rotation in that chapter
    2*np.arccos(np.sqrt(r/(noted_qubit + 1)

"""

    for n in range(n_qubits):
        ##Initial state
        if n < n_qubits - HWeight: 
            p += I(n)
        else:
            p += X(n)
        #When q2 is high, then for 40% of the time not q2 (to low), and q1 to hi
    for noted_qubit in range(n_qubits-1,0,-1) :
        #Max Number of rotations in first/any chapter
        #print('new chapter')
        n_rotations = n_qubits-1
        if n_rotations < 0: n_rotations=0
        for r in range(1,n_rotations+1):
            rotated_qubit = noted_qubit - r
            if rotated_qubit >=0:
                #print('rotated_qubit=', rotated_qubit, 'noted_qubit=',noted_qubit, 'rth rotation= ',r)
                p +=CNOT(rotated_qubit,noted_qubit )        #CNOT(cntrl, target)
                if noted_qubit == rotated_qubit+1:
                    #print('single control')
                    p += RY(2*np.arccos(np.sqrt(r/(noted_qubit + 1)) ), rotated_qubit).controlled(noted_qubit)
                else:
                    p += RY(2*np.arccos(np.sqrt(r/(noted_qubit + 1)) ), rotated_qubit).controlled(noted_qubit).controlled(rotated_qubit+1)
                p +=CNOT(rotated_qubit,noted_qubit)
    return p
def fixed_hw_states(init_prog0, qudit_start, qudit_end, HW, init_state='000111'):
    """
    return 'init_prog0' with a fixed hamming weight('HW') program state between qubits; qudit_start, qudit_end
    """
    num_q = qudit_start
    if len(init_state) != qudit_end - qudit_start +1:
        print('len(init_state) != qudit_end - qudit_start')
    for pos, bit in enumerate(init_state):
        if int(bit)==1 :
            #print(bit)
            init_prog0 += X(num_q+pos)    

    return init_prog0

def three_XYMixers_Dicke_states(n_qubits, n_destinations,init_state=None):
    """
    Return (type:paulisum) observeable that when exponentiated represents the Phase separating Hamiltonian of an Ansatz

    para 'init_state', when set there is a single hamming weight init_prog in state=(init_state) rather that a Dickie even superpostion of all similar HW states.
    type:string of a state. example '1110111'
        """
    init_prog0 = Program()
    for num_q in range(n_qubits):
        init_prog0 += I(num_q)
    XYmixer = 0
    # DESTINATIONS Make XYmixer0 and Dicke state for subcircuit:
    XYmixer_list = []
    qudit_start = 0
    qudit_end = qudit_start + 2 * n_destinations-1
    HW = n_destinations
    if init_state== None:
        init_prog1 = Dicke_state_local(init_prog0, qudit_start, qudit_end, HW)
    else:
        init_prog1 = fixed_hw_states(init_prog0, qudit_start, qudit_end, HW, init_state='11000011')
        
    # XY to constrain the number of destinations service to n_destinations
    # Create xyMixer0, which will be rotated/exponentiated by beta0 and will have a HW dickie state of n_destinations
    XYmixer_list.append( create_local_ring_mixer(qudit_start, qudit_end))

    ## START_A0 : Make XYmixer1 and Dicke state for subcircuit:

    qudit_start = qudit_end + 1
    qudit_end = qudit_start +  n_destinations-1
    HW=1
    #init_prog2 = Dicke_state_local(init_prog1, qudit_start, qudit_end, HW)
    if init_state== None:
        init_prog2 = Dicke_state_local(init_prog1, qudit_start, qudit_end, HW)
    else:
        init_prog2 = fixed_hw_states(init_prog1, qudit_start, qudit_end, HW, init_state='1000')
        
    # Create xyMixer1, which will be rotated/exponentiated by beta1 and will have a HW dickie state of 1
    XYmixer_list.append(create_local_ring_mixer(qudit_start, qudit_end))

    qudit_start = qudit_end + 1
    qudit_end = qudit_start +  n_destinations-1
    HW = 1
    #init_prog = Dicke_state_local(init_prog2, qudit_start, qudit_end, HW)
    if init_state== None:
        init_prog = Dicke_state_local(init_prog2, qudit_start, qudit_end, HW)
    else:
        init_prog = fixed_hw_states(init_prog2, qudit_start, qudit_end, HW, init_state='0001')
        
    ##  START_A1 :Make XYmixer2 and Dicke state for subcircuit:

    XYmixer_list.append(create_local_ring_mixer(qudit_start, qudit_end))
    
    return XYmixer_list, init_prog

def create_local_ring_mixer( qudit_start:'FIRST qubit in ring mixer'=0, qudit_end:'LAST qubit in ring mixer'=2):
    """
    returns (type: Pauli_sum) a sum of observables that correspond to a ring_mixer with qubits from qudit_start to qudit_end
    """
    XYmixer =0
    if 1:
        for a in range(qudit_start, qudit_end):
            XYmixer += 0.5*sX(a)*sX(a+1) + 0.5*sY(a)*sY(a+1)
    if 1:   # 1 for ring mixer, 0 for 'Line'
        XYmixer += 0.5*sX(qudit_start)*sX(qudit_end) + 0.5*sY(qudit_start)*sY(qudit_end)      
    
    return XYmixer
def create_ring_mixer(n_qubits):
    """
    returns (type: Pauli_sum) a sum of observables that correspond to a ring_mixer with n_qubits
    """
    
    return create_local_ring_mixer( qudit_start=0, qudit_end=n_qubits-1)
# Open each line in the file and format the line as a list, then add that list to the DataOut list [[line1],  [line2]]
def Adjacency_qubo_to_Regetti( Adjacency:"Adjacency table in form: qubo[(r,c)]", n_qubits):
    """
    Returns: ListPauli_termsMy (type:lists[pauli_sum]),and SumPauli_termsMy (type: pauli_sum)  that are required by the Rigetti function;

    expectation = WavefunctionSimulator().expectation(prep_prog=prog_init, pauli_terms=ListPauli_termsMy)
        or vqe_run( ListPauli_termsMy)
    
    In a qubo, as opposed to an Ising, The expectation of a node or edge of weighting W are:
        Relation        Output/expectation      Qubit condition                 Observable formula
        Node            0                       qubit == 0
        Node            W                       qubit == 1                      W * 0.5 * (1-sZ(r)
        Edge            0                       at least one qubit in the edge is zero
        Edge            W                       both qubits == 1                W * ( 0.25 * sZ(r) * sZ(c) + 0.25 - 0.25*sZ(r) - 0.25*sZ(c) )
    
    # When Program(X(0)) then Expecation of sZ(0) is -1
    # When Program(I(0)) then Expecation of sZ(0) is +1    
    """
    ListPauli_termsMy = []
    SumPauli_termsMy = 0
    EdgeEnergy = 0
    for (r,c),val in Adjacency.items():
        #Adjacency[(r,c)] the order of -1 and sZ(r) is key for later classical calculation
            if c == r:
                SumPauli_termsMy +=Adjacency[(r,c)] * 0.5 * (1-sZ(r))
                ListPauli_termsMy.append(Adjacency[(r,c)] * 0.5 * (1-sZ(r)))  
            if c > r :                                  
                EdgeEnergy = ( 0.5 * sZ(r) * sZ(c) + 0.5 )       # was 0.5 * sZ(r) * sZ(c) - 0.5 
                EdgeEnergy = ( 0.25 * sZ(r) * sZ(c) + 0.25 - 0.25*sZ(r) - 0.25*sZ(c) ) #Has expectation 1 when program(X(0),X(1)) otherwise 0
                ListPauli_termsMy.append(Adjacency[(r,c)]  * EdgeEnergy) 
                SumPauli_termsMy +=Adjacency[(r,c)]  * EdgeEnergy         
    return ListPauli_termsMy,SumPauli_termsMy 

 #   from pyquil.paulis import *
def Adjacency_Ising_to_Regetti( Adjacency: 'Adjacency table in form: qubo[(r,c)]', n_qubits):
    """
    Returns: ListPauli_termsMy (type:lists[Paulisum]) that are required by the Rigetti function;

    expectation = WavefunctionSimulator().expectation(prep_prog=prog_init, pauli_terms=ListPauli_termsMy)

    In a Ising , as opposed to an qubo, The expectation of a node or edge of weighting W are:
        Relation        Output/expectation      Qubit condition
        Node            -W                      qubit == 0
        Node            W                       qubit == 1
        Edge            -W                      both qubits are the DIFFERENT
        Edge            W                       both qubits are the SAME
    """
    # When Program(X(0)) then Expecation of sZ(0) is -1
    # When Program(I(0)) then Expecation of sZ(0) is +1   
    ListPauli_termsMy = []
    SumPauli_termsMy=0
    for (r,c),val in Adjacency.items():
            if c == r:
                ListPauli_termsMy.append(Adjacency[(r,c)] * -1*sZ(r) + 0)  
                SumPauli_termsMy += Adjacency[(r,c)] * -1*sZ(r) + 0
            if c > r :                                  
                ListPauli_termsMy.append(Adjacency[(r,c)] *sZ(r) *sZ(c) +0 )
                SumPauli_termsMy += Adjacency[(r,c)]  *sZ(r) *sZ(c) + 0 
    return ListPauli_termsMy,SumPauli_termsMy
def get_mixer(n_qubits,use_XY_mixer,n_ambulance,Ansatz_type,n_destinations, print_details=True ):
    """
    Return Rigetti observables that can be exponentiated, and a description 'Ansatz_type' (type:dict)
    type: Pauli Observables, for example sX(n) can be exponentiated e^{-beta sum_0^n(sX(n))} which causes a X axis rotation of beta for all n qubits.
    if use_XY_mixer==1 then the mixer created is an XY ring mixer otherwise it is an X mixer 
    if n_ambulance ==2 and use_XY_mixer==1, then 3 XY ring mixers are created otherwise it is an X mixer
    if,
        mixer_{j,k} = sX(j)sX(k) + sY(j)sY(k)
    we can define a ring mixer length n as 
        ring_mixer = mixer_{0,n} + sum_{j=0}^{n-1}(mixer_{j,j+1)
    and create the quantum circuit for an XY ring mixer
        e^{-gamma ring_mixer}
    """
    MyMixerHam = []
    MyMixerHam.append(0)                # Perhaps required by QAOA()
    if  use_XY_mixer:           
        XYmixer = create_ring_mixer(n_qubits)
        MyMixerHam[0] = XYmixer         #XY
        Ansatz_type['Mixer'] ='XY Ring mixer'
        if n_ambulance==2:
            MyMixerHam, prog_init_three_XYMixers = three_XYMixers_Dicke_states(n_qubits, n_destinations) 
            Ansatz_type['prog_init'] = 'Three Dicke states'
            Ansatz_type['Mixer'] ='Three XY Ring mixers' 
            labels= ['Destination XYmixer','Start0 XYmixer','Start1 XYmixer']
            for n_mixer, mixer in enumerate(MyMixerHam):
                if print_details: 
                    print(labels[n_mixer])
                    print(mixer)
    else:
        for n in range(n_qubits):       
            MyMixerHam[0] += 1*sX(n) 
        Ansatz_type['Mixer'] ='TRADITIONAL /  Farhi  mixer' 
        if print_details: 
            print(Ansatz_type['Mixer'])
            print(MyMixerHam[0])        
        
    return MyMixerHam#, Ansatz_type
def get_mixer( print_details=True,**Ansatz_type ):
    """
    Return Rigetti observables that can be exponentiated, and a description 'Ansatz_type' (type:dict)
    type: Pauli Observables, for example sX(n) can be exponentiated e^{-beta sum_0^n(sX(n))} which causes a X axis rotation of beta for all n qubits.
    if use_XY_mixer==1 then the mixer created is an XY ring mixer otherwise it is an X mixer 
    if n_ambulance ==2 and use_XY_mixer==1, then 3 XY ring mixers are created otherwise it is an X mixer
    if,
        mixer_{j,k} = sX(j)sX(k) + sY(j)sY(k)
    we can define a ring mixer length n as 
        ring_mixer = mixer_{0,n} + sum_{j=0}^{n-1}(mixer_{j,j+1)
    and create the quantum circuit for an XY ring mixer
        e^{-gamma ring_mixer}
    """
    #n_qubits,use_XY_mixer,n_ambulance,Ansatz_type,n_destinations,
    HW= Ansatz_type['HammingWeightOfConstraint']
    n_q = Ansatz_type['n_qubits']
    MyMixerHam = []
    MyMixerHam.append(0)                # Perhaps required by QAOA()
    if  Ansatz_type['use_XYMixer_constraints']:           
        XYmixer = create_ring_mixer(n_q)
        MyMixerHam[0] = XYmixer         #XY
        Ansatz_type['Mixer'] ='XY Ring mixer'
        if Ansatz_type['n_ambulance']==2:
            MyMixerHam, prog_init_three_XYMixers = three_XYMixers_Dicke_states(n_q, Ansatz_type['n_destinations']) 
            Ansatz_type['prog_init'] = 'Three Dicke states'
            Ansatz_type['Mixer'] ='Three XY Ring mixers' 
            labels= ['Destination XYmixer','Start0 XYmixer','Start1 XYmixer']
            for n_mixer, mixer in enumerate(MyMixerHam):
                if print_details: 
                    print(labels[n_mixer])
                    print(mixer)
    else:
        for n in range(n_q):       
            MyMixerHam[0] += 1*sX(n) 
        Ansatz_type['Mixer'] ='TRADITIONAL /  Farhi  mixer' 
        if print_details: 
            print(Ansatz_type['Mixer'])
            print(MyMixerHam[0])        
        
    return MyMixerHam#, Ansatz_type

def get_prog_init(n_qubits,use_XY_mixer,n_ambulance,HammingWeightOfConstraint,n_destinations, Ansatz_type=None ):
    """
    Return the Quantum circuit of the initial state (eg Hamamard for Xmixer)
    if use_XY_mixer==1 then the circuit created is a single Dicke_state otherwise it is an Hamamard
    if n_ambulance ==2  and use_XY_mixer==1, then 3  Dicke_states are created otherwise it is       a Hadamard.
    """
    prog_init = Program()

    if not use_XY_mixer:
        #A) 
        print('Hadamard initial state.   ', n_qubits, ' = n_qubits')
        state_init = ''
        for n in range(n_qubits):
            prog_init += H(n)
            state_init += 'H'
        if not Ansatz_type==None:Ansatz_type['prog_init'] = 'Pure start state |' + state_init +'>'
    else:
        if n_ambulance==2:
            #B)
            MyMixerHam, prog_init_three_XYMixers = three_XYMixers_Dicke_states(n_qubits, n_destinations)  
            prog_init = prog_init_three_XYMixers       
            print('The initial state is made of three separate Dicke_state s')
            if not Ansatz_type==None:Ansatz_type['prog_init'] = 'three separate Dicke_state s'
        else:
            #C)
            print('Dicke_state initial state',HammingWeightOfConstraint, ' = HW')
            prog_init = Dicke_state(n_qubits, HammingWeightOfConstraint)
            if not Ansatz_type==None:Ansatz_type['prog_init'] = 'Dicke_state of Hamming weight = ' + str(HammingWeightOfConstraint)
    return prog_init
def get_prog_init(prt_details = 0, **Ansatz_type ):
    """
    Return the Quantum circuit of the initial state (eg Hamamard for Xmixer)
    if use_XY_mixer==1 then the circuit created is a Dicke_state otherwise it is an Hamamard
    if n_ambulance ==2 then and use_XY_mixer==1, then 3  Dicke_states are created otherwise it is       a Hamamard.
    """
    HW= Ansatz_type['HammingWeightOfConstraint']
    n_q = Ansatz_type['n_qubits']
    prog_init = Program()
    #n_qubits,use_XY_mixer,n_ambulance,HammingWeightOfConstraint,n_destinations,
    if not Ansatz_type['use_XYMixer_constraints']:
        #A) 
        if prt_details: print('Hadamard initial state.   ',n_q , ' = n_qubits')
        state_init = ''
        for n in range(n_q):
            prog_init += H(n)
            state_init += 'H'
        Ansatz_type['prog_init'] = 'Pure start state |' + state_init +'>'
    else:
        if Ansatz_type['n_ambulance'] ==2:
            #B)
            MyMixerHam, prog_init_three_XYMixers = three_XYMixers_Dicke_states(n_q, Ansatz_type['n_destinations'])  
            prog_init = prog_init_three_XYMixers       
            if prt_details: print('The initial state is made of three separate Dicke_state s')
            Ansatz_type['prog_init'] = 'three separate Dicke_state s'
        else:
            #C)
            if prt_details: print('Dicke_state initial state',HW, ' = HW')
            prog_init = Dicke_state(Ansatz_type['n_qubits'], HW)
            Ansatz_type['prog_init'] = 'Dicke_state of Hamming weight = ' + str(HW )
    return prog_init

from pyaqc.RCModules.QAOARC import decimal_state_to_binary_reversed, Energy_of_binary_state
from pyaqc.RCModules.aqc_rigetti import min_energy
def get_QAOA_circuit(p_init, prt_details=1, state_feasible=None,**Ansatz_type):
    #n_qubits, n_destinations,n_ambulance,use_XYMixer_constraints,HammingWeightOfConstraint,Adjacency,Adjacency_constraint,Adjacency_feasible,
    """
    Return the quantum circuits that comprise the QAOA circuit and the observables that represent the edges and nodes of the problem.
    """
    # Get the Quantum circuit of the initial state (eg Hamamard for Xmixer)
    prog_init = get_prog_init( **Ansatz_type  )
    # Get Mixer quantum circuit
    MyMixerHam = get_mixer( print_details=prt_details,**Ansatz_type)
    # Get Pauli observables that represent an Ising problem, derived from an Adjacency table using a qubo definition of edge.
    ListPauli_terms,SumPauli_terms  = Adjacency_qubo_to_Regetti( Ansatz_type['Adjacency'], Ansatz_type['n_qubits'])
    #Create the QAOA program from prog_init,MyMixerHam,ListPauli_terms,SumPauli_terms 
    ansatz_prog = ansatz_prog_init(prog_init, ListPauli_terms, Ansatz_type['n_qubits'], Ansatz_type['n_destinations'],   p_init, MyMixerHam) 
    return prog_init,MyMixerHam,ListPauli_terms,SumPauli_terms,ansatz_prog


def get_gnd_state_probs_and_approx_ratio_simple_init(p_init, prt_details=1, state_feasible=None,**Ansatz_type):
    """
    Returns a function for calculating prob_gnd_state, approx_ratio, prob_feasible_state

    param Ansatz_type (type:dict) supplies the details necessary to create the QAOA circuit with get_QAOA_circuit(), 
        eg, Ansatz_type['Adjacency'],Ansatz_type['n_qubits'],Ansatz_type['use_XYMixer_constraints'],Ansatz_type['n_destinations']
    """
    #UseBetaGamma == { 'Betas':3, 'Gammas':3}
    # n_qubits, n_destinations,n_ambulance,use_XYMixer_constraints,HammingWeightOfConstraint,Adjacency,Adjacency_constraint,Adjacency_feasible,
    prog_init,MyMixerHam,ListPauli_terms,SumPauli_terms,ansatz_prog = get_QAOA_circuit(p_init, prt_details=0, state_feasible=state_feasible,**Ansatz_type)
    
    get_gnd_state_probs_and_approx_ratio = get_approx_ratio_init(prog_init,MyMixerHam,ListPauli_terms,SumPauli_terms,ansatz_prog,Ansatz_type['n_qubits'], Ansatz_type['n_destinations'],p_init,Ansatz_type['Adjacency_constraint'],Ansatz_type['Adjacency_feasible'], prt_details=prt_details, state_feasible=state_feasible)
    
    return get_gnd_state_probs_and_approx_ratio
def get_gnd_state_probs_and_approx_ratio_VQE_init(prog_VQE_ansatz, prt_details=1, state_feasible=None,**Ansatz_type):
    prog_init = None
    p_init = None
    MyMixerHam=None
    use_VQE=1
    
    ListPauli_termsMy,SumPauli_termsMy = Adjacency_qubo_to_Regetti( Ansatz_type['Adjacency'], Ansatz_type['n_qubits'])
    return get_approx_ratio_init(\
        prog_init,MyMixerHam,ListPauli_termsMy,SumPauli_termsMy,prog_VQE_ansatz,Ansatz_type['n_qubits'],\
        Ansatz_type['n_destinations'],p_init,Ansatz_type['Adjacency_constraint'],Ansatz_type['Adjacency_feasible'],\
        prt_details=prt_details, state_feasible=state_feasible, show_debug=False, Betas=1,Gammas=1, use_VQE=use_VQE)
    
    
def get_approx_ratio_init(prog_init,MyMixerHam,ListPauli_termsMy,SumPauli_termsMy,ansatz_prog,n_qubits, n_destinations,p_init,Adjacency_constraint,Adjacency_feasible, prt_details=1, state_feasible=None, show_debug=False, Betas=1,Gammas=1, use_VQE=0):
    """
    Returns get_gnd_state_probs_and_approx_ratio() 

    The energy of , Adjacency_feasible (the X-mixer constraint adjacency), in a state_feasible, will be the characteristic energy that all feasible states share. 
    
    The energy of Adjacency_constraint in its ground state, is the energy cost of complying with all constraints in the problem
    
    Initialise the function, 'get_gnd_state_probs_and_approx_ratio(opt_betagamma)' , in one of two ways
        1) Each call will be for the same ising problem. Hence, it will use the same E_min_feas,  E_max_feas, and min_energy_constraint
        2) each call of 'get_gnd_state_probs_and_approx_ratio' can use a different betagamma,  otherwise a p based on opt_betagamma
    """
    
    if isinstance(state_feasible,str):
        state_feasible=state_feasible.replace(' ','')
    "get  Minimum energy of the Adjacency by considering all 2^n_qubits states. Or more quickly using a suppied state_feasible"
    energy_of_feasible_state = min_energy( Adjacency_feasible, n_qubits, state_feasible=state_feasible)      ##use An adjacency that defines feasible space
    #min_energy_constraint_problem = min_energy( Adjacency_constraint, n_qubits, state_feasible=state_feasible)              ##use problem's actual mixer_constraints
    
    ListPauli_terms_constr,SumPauli_termsMy_constr  = Adjacency_qubo_to_Regetti( Adjacency_constraint, n_qubits)
    ### Calculate E_min_feas and E_max_fease  by a 2^n search ##############################################
    ListPauli_terms_feas,SumPauli_terms_feas  = Adjacency_qubo_to_Regetti( Adjacency_feasible, n_qubits)
    def energy_min_max_feasible():
        seen_first_compliant_state = 0
        for n in range(2**n_qubits):
            E_feas_of_state = Energy_of_binary_state( n, SumPauli_termsMy,n_qubits)#-min_energy_constraint_problem
            
            if Energy_of_binary_state( n, SumPauli_terms_feas, n_qubits) == energy_of_feasible_state:
                if seen_first_compliant_state==0:
                    Energy_frm_constraint = Energy_of_binary_state( n, SumPauli_termsMy_constr,n_qubits) # This is the same value every feasible state, as it is part of feas definition
                    E_max_feas = E_feas_of_state - Energy_frm_constraint
                    E_min_feas = E_feas_of_state -Energy_frm_constraint
                    seen_first_compliant_state = 1
                    
                E_max_feas = max(E_max_feas,E_feas_of_state - Energy_frm_constraint)
                E_min_feas = min(E_min_feas,E_feas_of_state - Energy_frm_constraint)
                
        return E_min_feas, E_max_feas,Energy_frm_constraint 

    E_min_feas, E_max_feas,Energy_frm_constraint = energy_min_max_feasible()
    ##############################################
    if prt_details==1:
        print('E_min_feas = ', E_min_feas,'\nE_max_feas = ', E_max_feas,'\n energy_of_feasible_state = ', energy_of_feasible_state,'\n Energy_frm_constraint in feas state = '
        ,Energy_frm_constraint  )
    
    def get_gnd_state_probs_and_approx_ratio(opt_betagamma,  prt_details=1):  
        """prob_gnd_state, approx_ratio, prob_feasible_state
         
         return A)'prob_gnd_state' the cumulative probability of the lowest energy state(s) sampled, 
                
                B) 'prob_feasible_state' The probability that one of the suggested states is a feasible state (defined as states that minimize 'Adjacency_feasible') 
         
                C) approx_ratio = (EV_feasible/prob_feasible_state-E_max_feas) /(E_min_feas-E_max_feas)
                    Where 
                    E_min_feas = the problem's minimum energy in the feasible space with all constraint energies set to 0.
                    E_max_feas = the problem's maximum energy in the feasible space with all constraint energies set to 0.
                    EV_feasible = the expected energy of suggested states that are feasible
                        As defined by Wang in arXiv:1904.09314v2 [quant-ph] 21 May 2020, page 3 col 2
        type: float,float,float

        param
            ansatz_prog a quantum circuit eg for p=1, e^{beta * H_mixer} e^{gamma * H_cost} prog_init
                type:function that returns a Rigetti.Program()
            'Adjacency_constraint'  The adjacency table which contains only the problem's constraint energies , and no Soft energies (eg distance)
            This allows the constraint energies to be calculated separately from the non-constraint energies.
                type: dict example {(0, 0): 5,  (1, 1): -1,  (2, 2): -1}
            
        kwargs
            
            'state_feasible' a string of a state that complies with all the constraints in Adjacency_constraint. For example;  '11011'. When supplied the minimum energy of                 all the constraints is calculated in one calculation, rather than 2^n_qubits calculations. This speeds up the approximation ratio calculation
                type: string of a binary state
            
            
        if p == len(opt_betagamma)//2=False, the ansatz_prog is re initiated each call of get_gnd_state_probs_and_approx_ratio, to reflect the p implied by opt_betagamma (eg when creating the approx ratio for many different p s of the same problem)
        if p == len(opt_betagamma)//2=True or use_VQE==1 , execution time is saved by using the ansatz_prog that was supplied by get_approx_ratio_init (eg when looking for betagamma minimum)
        In the feasible space, what is the average energy found after substituting the actual constraint energy with 0? 
        """
        prob_gnd_state = 0
        prob_feasible_state = 0
        EV_feasible = 0 
        EV_not_feasible = 0 
        n_feasible = 0
        approx_ratio = -1
        EV_all_suggestions = 0
        Expectation_of_state = 0
        # UseBetaGamma == { 'Betas':3, 'Gammas':3}
        #print(Betas, Gammas, 'from init')
        if p_init == len(opt_betagamma)//(Betas+ Gammas) or use_VQE:
            # calculate solution_ansatz with latest opt_betagamma, using the same p as when get_approx_ratio_init first called
            # or if using a VQE ansatz keep the ansatz program first suggested
            solution_ansatz = ansatz_prog(opt_betagamma)        #solution_ansatz is superposition of many states.
        else:
            # calculate solution_ansatz with latest opt_betagamma and its implied p
            p_new = len(opt_betagamma)//2
            ansatz_prog_flexible_p = ansatz_prog_init(prog_init, ListPauli_termsMy, n_qubits, n_destinations,   p_new, MyMixerHam)
            solution_ansatz = ansatz_prog_flexible_p(opt_betagamma)
            
        
        probsAbsolute  = WavefunctionSimulator().wavefunction(solution_ansatz).get_outcome_probs()
        opt_probs_absolute = [x for x in probsAbsolute.values()]
        
        for n in range(2**n_qubits):
        
            if opt_probs_absolute[n] > 0.00001 :    #|qn...q0>
                # calculation of energy, using n as the state, NOT recalculating the state by building its quantum circuit eg Program(X(0),I(1))
                Expectation_of_state = Energy_of_binary_state( n, SumPauli_termsMy,n_qubits)    

            if Expectation_of_state <= E_min_feas +  Energy_frm_constraint:   #test if a ground state
                
                #psi_opt, state_string = decimal_state_to_binary_reversed(n,n_qubits)
                prob_gnd_state += opt_probs_absolute[n]
            
            # If the Adjacency_constraint has an energy that is the minimium possible for a feasible state, the x mixer constraints have been met, ie space is 'feasible'..
            #...then the energy, without any problem_constraint energy, is averaged.
            # Energy_frm_constraint
            if Energy_of_binary_state( n, SumPauli_terms_feas, n_qubits) == energy_of_feasible_state: #use X mixer_c
                EV_feasible += opt_probs_absolute[n] * (Expectation_of_state - Energy_frm_constraint )  #use problem's mixer_c
                prob_feasible_state +=opt_probs_absolute[n]
                n_feasible +=1
            elif show_debug:
                EV_not_feasible += opt_probs_absolute[n] * Expectation_of_state 
            if show_debug:EV_all_suggestions += opt_probs_absolute[n] * Expectation_of_state 
        if show_debug:
            approx_ratio_old = EV_feasible/(prob_feasible_state*E_min_feas)
            print('approx_ratio_old = ',approx_ratio_old,  '\nEV_feasible/prob_feasible_state = ', EV_feasible/prob_feasible_state, '\nEV_not_feasible =',EV_not_feasible,'\nEV_all_suggestions( P>0.00001) = ', EV_all_suggestions,'\nn_feasible =',n_feasible)

        approx_ratio = (EV_feasible/prob_feasible_state-E_max_feas) /(E_min_feas-E_max_feas)
        return prob_gnd_state, approx_ratio, prob_feasible_state
    return get_gnd_state_probs_and_approx_ratio
