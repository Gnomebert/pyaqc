# QAOARC.py
# For a measured waveform calculate its energy from the PauliSum of observables
from pyquil.paulis import ID,sZ,sX,sY, PauliSum, PauliTerm
from pyquil.api import * #
from pyquil.api import local_forest_runtime
from pyquil.gates import CNOT, I,Z,X,Y,RZ,RY,RX,H, SWAP, CSWAP, PHASE
from pyquil import Program
import numpy as np
__all__=['Energy_of_binary_state', 'parity_even_p', 'decimal_state_to_binary_reversed']
def Energy_of_binary_state(binary_state:'defined as list[qn...q0] or an interger', pauli_sum:' Pauli_sum eg 0.5-sZ(10',n_qubits):
    """
    Returns (type:float) the energy of the binary_state operating on the sum of Pauli terms,
            (type:string), state_string, a string of the binary state defined as q0...qn
    params
        'binary_state' either 1) a list [qn...q0] of the state whose energy is to be calculated or 2) an integer whose binary value represents that state
    """
    state_string = ''
    #Convert state, (type:list), to an integer if the argument is not already an integer
    if isinstance(binary_state, list):
        b_state = 0
        for n, elem in enumerate( reversed(binary_state) ):
            if elem:
                b_state |= 1 << n       #bitwise operator qn...q0 
    else:
        b_state = binary_state      # qn...q0 

    samples = 1 
    expectation = 0
    for j, term in enumerate(pauli_sum.terms):
            #meas_basis_change = pq.Program()    # start building a circuit to rotate if measuring a X or Y observable
            qubits_to_measure = []          #RC list of qubit that a Pauli_term/observable operates upon
            meas_outcome = 0
            # RC when there is edges between the term's qubits and another terms qubits these are ????
            #if term.id() == "":     #RC eg within sZ(0)+ 0.5 then the 0.5 would be expectation= 1 * coef, where coef=0.5
            if len( term.operations_as_set())==0:
                meas_outcome = 1.0
                #print( 'id term * np.real(term.coefficient) =',np.real(term.coefficient)*meas_outcome)
            else:                           # eg 3*sZ(0) * sZ(2) , qubits_to_measure[0,2], term.coefficient =3.0, no RY or RX
                for index, gate in term:    # from the term the type of gate and the qubit number it operates on is extracted
                    qubits_to_measure.append(index)     #record all qubits operated upon by this Pauli_term
                    
                    #I think this should NOT be indented so the expectation is only evaluated when all the qubits in the observable have been identified
                expectation_of_single_term = 0
                if 1:       #for bitstring, count in freq.items():
                    # following the structure of pyvqe\vqe.py\expectation, this 'if' statement could be a repeated measurement loop, the count would
                    # be the frequency of the state after samples measurements
                    count=1
                    # logical AND it with the solution, 'binary_state', defined as |qn...q0>
                    even = parity_even_p(b_state, qubits_to_measure)
                    if even:
                        expectation_of_single_term += float(count)/samples             #odd number of z gates in the observable have excited qubits in the ansatz
                    else:
                        expectation_of_single_term -= float(count)/samples
                meas_outcome = expectation_of_single_term       # sampling before applying the weighting/coefficient                 
            expectation += term.coefficient * meas_outcome

    return expectation.real
def parity_even_p(state, marked_qubits):
    """
    Calculates the parity of elements at indexes in marked_qubits

    Parity is relative to the binary representation of the integer state.

    :param state: The wavefunction index that corresponds to this state.
    :param marked_qubits: The indexes to be considered in the parity sum.
    
    :returns: A boolean corresponding to the parity.
    """
    assert isinstance(state, int), "{} is not an integer. Must call " \
                                   "parity_even_p with an integer " \
                                   "state.".format(state)
    mask = 0
    for q in marked_qubits:     # example observable  marked_qubits= [0,2] so mask = 0b101, the binary int that represent a sZ(2) and sZ(0)
        mask |= 1 << q          # mask become a binary number |100> would mean q2 is part of the observable
    # RC convert to a string to count the number of excited qubits common to the ansatz (state) and the observable (marked_aubits)
    #Ansatz cases: a) 20% of time|101> or b) 80% of time|011>
    #RC case a) mask & state =  0b101 ie even so returns 0 case b) mask & state = 0b001 ie odd so returns 1
    # this is the correct outcome from an observable that is an edge
    return bin(mask & state).count("1") % 2 == 0

def decimal_state_to_binary_reversed(state:'(type:int) decimal ref to a state qn...q0', n_qubits):
    """

    Returns type: List[int] 'state' represented as the digits of its binary representation [qn,...,q0] then reverse the order to [q0...qn] ,
            type:string  state defined as 'q0...qn'
    
    Arguments
        state, type:int a decimal integer that when expressed as n_qubit binary represents the state of each qubit qn...q0
        n_qubits (type:int), the number of qubits in the system being represented
    """        
    # convert Decimal to binary representation in a list of length n_qubits
    psi_opt = []
    state_string = ''
    zeros =  n_qubits - len (bin(state)[2:])        # missing 0s to ensure the state description is n_qubits long
    for p in range(0,zeros):
        psi_opt.append(0)
    psi_opt.extend([int(n) for n in bin(state)[2:]]) 
    psi_opt.reverse()  # defined as [q0...qn] this is the list to be returned
    # and into a string defined as 'q0...qn'
    for elem in (psi_opt):
        state_string += str(elem)
    state = psi_opt 
    return psi_opt, state_string
