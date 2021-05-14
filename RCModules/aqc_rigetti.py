# aqc_rigetti.py
# aqc functions that rely on the rigetti pyquil funtions
from pyquil import *
from pyquil.gates import *
import numpy as np
from pyquil.api import WavefunctionSimulator
from pyquil.paulis import ID, sZ, sX, sY , exponential_map, exponentiate_commuting_pauli_sum
"""
def Dicke_state_local( 
p: 'Type:Program()',
qudit_start:'FIRST qubit in ring mixer(Type: int)'=0, 
qudit_end: 'LAST qubit in ring mixer (Type: int)'=2,  HWeight: 'int'=1):
def Dicke_state(n_qubits: int , HWeight: int):

def create_local_ring_mixer( qudit_start:'FIRST qubit in ring mixer'=0, qudit_end:'LAST qubit in ring mixer'=2):
def create_ring_mixer(n_qubits):

def Adjacency_Ising_to_Regetti( Adjacency: 'Adjacency table in form: qubo[(r,c)]', n_qubits):
def Adjacency_qubo_to_Regetti( Adjacency: 'Adjacency  table in form: qubo[(r,c)]', n_qubits):
 
"""
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
    #XYmixer = 0.5*sX(0)*sX(n_qubits) + 0.5*sY(0)*sY(n_qubits)           #old seemed to work but adds XY0,n_qubits
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
