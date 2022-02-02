#testDicke_state_local.py

from states import  Dicke_state,Dicke_state_local
from pyquil import *
from pyquil.gates import *
n_qubits =12
prog = Program()
for num in range(n_qubits):
    prog += I(num)

if 1:
    # Make qubits 0 to 2 in a HW=1 state
    qudit_start = 0
    qudit_end=2
    HW=1
    init_prog = Dicke_state_local(prog,qudit_start, qudit_end, HW)

    # Make qubits 3 to 5 in a HW=1 state
    
    qudit_start = 3
    qudit_end=5
    HW=1
    init_prog = Dicke_state_local(init_prog, qudit_start, qudit_end, HW)

    # Make qubits 6 to 11 in a HW=3 state
    qudit_start = 6
    qudit_end=11
    HW = 3
    init_prog = Dicke_state_local(init_prog, qudit_start, qudit_end, HW)



wfn = WavefunctionSimulator().wavefunction(init_prog)
# print states with non-zero probability
for state, prob in wfn.get_outcome_probs().items():
    if(prob>0):
        print(state,  '%3.2f'%prob)
print(wfn.get_outcome_probs() )

print(init_prog)
from pyquil.api import WavefunctionSimulator

#Test senarios
# test Dicke_state_local(prog, 0, 2, HW) output is the same as Dicke_state(3, HW)
if 0:
    qudit_start = 0
    qudit_end=2
    HW=1
    init_prog = Dicke_state(3, HW)
    wfn = WavefunctionSimulator().wavefunction(init_prog)
    print(wfn.get_outcome_probs() )
    prog = Dicke_state_local(prog,qudit_start, qudit_end, HW)
    wfn = WavefunctionSimulator().wavefunction(prog)
    print(wfn.get_outcome_probs() )

if 0:   #test if  qudit_end is  not a qubit number in 'prog'
    prog = Program(I(0),I(3))
    qudit_start = 0
    qudit_end=4
    HW=3
    init_prog = Dicke_state_local(prog,qudit_start, qudit_end, HW)
prog = Program(I(1),I(4))
qudit_start = 0
qudit_end=4
HW=3
if 0:      #test if  qudit_start is  not a qubit number in 'prog'
    init_prog = Dicke_state_local(prog,qudit_start, qudit_end, HW)