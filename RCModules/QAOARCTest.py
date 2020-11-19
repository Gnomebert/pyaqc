#QAOARCTest.py
from QAOARC import  Energy_of_binary_state
from pyquil.paulis import ID,sZ,sX,sY, PauliSum, PauliTerm
#state_to_measure = #'defined as list[qn...q0] or an interger'
EV = Energy_of_binary_state( [1,1,0],sZ(1)*sZ(2)+0,3) #'defined as list[qn...q0] or an interger'
assert(EV==1)
print(EV)
EV = Energy_of_binary_state( [0,1,0], 2*sZ(1)*sZ(2)+0,3) #'defined as list[qn...q0] or an interger'
assert(EV==-2)

print(EV)
