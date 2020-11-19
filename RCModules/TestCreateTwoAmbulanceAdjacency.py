# test PlatformConversion.py functions.
from PlatformConversion import *# print_QUBOdetails, CreateTwoAmbulanceAdjacency

n_destinations = 5  
gridWidth = 5

result = CreateTwoAmbulanceAdjacency(gridWidth, n_destinations, Adddistance = 1, remove_constraint=False)

filename = 'Twoambulances.txt'

# Output


qubo = result['qubo']
quboHybrid = result['quboHybrid']
n_qubits = result['n_qubits']
ConstraintMultiplier = result['ConstraintMultiplier']
max_distance = result['max_distance']

print_QUBOdetails(result['qubo'],n_qubits,filename)