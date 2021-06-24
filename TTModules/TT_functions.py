# This is a file of functions created by TT and saved so that they can be imported as a module #
# The function list in this file is
# CreateTwoAmbulanceAdjacency_TT
# distance_TT -     produces the square of the distance calculated by the RC function distance
# co_ord_TT
# SaveListTxt_TT
# OpenListTxt_TT    Assumes a Dict as the first line
# OpenListTxt       original RC open file function - does not allow for a Dict as the first line
# AppendListTxt_TT
# AppendListTxt     Original RC open file Append function for saving a new result to an existing file
# print_details     Original RC print details function
# Ham_wt_of_string(string)  RC function to calculate the Hamming weight of a string
# Expect_of_state(state:int, MyCostHam,n_qubits)  RC function which calculates the Expectation value for a given state for a cost function: MyCostHam

import numpy as np
from pyquil.api import WavefunctionSimulator
import ast
from pyquil import Program
from pyquil import *
from pyquil.gates import *
import numpy as np
from pyquil.paulis import ID, sZ, sX, sY , exponential_map, exponentiate_commuting_pauli_sum

def CreateTwoAmbulanceAdjacency_TT(gridWidth,n_destinations: 'Number of locations servied by ambulance fleet'= 5, Adddistance = 1
,ConstraintMultiplier:'If ConstraintMultiplier= 0, it is calculated automatically'=0,  qubo_model:'qubo_model definition used in Adjacency table'=True):
    """
    This is a variation of the original Adjacency function but uses the square of the previous distance function instead, resulting in bigger differences between the location positions and may result in better optimiser performance in finding low energy states
   
    Returns: type dict result={'qubo':qubo, 'quboHybrid':quboHybrid, 'n_qubits':n_qubits, 'ConstraintMultiplier':ConstraintMultiplier,'max_distance': max_distance, 'sum_distance':sum_distance}
    Returns: type dict eg  result['qubo']={(0,3):4}, is an edge of weight 4 between qubits 0 and 3.  QUBO to represent the ising ambulance problem, it involves two constraints and minimising distance driven.
################# DESCRIPTION OF THE ISING PROBLEM #################
BINARY/qubo MODEL  
This splits the geography into Width*Height (W*H) qubits, and mininises the distance from...
        A0 starting location (described by a W*H qubits labelled 'A0_Start')...
       to destination location (W*H qubits labelled collectively A0_Des )
With the contraints that;
       A0 has only one starting location (ie there is a single one in the W*H qubits collectively labelled 'A0_Start')
       A0 only services the destinations in A0_Des that have a 1
       A0 does not share any destinations with A1
Energy is at a min when  the distance is minimised and the constraints are met.
For example; location grids of 2*2, or 3*2, or 3*3 . These use; 2*2*4 = 16 , or 3*2*4 = 24, or 3*3*4 = 36   qubits respectively.  
"""

    ################# The remainder of this cell creates the adjaceny table ('qubo') that represents the Ising problem 
            


    #Specify the dimensions of the location grid - use 2*2 grid to understand method of creation

    Width  = gridWidth
    Height = n_destinations // gridWidth
    # Each ambulance needs a feature to describe it;
        # a location is to be visited by that Amb ('A0_Des' or A1_Des) and 
        # where that Amb is kept awaiting call out ( 'A0_Start').  
    Feature_options = {0:'A0_Des', 1:'A1_Des',2: 'A0_Start',3:'A1_Start'}
    #Feature_options = {0:'A0_Des', 1:'A1_Des'}      # to test destination constraint
    #Feature_options = {0: 'A0_Start'}#,1:'A1_Start'}  # to test START constraint
    n_grids = len(Feature_options)
    n_qubits =  Width * Height * n_grids
    n_qubitsHybrid =  Width * Height * 2 + 2
    num_of_ones_constraint = 1                 # Defined as number of ambulances that service each destination (normally 1)
    #Weighting that control labelling which ambulance will use which DESTINATIONS in the,  Width * Height, grid 
    # These constraints are based on a qubo model when a qubit = 0 or 1 not Ising def of -1 and 1  
    if 1:
        Des_quadWt = 2                  # adds energy if A1_des is same a A0_des
        Des_linearWt = -1               # reduces energy if a destination is removed
    else:
        Des_quadWt = 0                  # adds energy if A1_des is same a A0_des
        Des_linearWt = 0               # reduces energy if a destination is removed
    if 1:
        Start_quadWt = 2                  # adds energy if A1_Start is same a A0_Start
        Start_linearWt = 1-2*num_of_ones_constraint               # this is the constraint condition for a binary/qubo model (x= 0 or 1)
    else:
        Start_quadWt = 0                  # adds energy if A1_Start is same a A0_Start
        Start_linearWt = 0               # reduces energy if a Start is removed


    ################# Create parameters for two different binary/qubo model. a) Complete Ising b) Ising where the starting position of each ambulance is fixed and known; 'quboHybrid'. 
    qubo   = {(0,0):0}
    quboHybrid   = {(0,0):0}
    
    # Create adjacency table of linear weights on the diagonal and edge weight in the top triangle above the diagonal
    ################# FIRST calculate the sum of all the distances between destinations to scale up the energy of constraint energies
    sum_distance = 0
    max_distance = 0
    
    ### Note that functions co_ord_TT and distance_TT are used
    for q1 in range(0,n_qubits):
        feature1,r1,c1 = co_ord_TT(q1,Height, Width,Feature_options)
        for q2 in range(q1,n_qubits):
            feature2,r2,c2 = co_ord_TT(q2,Height, Width,Feature_options)
            if feature1 == 'A0_Des' and feature2 == 'A0_Start':
                if  (r2 * Width )+ c2 > (r1 * Width) + c1 or 1:
                    sum_distance += distance_TT(r1,c1,r2,c2)
                    max_distance = max(max_distance, distance_TT(r1,c1,r2,c2))
    
    #Use sum of distances to ensure any constraint energy is larger than any distance energy
    #sum_distance = 0
    if ConstraintMultiplier == 0:
        if 0 :
            ConstraintMultiplier = sum_distance 
        else:
            ConstraintMultiplier = max_distance * 10#0      #Lagrange parameter
    ################# SECOND create adjaceny maxtrix
    for q1 in range(0,n_qubits):
        feature1,r1,c1 = co_ord_TT(q1,Height, Width,Feature_options)
        for q2 in range(q1,n_qubits):
            feature2,r2,c2 = co_ord_TT(q2,Height, Width,Feature_options)
            ##################  Constraint that every DESTINATION is visited by exactly one ambulance  #########################  
            if c1 == c2 and r1==r2 and feature1 == 'A0_Des' and feature2 == 'A1_Des':
                qubo[(q1,q2)] = ConstraintMultiplier * (Des_quadWt )
                if q1 < n_qubitsHybrid and q2 < n_qubitsHybrid:
                    quboHybrid[(q1,q2)] = qubo [(q1,q2)]
            if q1 == q2 and (feature1 == 'A0_Des' or feature1 == 'A1_Des'):
                qubo[(q1,q2)] = ConstraintMultiplier * Des_linearWt
                if q1 < n_qubitsHybrid and q2 < n_qubitsHybrid:
                    quboHybrid[(q1,q2)] = qubo [(q1,q2)]
            ##################  Constraint that each ambulance has just one START location #########################  
            if  (feature1 == 'A0_Start' or feature1 == 'A1_Start') and feature1==feature2:            
                if(c1 == c2 and r1 == r2):
                    qubo[(q1,q2)] = ConstraintMultiplier * Start_linearWt
                    if q1 < n_qubitsHybrid and q2 < n_qubitsHybrid:
                        quboHybrid[(q1,q2)] = qubo [(q1,q2)]
                else:
                    qubo[(q1,q2)] = ConstraintMultiplier * Start_quadWt
            ##################  Constraint that ambulance 0 minimises START to destination distance ######################### 
            
            ### Note that function distance_TT is used ###
            if feature1 == 'A0_Des' and feature2 == 'A0_Start' or feature1 == 'A1_Des' and feature2 == 'A1_Start':
                if  (r2 * Width )+ c2 > (r1 * Width) + c1 or 1:
                    d = 0 
                    if Adddistance:
                        d = distance_TT(r1,c1,r2,c2)
                    if d!=0:
                        qubo[(q1,q2)] = d
                #else:
                 #   qubo[(q1,q2)] = 0 
    result={'qubo':qubo, 'quboHybrid':quboHybrid, 'n_qubits':n_qubits, 'ConstraintMultiplier':ConstraintMultiplier,'max_distance': max_distance, 'sum_distance':sum_distance}
    return result

def distance_TT(r1,c1,r2,c2):
    #Calc the distance between node i and node j. ASSUME, square gird with n nodes, and distance of 1 between neighbouring nodes
    """
    THIS IS A VARIATION ON THE PREVIOUS DISTANCE FUNCTION AND GENERATES (distance)**2 which may help with optimiser performance for this 2 ambulance problem
    """
    y = (r1-r2)**2
    x = (c1-c2)**2
    return (x+y)**2             #distance squared

def co_ord_TT(qubit,Height, Width,Feature_options):
    """
    THIS IS A COPY OF THE co_ord function CREATED BY RC
    return type:int,int  row, column in the position grid of dimension Width * Height from the qubit number , and which of 4 features is stored there
    r,c describes what location the feature_num is attached to.
    """
    feature_num = 0
    feature_num = qubit//( Width * Height)
    r = (qubit//Width) %  Height      # // means division followed by integer truncation
    c= qubit % Width

    return Feature_options[feature_num],r,c


def SaveListTxt_TT(filename: '.txt file', file_details: 'Dict specifying details being saved', ListToSave: 'A list of lists of results for each starting betagamma'):
    """
    Create a line string populated by  1) a Dict which specifies key variables which have been used to generfate the results and 2) the sequential elements of ListToSave,
    separated by ' ' . Store these line strings in a .txt file
    
    For example
    SaveListTxt('temp.txt', Ansatz_type, [
                            [-1050.451, 0.195, 4.0, [0.868, 0.876], [0.131, 0.932]  ],      #FIRST line in the .txt
                            [-1050.451, 0.195, 4.0, [0.868, 1,0.876], [0.131,2, 0.932]  ]   #SECOND line in the .txt
                                                                                    ])
    creates a file called 'temp.txt' whose first line is a Dict 
    Ansatz_type, e.g. {'n_qubits': 16, 'p': 2, 'Mixer': 0, 'prog_init': 'Pure start state |HHHHHHHHHHHHHHHH>', 'ConstraintMultiplier': 90, 'remove_constraint': 0, 'Optimizer_maxiter': 10}
    and second line is
    -1050.451, 0.195, 4.0, [0.868, 0.876], [0.131, 0.932]
    and third line is
    -1050.451, 0.195, 4.0, [0.868, 1,0.876], [0.131,2, 0.932]
    """
    f = open(filename, "w+")    #erases files if it already exists
    f.writelines(str(file_details)+' ')
    f.writelines('\n') 
    for d in range(len(ListToSave)):
        for n, prob in enumerate(ListToSave[d] ):
            f.writelines(str(prob)+' ')
        f.writelines('\n')              #marks the end of a line
    f.close()


def SaveListTxt(filename: '.txt file', ListToSave: 'A list of lists of results for each starting betagamma'):
    """
    Create a line string populated by the sequential elements of ListToSave,  separated by ' ' . Store these line strings in a .txt file
    
    For example
    SaveListTxt('temp.txt',[
                            [-1050.451, 0.195, 4.0, [0.868, 0.876], [0.131, 0.932]  ],      #FIRST line in the .txt
                            [-1050.451, 0.195, 4.0, [0.868, 1,0.876], [0.131,2, 0.932]  ]   #SECOND line in the .txt
                                                                                    ])
    creates a file called 'temp.txt' whose first line is
    -1050.451, 0.195, 4.0, [0.868, 0.876], [0.131, 0.932]
    and second line is
    -1050.451, 0.195, 4.0, [0.868, 1,0.876], [0.131,2, 0.932]
    """
    f = open(filename, "w+")    #erases files if it already exists
    for d in range(len(ListToSave)):
        for n, prob in enumerate(ListToSave[d] ):
            f.writelines(str(prob)+' ')
        f.writelines('\n')              #marks the end of a line
    f.close()

def is_number(s):
    """" Return True if s is a number not a string
        type:bool"""
    try:
        float(s)
        return True
    except ValueError:
        return False

def OpenListTxt_TT(filename, decimal_places=3):
    """
    Filename is presumed to be headed by a Dict followed by a list containing numbers and lists
    """
    import ast
    f = open(filename, "r+")        #does not erase file just reads it.
    DataOut = []
    sub_list = []
    liststr = ''
    string = f.readlines()              #read all the lines in the file
    #print(string[0])                               # line by line  and create a list using the space delimiter
    dict_str = string[0]
    convertedDict = ast.literal_eval(dict_str)
    string1 = string[1:]
    print(len(string1), 'string1')
    for n, line in enumerate(string1):
        #print('n', n, line, 'line')#, end='')
        line.rstrip()
        list_of_new_line = []
        listfnd = 0
        for x in line.split(' '):  #A)strip end white space B) split the string by the delimiter C) convert to float list:
            if x != '\n':
                if  x.rfind('[')>=0 or listfnd  > 0 :   #start of a sub_list
                    
                    listfnd += 1            # ADD number of [
                    y = x.strip('[')
                    y = y.strip(']')
                    y = y.strip(',')
                    if is_number(y):sub_list.append( float(('%3.'+str(decimal_places)+'f')% float(y)) )
                    if x.rfind(']') >=0:                #end of a sub_list
                        listfnd -= 1        # SUB number of ]
                        list_of_new_line.append(sub_list)
                        sub_list = []
                else:
                    if is_number(x):list_of_new_line.append( float(('%3.'+str(decimal_places)+'f')% float(x)) )
            else:
                DataOut.append( list_of_new_line)
    f.close()
    return convertedDict, DataOut


def OpenListTxt(filename, decimal_places=3):
    """
    filename is presumed to be made up of a list containing numbers and lists 
    This is RC's original OpenListTxt function and has not been amended
    """
    f = open(filename, "r+")        #does not erase file just reads it.
    DataOut = []
    sub_list = []
    liststr = ''
    string = f.readlines()              #read all the lines in the file
                                     # line by line  and create a list using the space delimiter
    for n, line in enumerate(string):
        #print(line, 'line')#, end='')
        line.rstrip()
        list_of_new_line = []
        listfnd = 0
        for x in line.split(' '):  #A)strip end white space B) split the string by the delimiter C) convert to float list:
            if x != '\n':
                if  x.rfind('[')>=0 or listfnd  > 0 :   #start of a sub_list
                    
                    listfnd += 1            # ADD number of [
                    y = x.strip('[')
                    y = y.strip(']')
                    y = y.strip(',')
                    sub_list.append( float(('%3.'+str(decimal_places)+'f')% float(y)) )
                    if x.rfind(']') >=0:                #end of a sub_list
                        listfnd -= 1        # SUB number of ]
                        list_of_new_line.append(sub_list)
                        sub_list = []
                else:
                    list_of_new_line.append( float(('%3.'+str(decimal_places)+'f')% float(x)) )
            else:
                DataOut.append( list_of_new_line)
    f.close()
    return DataOut



def AppendListTxt_TT(filename:'file to be appended to', ListToSave: 'A list of lists of results for each starting betagamma'):
    """
    Add on to the end of an existing file
    Create a line string, terminated with \n, populated with the elements of ListToSave,  separated by ' ' . Store these line strings in a .txt file
    Same as AppendListTxt function in terms of code
    """
    f = open(filename, "a+")    #"w+" erases files if it already exists
    for d in range(len(ListToSave)):
        for data in ListToSave[d] :
            f.writelines(str(data)+' ')
        f.writelines('\n')              #marks the end of a line
    f.close()


def AppendListTxt(filename:'file to be appended to', ListToSave: 'A list of lists of results for each starting betagamma'):
    """
    Add on to the end of an existing file
    Create a line string, terminated with \n, populated with the elements of ListToSave,  separated by ' ' . Store these line strings in a .txt file
    """
    f = open(filename, "a+")    #"w+" erases files if it already exists
    for d in range(len(ListToSave)):
        for data in ListToSave[d] :
            f.writelines(str(data)+' ')
        f.writelines('\n')              #marks the end of a line
    f.close()





def print_details(quboUnordered,n_qubits,name):
    #Print adjacency graph of weights Wij between one qubit index (Qi)  and another (Qj)
    i = 0                                   # index of row of adjacency matrix
    j = 0                                   # index of column
    print( 'ADJACENCY MATRIX',name,' \n\n\t', end='')
    qubo = {}

    #order the dicts by the first element of each tuple, this helps the printing order
    for key in sorted( quboUnordered.keys() ):
        qubo[key] = quboUnordered[key]
    for n in range(0,n_qubits):
        print('q{}'.format(n),'\t', end='') #'Qubit index'
    print('\n\n', end='')
    j=-1
    for pos, weight in qubo.items():        # 'weight' is the weight between qubit i and qubit j pos(i,j) 
        if pos[0] > i :                     #row
            j=-1        
        for m in range(i,pos[0]):
            print('\n', end='')
        i = pos[0]
        if j==-1:
            print('q{}'.format(i),'\t', end='')
            j=0
        for n in range(j,pos[1]):
            print('\t', end='')
        j = pos[1]
        print(weight, end='')   

def Ham_wt_of_string(string):
    h=0
    for m in range(len(string)):
        h +=int(string[m])
    return h

def Expect_of_state(state:int, MyCostHam, n_qubits):
    """
    #return (type: float) expectation of a single binary quantum state applied to MyCostHam (type:list of pauli_term or pauli_sum)

    # state int      such that bin(state)=          |qn...q0>  expressed as a string eg '0b111'
    #Method uses the analytic WavefunctionSimulator().expectation() not a sampling approach
    """
    prog_init = Program()
    # From decimal generate sequential binarys to describe a state 001,010,011..etc
    
    state_string = ''
    
    if isinstance(state, int):
        # convert Decimal to binary representation in a list of integers of length n_qubits
        zeros =  n_qubits - len (bin(state)[2:])        # missing 0s to ensure the state description is n_qubits long
        psi_opt = []
        for p in range(0,zeros):
            psi_opt.append(0)
        psi_opt.extend([int(n) for n in bin(state)[2:]])
        for elem in reversed(psi_opt):
            state_string += str(elem)
        state = psi_opt
    #print(psi_opt)
    for n, q in enumerate(reversed(state)) :
        if q==1:
            prog_init += X(n)
        
    expectation_weighted_sum = 0
    #print(prog_init)
    for nth_obs in range( len(MyCostHam)):    
            expectation = WavefunctionSimulator().expectation(prep_prog=prog_init, pauli_terms=MyCostHam[nth_obs])
            expectation_weighted_sum += expectation 
            if state==0 and 0:
                print('{0:18} {1:}'.format(expectation,expectation_weighted_sum), end=' ')
    return np.real(expectation_weighted_sum)