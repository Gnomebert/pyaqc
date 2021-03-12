# 

# SAME adjacency array as Dwave equivalent Ising model:  IsingDwaveAmbulance.ipynb
# CREATE n_qubits * n_qubits adjacency table, based on a constraint of the number of 1s in the solution
import numpy as np

#__all__ = ['testfunc','distance','print_details', 'CreateAmbulanceAdjacency']
#My Methods used:
"""
def distance(i:'first qubit',j:'second qubit',n_qubits,gridWidth:'use a grid width to rank the qubits'=0):
def print_details(quboUnordered,n_qubits,name)
def CreateAmbulanceAdjacency(gridWidth,n_qubits = 5, ConstraintMultiplier = 1, Adddistance = 1,remove_constraint: bool =True,
HammingWeightOfConstraint=1, qubo_model:'qubo_model definition used in Adjacency table'=True ):

        pyquil specific: should be aqc_rigetti.py TRANSFERRED to aqc_rigetti.py
def Dicke_state_local( 
p: 'Type:Program()',
qudit_start:'FIRST qubit in ring mixer(Type: int)'=0, 
qudit_end: 'LAST qubit in ring mixer (Type: int)'=2,  HWeight: 'int'=1):
def Dicke_state(n_qubits: int , HWeight: int):

    Pure python files - should be platform conversion
def Adjacency_to_Regetti_qubo( Adjacency: 'Adjacency  table in form: qubo[(r,c)]', n_qubits):
def Adjacency_to_Regetti_Ising( Adjacency: 'Adjacency table in form: qubo[(r,c)]', n_qubits):
    pyquil specific:    should be aqc_rigetti.py TRANSFERRED to aqc_rigetti.py
def create_local_ring_mixer( qudit_start:'FIRST qubit in ring mixer'=0, qudit_end:'LAST qubit in ring mixer'=2):
def create_ring_mixer(n_qubits):
    Pure python files
def OpenListTxt(filename, decimal_places=3):
def SaveListTxt(filename: '.txt file', ListToSave: 'A list of lists of results for each starting betagamma'):
def AppendListTxt(filename:'file to be appended to', ListToSave: 'A list of lists of results for each starting betagamma'):
"""
def testfunc(i):
    print(i)
    return i
def distance(i:'first qubit',j:'second qubit',n_qubits,gridWidth:'use a grid width to rank the qubits'=0):
    """
    Return a distance based on the numbers of two different qubits
    """
    if gridWidth ==0:
        #Calc the distance^2 between node i and node j. ASSUME, square gird with n nodes, and distance of 1 between neighbouring nodes
        gridWidth = np.sqrt(n_qubits)
    y = (i//gridWidth - j//gridWidth)**2
    x = (i%gridWidth - j%gridWidth)**2
    return (x+y)                    #distance squared
"""def print_details(quboUnordered,n_qubits,name):
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
"""
#'type:int number of qubits in the system '

#from pyquil.paulis import *




def OpenListTxt(filename, decimal_places=3):
    """
    filename is presumed to be made up of a list containing numbers and lists
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
