# PlatformConversion.py
#All functions in this file should not include any quantum computer specific code, eg not import qiskit nor import pyquil, nor import dimod
def test1(x):
    if x:
        print(x,' != 0 hi')
    else:
        print(x,' == 0')
"""
Adjaceny_to_DWave_ising( qubo: 'Adjaceny table in dict form: qubo[(r,c)]', n_qubits):
def SaveListTxt(filename: '.txt file', ListToSave: 'A list of lists of results for each starting betagamma'):
def SaveList_as_txt_file(ListToSave, filename):
def CreateTwoAmbulanceAdjacencyV2(gridWidth,n_destinations: 'Number of locations servied by ambulance fleet'= 5, Adddistance = 1
,ConstraintMultiplier:'If ConstraintMultiplier= 0, it is calculated automatically'=0, use_XYMixer_constraints = 0 ):
def CreateTwoAmbulanceAdjacencyV1(gridWidth,n_destinations: 'Number of locations servied by ambulance fleet'= 5, Adddistance = 1
,ConstraintMultiplier:'If ConstraintMultiplier= 0, it is calculated automatically'=0,  qubo_model:'qubo_model definition used in Adjacency table'=True):
def CreateTwoAmbulanceAdjacency(gridWidth,n_destinations: 'Number of locations servied by ambulance fleet'= 5, Adddistance = 1
,ConstraintMultiplier:'If ConstraintMultiplier= 0, it is calculated automatically'=0,  qubo_model:'qubo_model definition used in Adjacency table'=True):
def distance(r1,c1,r2,c2):
def co_ord(qubit,Height, Width,Feature_options):
"""
def Adjaceny_to_DWave_ising( qubo: 'Adjaceny table in dict form: qubo[(r,c)]', n_qubits):
    """
    Returns: linear, quadratic (type: dicts) that are required by the DWave function BinaryQuadraticModel.from_ising(linear, quadratic)
    DWave Ising model needs: linear, quadratic;
    eg:     linear =  {0: -1, 1: -1, 2: -1}
    and:    quadratic = {(0, 1): 1, (0, 2): 1, (1, 2): 1}
    DWave qubo model needs: qubo            eg qubo = {(0, 0): -1, (0, 1): 1, (0, 2): 1, (1, 1): -1, (1, 2): 1, (2, 2): -1}
    """
    #Create parameters for the Ising model in Dwave
    linear = {0:1.0}#   Node strength
    quadratic = {(0,1):0}
    for r in range(0, n_qubits):
        for c in range(0, n_qubits):
            if c == r:
                linear[c] = qubo[(r,c)] 
            if c > r :
                quadratic[(r,c)]        = qubo[(r,c)] 
    return linear, quadratic

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
def SaveList_as_txt_file(ListToSave, filename):
    """
    Create a line of strings populated by the sequential elements of ListToSave,  separated by '\t' . Store these in a .txt file
    eg filename = 'myfile.txt'
    """
    f = open(filename, "w+")    #erases files if it already exists
    for d in range(len(ListToSave)):
        for n, prob in enumerate(ListToSave[d] ):
            f.writelines(str(prob)+'\t')
        f.writelines('\n')              #marks the end of a line
    f.close()

def print_QUBOdetails(quboUnordered: 'dict of edge and nodes', n_qubits, filename:'name of .txt file to create'):
    """
    Print adjacency MATRIX of weights Wij between one qubit index (Qi)  and another (Qj)
    Save this as a text file with filename.
    eg if quboUnordered = {(0, 0): -160, (0, 5): 320, (0, 10): 1} that is a node of -160 at Q0, a 320 edge at Q0 to Q5, and a 1 edge Q0 to Q10
    """
    r = 0                                   # index of row of adjacency matrix
    c = 0                                   # index of column of adjacency matrix
    print( 'ADJACENCY MATRIX',filename,' \n\n', end='')
    qubo = {}
    qubolist = [[]]
    line = 0
    #order the dicts by the first element of each tuple (its row), this helps the printing order
    for key in sorted( quboUnordered.keys() ):
        qubo[key] = quboUnordered[key]
    qubolist[line].append('Qubit')
    for n in range(0,n_qubits):
        qubolist[line].append('q{}'.format(n))
    # FIRST row complete
    # Start list of second line
    qubolist.append([])
    line +=1
    c=-1
    for pos, weight in qubo.items():        # 'weight' is the weight between qubit r and qubit c pos(r,c) 
        if pos[0] > r :                     #row
            c=-1    
        #1) New list for each blank line    
        for m in range(r,pos[0]):
            qubolist.append([])
            line +=1
        r = pos[0]
        #2) first element of a line
        if c==-1:       
            qubolist[line].append('q{}'.format(r))
            qubolist[line].append(' ')
            c=0
        #3) Rest of line
        for n in range(c,pos[1]):     # (0,0) c=0
            qubolist[line].append(' ')     
        c = pos[1]
        #print(weight, end='') 
        qubolist[line].pop()
        qubolist[line].append(weight)
    SaveList_as_txt_file(qubolist, filename)
    #print out the qubo
    maxlen =   0
    # find the longest string... 
    for linelist in qubolist:
        for s in linelist:
            maxlen = max(maxlen, len(str(s)))
    #...then format the column width to this
    for linelist in qubolist:
        for s in linelist:
            print(str(s).center(maxlen), end='')
        print('\n')
 
   
def CreateAmbulanceAdjacency(gridWidth,n_qubits = 5, ConstraintMultiplier = 1, Adddistance = 1,remove_constraint =True,
HammingWeightOfConstraint=1, qubo_model=True ):
    """
    
    Returns a 2d dict: type {(0,0):0}. n_qubits * n_qubits,  of an adjacency table for the ambulance problem. 
     
    HammingWeightOfConstraint = 2 sets as a constraint the number of 1s in the solution to 2.
    remove_constraint: type bool = True , removes the constraint(number of ambulances) from the adjacency table.
    Adddistance: type bool = True, removes the distance between nodes from the adjacency table
    qubo_model:type bool = True, when True the adjacency is based on a qubo model of edges and nodes, otherwise Ising model
    """
    
    # CONSTRAINTS are calculated differently for a binary/qubo model (x= 0 or 1), than for a spin/Ising model (x= -1 or 1)
    if qubo_model :       #QUBO
        linearWt = 1-2* HammingWeightOfConstraint               # this is the constraint condition for a binary/qubo model (x= 0 or 1)
        quadWt = 2
    else:       #ISING
        linearWt = n_qubits-2* HammingWeightOfConstraint          # this is the constraint condition for a spin/Ising model (x= -1 or 1)
        quadWt = 1                                              

    Adjacency   = {(0,0):0}
    linearWt = linearWt * ConstraintMultiplier
    for r in range(0,n_qubits):
        for c in range(0,n_qubits):
            if c == r:
                if remove_constraint==0:
                    Adjacency[(r,c)] = linearWt
            elif c > r :
                if Adddistance :
                        coeff = - distance(r,c, n_qubits,gridWidth)  # distance assumes all qubits are evenly space in a grid structure
                        if remove_constraint==0:
                            coeff +=  ConstraintMultiplier*quadWt
                else:
                    if remove_constraint==0:
                        coeff = ConstraintMultiplier*quadWt       # use this to calculate solely the constraint weights 
                    else:
                        coeff = 0          
                Adjacency[(r,c)]             = coeff
    return Adjacency

   
    
def CreateTwoAmbulanceAdjacencyV2(gridWidth,n_destinations: 'Number of locations servied by ambulance fleet'= 5, Adddistance = 1
,ConstraintMultiplier:'If ConstraintMultiplier= 0, it is calculated automatically'=0, use_XYMixer_constraints = 0 ):
    """
    Returns: type dict result={'qubo':qubo, 'quboHybrid':quboHybrid, 'n_qubits':n_qubits, 'ConstraintMultiplier':ConstraintMultiplier,'max_distance': max_distance, 'sum_distance':sum_distance}
    Returns: type dict eg  result['qubo']={(0,3):4}, is an edge of weight 4 between qubits 0 and 3.  QUBO to represent the ising ambulance problem, it involves two constraints and minimising distance driven.
    params use_XYMixer_constraints (type:bool) uses fewer constraints see below
################# DESCRIPTION OF THE ISING PROBLEM #################
BINARY/qubo MODEL  
This splits the geography into Width*Height (W*H) qubits, and adds the distance^2 from...
        A0 starting location (described by a W*H qubits labelled 'A0_Start')...
        to destination locations (W*H qubits labelled collectively A0_Des )
        and
        A1 starting location (described by a W*H qubits labelled 'A1_Start')...
        to destination locations (W*H qubits labelled collectively A1_Des )
With the contraints that if use_XYMixer_constraints == False;
       A0 and A1 have only one starting location (ie there is a single one in the W*H qubits collectively labelled 'A0_Start', and 'A1_Start')
       A0 only services the destinations in A0_Des that have a 1
       A0 does not share any destinations with A1
Else   if use_XYMixer_constraints== True:
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
    
    Des_quadWt = 2                  # adds energy if A1_des is same a A0_des
    Des_linearWt = -1               # reduces energy if a destination is removed
    
    
    Start_quadWt = 2                  # adds energy if A1_Start is same a A0_Start
    Start_linearWt = 1-2*num_of_ones_constraint               # this is the constraint condition for a binary/qubo model (x= 0 or 1)
    

    ################# Create parameters for two different binary/qubo model. a) Complete Ising b) Ising where the starting position of each ambulance is fixed and known; 'quboHybrid'. 
    qubo   = {(0,0):0}
    quboHybrid   = {(0,0):0}
    
    # Create adjaceny table of linear weights on the diagonal and edge weight in the top triangle above the diagonal
    ################# FIRST calculate the sum of all the distances between destinations to scale up the energy of constraint energies
    sum_distance = 0
    max_distance = 0
    for q1 in range(0,n_qubits):
        feature1,r1,c1 = co_ord(q1,Height, Width,Feature_options)
        for q2 in range(q1,n_qubits):
            feature2,r2,c2 = co_ord(q2,Height, Width,Feature_options)
            if feature1 == 'A0_Des' and feature2 == 'A0_Start':
                if  (r2 * Width )+ c2 > (r1 * Width) + c1 or 1:
                    sum_distance += distance(r1,c1,r2,c2)
                    max_distance = max(max_distance, distance(r1,c1,r2,c2))
    
    #Use sum of distances to ensure any constraint energy is larger than any distance energy
    #sum_distance = 0
    if ConstraintMultiplier == 0:
        if 0 :
            ConstraintMultiplier = sum_distance 
        else:
            ConstraintMultiplier = max_distance * 10#0      #Lagrange parameter
    ################# SECOND create adjaceny maxtrix
    for q1 in range(0,n_qubits):
        feature1,r1,c1 = co_ord(q1,Height, Width,Feature_options)
        for q2 in range(q1,n_qubits):
            feature2,r2,c2 = co_ord(q2,Height, Width,Feature_options)
            
            #########################   Constraint that NO DESTINATION is visited by both ambulances #########################
            if c1 == c2 and r1==r2 and feature1 == 'A0_Des' and feature2 == 'A1_Des' :
                qubo[(q1,q2)] = ConstraintMultiplier * (Des_quadWt )
                if q1 < n_qubitsHybrid and q2 < n_qubitsHybrid:
                    quboHybrid[(q1,q2)] = qubo [(q1,q2)]
            #########################   Constraint that EVERY DESTINATION is visited by at least one ambulances #########################
            if q1 == q2 and (feature1 == 'A0_Des' or feature1 == 'A1_Des') and use_XYMixer_constraints == False:
                qubo[(q1,q2)] = ConstraintMultiplier * Des_linearWt
                if q1 < n_qubitsHybrid and q2 < n_qubitsHybrid:
                    quboHybrid[(q1,q2)] = qubo [(q1,q2)]
            ##################  Constraint that each ambulance has just one START location #########################  
            if  (feature1 == 'A0_Start' or feature1 == 'A1_Start') and feature1==feature2 and use_XYMixer_constraints == False:            
                if(c1 == c2 and r1 == r2):
                    qubo[(q1,q2)] = ConstraintMultiplier * Start_linearWt
                    if q1 < n_qubitsHybrid and q2 < n_qubitsHybrid:
                        quboHybrid[(q1,q2)] = qubo [(q1,q2)]
                else:
                    qubo[(q1,q2)] = ConstraintMultiplier * Start_quadWt
            ##################  Constraint that ambulance 0 minimises START to destination distance #########################  
            if feature1 == 'A0_Des' and feature2 == 'A0_Start' or feature1 == 'A1_Des' and feature2 == 'A1_Start':
                if  (r2 * Width )+ c2 > (r1 * Width) + c1 or 1:
                    d = 0 
                    if Adddistance:
                        d = distance(r1,c1,r2,c2)
                    if d!=0:
                        qubo[(q1,q2)] = d
                #else:
                 #   qubo[(q1,q2)] = 0 
    result={'qubo':qubo, 'quboHybrid':quboHybrid, 'n_qubits':n_qubits, 'ConstraintMultiplier':ConstraintMultiplier,'max_distance': max_distance, 'sum_distance':sum_distance}
    return result  
def CreateTwoAmbulanceAdjacencyV1(gridWidth,n_destinations: 'Number of locations servied by ambulance fleet'= 5, Adddistance = 1
,ConstraintMultiplier:'If ConstraintMultiplier= 0, it is calculated automatically'=0,  qubo_model:'qubo_model definition used in Adjacency table'=True):
    """
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
    
    # Create adjaceny table of linear weights on the diagonal and edge weight in the top triangle above the diagonal
    ################# FIRST calculate the sum of all the distances between destinations to scale up the energy of constraint energies
    sum_distance = 0
    max_distance = 0
    for q1 in range(0,n_qubits):
        feature1,r1,c1 = co_ord(q1,Height, Width,Feature_options)
        for q2 in range(q1,n_qubits):
            feature2,r2,c2 = co_ord(q2,Height, Width,Feature_options)
            if feature1 == 'A0_Des' and feature2 == 'A0_Start':
                if  (r2 * Width )+ c2 > (r1 * Width) + c1 or 1:
                    sum_distance += distance(r1,c1,r2,c2)
                    max_distance = max(max_distance, distance(r1,c1,r2,c2))
    
    #Use sum of distances to ensure any constraint energy is larger than any distance energy
    #sum_distance = 0
    if ConstraintMultiplier == 0:
        if 0 :
            ConstraintMultiplier = sum_distance 
        else:
            ConstraintMultiplier = max_distance * 10#0      #Lagrange parameter
    ################# SECOND create adjaceny maxtrix
    for q1 in range(0,n_qubits):
        feature1,r1,c1 = co_ord(q1,Height, Width,Feature_options)
        for q2 in range(q1,n_qubits):
            feature2,r2,c2 = co_ord(q2,Height, Width,Feature_options)
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
            if feature1 == 'A0_Des' and feature2 == 'A0_Start' or feature1 == 'A1_Des' and feature2 == 'A1_Start':
                if  (r2 * Width )+ c2 > (r1 * Width) + c1 or 1:
                    d = 0 
                    if Adddistance:
                        d = distance(r1,c1,r2,c2)
                    if d!=0 or 0:
                        qubo[(q1,q2)] = d
                #else:
                 #   qubo[(q1,q2)] = 0 
    result={'qubo':qubo, 'quboHybrid':quboHybrid, 'n_qubits':n_qubits, 'ConstraintMultiplier':ConstraintMultiplier,'max_distance': max_distance, 'sum_distance':sum_distance}
    return result
  
def CreateTwoAmbulanceAdjacency(gridWidth,n_destinations: 'Number of locations servied by ambulance fleet'= 5, Adddistance = 1
,ConstraintMultiplier:'If ConstraintMultiplier= 0, it is calculated automatically'=0,  qubo_model:'qubo_model definition used in Adjacency table'=True):
    """
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
    
    # Create adjaceny table of linear weights on the diagonal and edge weight in the top triangle above the diagonal
    ################# FIRST calculate the sum of all the distances between destinations to scale up the energy of constraint energies
    sum_distance = 0
    max_distance = 0
    for q1 in range(0,n_qubits):
        feature1,r1,c1 = co_ord(q1,Height, Width,Feature_options)
        for q2 in range(q1,n_qubits):
            feature2,r2,c2 = co_ord(q2,Height, Width,Feature_options)
            if feature1 == 'A0_Des' and feature2 == 'A0_Start':
                if  (r2 * Width )+ c2 > (r1 * Width) + c1 or 1:
                    sum_distance += distance(r1,c1,r2,c2)
                    max_distance = max(max_distance, distance(r1,c1,r2,c2))
    
    #Use sum of distances to ensure any constraint energy is larger than any distance energy
    #sum_distance = 0
    if ConstraintMultiplier == 0:
        if 0 :
            ConstraintMultiplier = sum_distance 
        else:
            ConstraintMultiplier = max_distance * 10#0      #Lagrange parameter
    ################# SECOND create adjaceny maxtrix
    for q1 in range(0,n_qubits):
        feature1,r1,c1 = co_ord(q1,Height, Width,Feature_options)
        for q2 in range(q1,n_qubits):
            feature2,r2,c2 = co_ord(q2,Height, Width,Feature_options)
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
            if feature1 == 'A0_Des' and feature2 == 'A0_Start' or feature1 == 'A1_Des' and feature2 == 'A1_Start':
                if  (r2 * Width )+ c2 > (r1 * Width) + c1 or 1:
                    d = 0 
                    if Adddistance:
                        d = distance(r1,c1,r2,c2)
                    if d!=0 :
                        qubo[(q1,q2)] = d
                #else:
                 #   qubo[(q1,q2)] = 0 
    result={'qubo':qubo, 'quboHybrid':quboHybrid, 'n_qubits':n_qubits, 'ConstraintMultiplier':ConstraintMultiplier,'max_distance': max_distance, 'sum_distance':sum_distance}
    return result

    import numpy as np
    #METHODS USED
def distance(r1,c1,r2,c2):
    #Calc the distance between node i and node j. ASSUME, square gird with n nodes, and distance of 1 between neighbouring nodes
    y = (r1-r2)**2
    x = (c1-c2)**2
    return (x+y)                    #distance squared
# return (np.sqrt(x+y))


def co_ord(qubit,Height, Width,Feature_options):
    """return type:int,int  row, column in the position grid of dimension Width * Height from the qubit number , and which of 4 features is stored there
    r,c describes what location the feature_num is attached to.
    """
    feature_num = 0
    feature_num = qubit//( Width * Height)
    r = (qubit//Width) %  Height      # // means division followed by integer truncation
    c= qubit % Width

    return Feature_options[feature_num],r,c
