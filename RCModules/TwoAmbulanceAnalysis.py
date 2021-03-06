#TwoAmbulanceAnalysis.py
#Methods used in the cell below
# These methods were first written to display Dwave solutions
import matplotlib.pyplot as plt
from datetime import datetime

def qubo_energy_value(ψ, qubo):
    E=0
    for qubo_pos, weight in qubo.items():   #qubo_pos(q1,q2)
        E += ψ[qubo_pos[0]] * ψ[qubo_pos[1]] * weight                               #Calc energy of the Ising nodes with weights contained in linear
    return E
    
def Grid_Destinations(psi_opt, H ,W):
    """
    Print 1 to indicate the position of AO destinations, or 0 for A1 destinations.
    """
    idx = 0
    Grid_Locations(psi_opt, H ,W, idx)
def Prt_A0_A1_Destinations(psi_opt,Height,Width):
    """Print Ambulances destinations
    param H 'grid height (Type:int)'

    param W 'grid width (Type:int)'

    param psi_opt 'Potential binary solution to objective function (Type:list int)'

    Return (Type:bool) True if the ambulances do not share a destination, otherwise false
    """
    ConstraintMet = 1
    for q0, q1 in zip(psi_opt[:Width * Height],  psi_opt[Width * Height:2*Width * Height]):
        if q0 + q1 >1:
            ConstraintMet = 0
        
    if ConstraintMet:
        print('Destination Constraints were Met.','\n Map of destinations of each ambulance, 1 for A0, 0 for A1')
        Grid_Destinations(psi_opt,Height,Width)
    else:
        print('ERROR: Destination Constraints not Met ')
        print(psi_opt[:Width * Height],              '= psi_opt first ambulance destinations')
        print(psi_opt[Width * Height:2*Width * Height], '= psi_opt 2nd ambulance destinations \n')
    return ConstraintMet

def Grid_StartPositionsAO(psi_opt, H ,W):
    """
    param H 'grid height (Type:int)'

    param W 'grid width (Type:int)'

    param psi_opt 'Potential binary solution to objective function (Type:list int)'

    Output prints starting positions of A0  in a grid measuring H by W
    """
    idx = 2
    print('StartPositionsA0')
    NLocs, Apos = Grid_Locations(psi_opt, H ,W, idx)
    if  NLocs==1 :
        print('\t\t\t Start Position A0 Constraints met')
    else:
        print('\t\t\tERROR Start Position A0 Constraints not met')
    return Apos

def Grid_StartPositionsA1(psi_opt, H ,W):
    """
    param H 'grid height (Type:int)'

    param W 'grid width (Type:int)'

    param psi_opt 'Potential binary solution to objective function (Type:list int)'

    Output prints starting positions of A1  in a grid measuring H by W

    Return qubit number representing the starting psoition of A0 (Type:int)
    """
    idx = 3
    print('StartPositionsA1')
    NLocs, Apos = Grid_Locations(psi_opt, H ,W, idx)
    if  NLocs==1 :
        print(' \t\t\tStart Position A1 Constraints met')
    else:
        print(' \t\t\tERROR Start Position A1 Constraints not met')
    return Apos
def Grid_Locations(psi_opt, H ,W,idx:'Feature, eg AO_Dest or A)_Start'):
    """
    Print the rows of a H*W grid that refers to feature number idx, 
    eg idx=0 is the grid of A0_Destination
    eg idx=1 is the grid of A1_Destination
    eg idx=2 is the grid of A0_Start
    """
    i=0
    sumlocs=0
    if len(psi_opt) >= (idx+1)* W*H:  
        StartGrid =  int(idx * W*H )         
        for i in range(0,int(H)):
            s = int( W * i + idx * W*H )        #Start of row
            e = int( s + W  )                   #End of row
            #Print each row, i, of the grid feature idx
            print (psi_opt[s:e])    
            
            sumlocs += sum(psi_opt[s:e])
        if sumlocs > 0:
            Apos = psi_opt[StartGrid:e].index(1)
        else:
            Apos = -1
    else:
        print('ERROR: psi_opt not large enough to print grid')
    return sumlocs,Apos 

def SolutionDistance(psi_opt,A0_sol,A1_sol,qubo,Nlocs):
    """
    params:
    psi_opt (Type:list int, length 4*Nlocs)
    A0_sol  (Type:int) value is 0 <= A0_sol <= Nlocs
    A1_sol  (Type:int) value is 0 <= A0_sol <= Nlocs
    qubo    (Type:dict (x,y):float)   When x!=y this is the edge weight between qubit_x and qubit_y
    Nlocs   (Type:int) Number of possible locations to be served by the ambulances.

    """
    d0 = 0
    d1 = 0

    for  qubo_pos, elem  in qubo.items():   #qubo_pos(y,x)
        if qubo_pos[1] == A0_sol + 2 * Nlocs:
            if qubo_pos[0] < 2*Nlocs :
                d0 += elem * psi_opt[ qubo_pos[0] ] 
        if qubo_pos[1] == A1_sol + 3 * Nlocs:
            if qubo_pos[0] < 2*Nlocs :
                d1 += elem * psi_opt[ qubo_pos[0] ]
    return d0 + d1
    
def Print_Destinations(response, Height,Width, sampler_used='Sampler()', A0_sol=-1 , A1_sol=-1, show_start_pos=0):
    print(response.first[1], ' = Lowest energy found by',sampler_used,' in a grid',Width,'(w) by',Height,'(h)' )
    ConstraintMet = 1
    psi_opt =list( response.first[0].values())
    display_psi = psi_opt
    #Check that Ambulances do not share a destination
    for q0, q1 in zip(psi_opt[:Width * Height],  psi_opt[Width * Height:2*Width * Height]):
        if q0 + q1 >1:
            ConstraintMet =0
    #Add 10 or 20 to indicate the starting position of A0 and A1
    if A0_sol != -1 and show_start_pos:
            display_psi[A0_sol] = display_psi[A0_sol] + 10
    if A1_sol != -1 and show_start_pos:
        display_psi[A1_sol] = display_psi[A1_sol] + 20
        
    if ConstraintMet:
        print('Destination Constraints were Met.')
        print(' Map of destinations of each ambulance, 1 for A0, 0 for A1. Where A0 starts add 10, A1 add 20')
        Grid_Destinations(display_psi,Height,Width)
    else:
        print('ERROR: Destination Constraints not Met ')
        print(display_psi[:Width * Height],              '= psi_opt first ambulance destinations')
        print(display_psi[Width * Height:2*Width * Height], '= psi_opt 2nd ambulance destinations \n')
    return psi_opt