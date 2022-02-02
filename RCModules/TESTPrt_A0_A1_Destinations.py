#TESTPrt_A0_A1_Destinations
from TwoAmbulanceAnalysis import *
if 1:
    psi_opt= [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0]
    Height=7 
    Width =7
    UseHybrid = 0
    if Height==7 and Width==7:      #I need to check this is the optimum position of the ambulances
        A0_sol = 22
        A1_sol = 26  
if 0:
    Distance = SolutionDistance(psi_opt,A0_sol,A1_sol,Adjacency,Nlocs) #(27, 82, 121, 1276) Asymmetric orthogonal 
    print(Distance,' = Distance')
if 0 :
    Prt_A0_A1_Destinations(psi_opt,Height,Width) 
    print('constraint met')
if 0:Print_Destinations(response, Height,Width,Method, A0_sol, A1_sol,show_start_pos=1)


if 1:
     Prt_A0_A1_Destinations(psi_opt,Height,Width,A0_sol=A0_sol , A1_sol=A1_sol,show_start_pos=1)

idx = 0 # 'Feature, eg AO_Dest or A)_Start'
print(len(psi_opt))
if not UseHybrid and 1:
    Grid_StartPositionsA0(psi_opt,Height,Width, prt_locs=1)
    Grid_StartPositionsA1(psi_opt,Height,Width, prt_locs=1)
#Grid_Locations(psi_opt,Height,Width,idx)
print(psi_opt[:Width * Height],              '= psi_opt first ambulance destinations')
print(psi_opt[Width * Height:2*Width * Height], '= psi_opt 2nd ambulance destinations \n')
A0_sol,A1_sol,Height,Width,
