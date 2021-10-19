# draw_rigetti.py

"""
def gate_update(gate_type= 'q',new_depth=0,gate_param=0,gate_cntrl=0,qubit_label=0,gate_fontsize = 9,wire_gap= 0.2,gate_spacing= 0.05)
def instruction_generator(prog):
def is_spanned(gate_kwargs_compare,gate_CNOT):
def position_CNOT(circuit_matrix,gate_kwargs_latest):
def reposition_old_gates(circuit_matrix,gate_kwargs_latest):  
def is_proposed_scanned_by_CNOT(gate_kwargs_latest):
    Return Is an existing CNOT gate, spanning this proposed gate ?
    type bool

    param gate_kwargs_latest is the proposed gate, to be added to the circuit_matrix
      
def draw(prog)
"""  
from matplotlib import axes, figure
from pyquil.paulis import *  #  exponential_map, exponentiate_commuting_pauli_sum
from pyquil.api import *
from pyquil.gates import *
from pyquil.gates import Z,CNOT
from pyquil.api import local_forest_runtime, WavefunctionSimulator
import matplotlib.pyplot as plt
# https://matplotlib.org/stable/api/_as_gen/atplotlib.axes.Axes.annotate.html#matplotlib.axes.Axes.annotate
import numpy as np
# globals used

fig = figure 
ax = axes 
circuit_matrix = []
circuit_width = 32
circuit_height = 7
wire_gap= circuit_height*0.012
gate_spacing= 0.05
gate_fontsize = 9

idx_frm_label = {}
gate_depth_by_qubit_idx = {}

def draw(prog, first_gate_to_draw=0,max_gates_to_draw=30):
    del circuit_matrix[:]
    idx_frm_label.clear()
    gate_depth_by_qubit_idx.clear()
    #fig1, ax1 = plt.subplots(figsize=(circuit_width, circuit_height))
    #fig=fig1
    #ax=ax1
    plt.plot([0],[0],'w')
    plt.plot([1],[1],'w')
    #my_gate_spacing= 0.03
    ############## Add and label a wire for each qubit ##############

    for idx,q in enumerate(prog.get_qubits()):
        idx_frm_label.update({q:idx})      #get_index
        gate_kwargs ={'gate_type':'q','gate_param':0,'gate_cntrl':0,'qubit_label':q,'modifiers':''} 
        gate_update(**gate_kwargs)
        circuit_matrix.append([gate_kwargs])
    #gate_update(['RX',.3,3],wire_gap=wire_gap_my)
    g= instruction_generator(prog)
    ############## Add gates to memory matrix ############## 
    #Iterate over prog
    for n,gate_kwargs_latest in enumerate(g):
        if n >= first_gate_to_draw and n< first_gate_to_draw+max_gates_to_draw:
            if gate_kwargs_latest['gate_type'][0]=='C' or gate_kwargs_latest['modifiers']=='C':
                #### find depth to position this CNOT
                deepest_spanned_CNOT_depth = position_CNOT(circuit_matrix,gate_kwargs_latest)
                add_depth_CNOT_cntrl_latest(deepest_spanned_CNOT_depth,circuit_matrix,gate_kwargs_latest)
                # Reposition any gates (already in the circuit_matrix) that are overwritten by a cnot gate connection
                #if gate_kwargs_latest=={'gate_type': 'CNOT', 'gate_param': 0, 'gate_cntrl': 0, 'qubit_label': 4}: print('CNOT(0,4) positioned')
                reposition_old_gates(circuit_matrix,gate_kwargs_latest)                
            else:
                """ check there isnt an existing CNOT gate spanning this proposed gate"""
                
                existing_wire = circuit_matrix[idx_frm_label[gate_kwargs_latest['qubit_label'] ]  ]
                #if gate_kwargs_latest['gate_type']=='Z':
                #    print(is_proposed_scanned_by_CNOT(gate_kwargs_latest), ' = cnot at the proposed depth.')
                if is_proposed_scanned_by_CNOT(gate_kwargs_latest):
                    existing_wire.append( 0)
                
                existing_wire.append( gate_kwargs_latest)
                #if gate_kwargs_latest['qubit_label'] ==0: print(len(circuit_matrix[0]), ' q_0 len()')
            if n > 11 and 0: print(gate_kwargs_latest)
    ############## DRAW gates from matrix ############## 
    #print('DRAW gates from matrix ')
    len_deepest_wire =0
    for wire in range(len(circuit_matrix)):   
        len_deepest_wire = max(len_deepest_wire, len(circuit_matrix[wire]) )
    for depth  in range(len_deepest_wire):
        for wire in range(len(circuit_matrix)):
            if  len(circuit_matrix[wire]) > depth:  
                gate_kwargs = circuit_matrix[wire][depth]
                if type(gate_kwargs) == dict:
                    gate_update(**gate_kwargs,new_depth=depth )
    #ax.axis(False)                    # hide axies
    plt.show()

def gate_update(gate_type= 'q',new_depth=0,gate_param=0,gate_cntrl=0,qubit_label=0,modifiers=''):
    """
    Create the axis command to draw annotations to represent a series of gate.
    gate_kwargs ={'gate_type':'q','gate_param':0,'gate_cntrl':0,'qubit_label':q,modifiers=''}
  
    """
    qubit_idx = idx_frm_label[qubit_label]
    cntrl_idx =  idx_frm_label[gate_cntrl]
    #{'gate_type':0,'gate_param':0,'gate_cntrl':0,'qubit_label':0}
    if not gate_depth_by_qubit_idx.__contains__(qubit_idx): 
        gate_depth_by_qubit_idx.update({qubit_idx:0})
    
    gate_color={
                        'id': '#F0E442',
                        'u0': '#E7AB3B',
                        'u1': '#E7AB3B',
                        'u2': '#E7AB3B',
                        'u3': '#E7AB3B',
                        'X': '#58C698',
                        'Y': '#58C698',
                        'Z': '#58C698',
                        'H': '#70B7EB',
                        's': '#E0722D',
                        'sdg': '#E0722D',
                        't': '#E0722D',
                        'tdg': '#E0722D',
                        'RX': '#b342f5',
                        'rx': '#b342f5',
                        'RY': '#b342f5',
                        'RZ': '#b342f5',
                        'CZ': '#b342f5',
                        'CNOT':'lightblue',
                        'reset': '#D188B4',
                        'target': '#70B7EB',
                        'meas': '#D188B4',
                        'q' :'w'
                    }
    gate_color_used = 'yellow'      #default color of box representing a gate
    if list(gate_color.keys()).count(gate_type)>0: gate_color_used = gate_color[gate_type]
    gate_string= {'q':f'q_{int(qubit_label):d}',
                                        'RX':r'$R_x$'+' \n'+f'{gate_param:3.1f}',
                                        'RY':r'$R_y$'+ '\n'+f'{gate_param:3.1f}',
                                        'RZ':r'$R_z$'+ '\n'+f'{gate_param:3.1f}',
                                        'PHASE':'Pha'+'\n'+f'{gate_param:3.1f}',
                                        'CNOT':'+',
                                        'CZ':r'CZ'+'\n'+'  ',
                                         }                   
        
    if list(gate_string.keys()).count(gate_type)==0:
        gate_string_used = '  '+gate_type +'  \n'
    else:
        gate_string_used = gate_string[gate_type]
    
    gate_pos = (gate_spacing*(new_depth), 1-(qubit_idx)*wire_gap)  
     
    
    if gate_type[0] == 'C' or modifiers=='C':
          
        ###### arrows cntl to 1) CNOT 2) prior gate 
        cntrl_pos =   (gate_spacing*(new_depth), 1-(cntrl_idx)*wire_gap)
        
        prior_gate_depth = gate_depth_by_qubit_idx[cntrl_idx]
        #arrow_to_prior_gate = (gate_spacing*(3),cntrl_pos[1])
        arrow_to_prior_gate = (gate_spacing*(prior_gate_depth),cntrl_pos[1])
        #if cntrl_idx==1: print(gate_depth_by_qubit_idx[1],'gate_depth_by_qubit_idx[1]. pos = ', new_depth, arrow_to_prior_gate, ' =arrow_to_prior_gate')
        #print(arrow_to_prior_gate , cntrl_idx)
        #arrow_to_prior_gate = (0,cntrl_pos[1])
        if gate_type == 'CNOT':
            plt.plot(cntrl_pos[0],cntrl_pos[1],'o',color='lightblue')                   #plot CNOT cntrl 
            #1) cntrl to CNOT
            my_cntrl = ax.annotate(gate_string_used, xy=cntrl_pos, xycoords="data",     #arrow CNOT target to cntrl 
                        xytext=gate_pos, textcoords="data",                             #CNOT target
                        arrowprops=dict(arrowstyle="-",color='lightblue'),
                        va="center", ha="center",
                        bbox=dict(boxstyle='circle',edgecolor='w', fc=gate_color_used)#,color='w')
                        ) 
            my_cntrl.set_color('w') #color of text in box
        else:
            plt.plot(cntrl_pos[0],cntrl_pos[1],'o',color='lightblue')                   #plot CZ cntrl 
            #1) cntrl to CNOT
            my_cntrl = ax.annotate(gate_string_used, xy=cntrl_pos, xycoords="data",     #arrow CZ target to cntrl 
                        xytext=gate_pos, textcoords="data",                             #CZ target
                        arrowprops=dict(arrowstyle="-",color='lightblue'),
                        va="center", ha="center",
                        bbox=dict(boxstyle='square',edgecolor='w', fc=gate_color_used)#,color='w')
                        ) 
            my_cntrl.set_color('b') #color of text in box
        #2) cntrl to prior gate
       
        arrow_cntrl_prior_gate= ax.annotate('', xy=arrow_to_prior_gate, xycoords="data", color='y',    #arrow cntrl to prior gate 
                    xytext=cntrl_pos, textcoords="data",                             #CNOT target
                    arrowprops=dict(arrowstyle="-",color='lightblue'),
                    va="center", ha="right",
                    bbox=dict(boxstyle='circle',edgecolor='w', fc=gate_color_used)#,color='w')
                    )
        
        prior_gate_depth = gate_depth_by_qubit_idx[qubit_idx]
        arrow_to_prior_gate = (gate_spacing*(prior_gate_depth),gate_pos[1])

        ###### arrow: CNOT to arrow_to_prior_gate
        my_gate = ax.annotate('', xy=arrow_to_prior_gate, xycoords="data",    #arrow to arrow_to_prior_gate
                    xytext=gate_pos, textcoords="data",                             #xy positions bbox relative to axies
                    arrowprops=dict(arrowstyle="-",color='lightblue'),
                    va="center", ha="right",
                    bbox=dict(boxstyle='circle',edgecolor='w', fc=gate_color_used)#,color='w')
                    )  
        my_gate.set_color('w') #color of text in box
        
    else:
        prior_gate_depth = gate_depth_by_qubit_idx[qubit_idx]
        arrow_to_prior_gate = (gate_spacing*(prior_gate_depth),gate_pos[1])
      
        my_gate = ax.annotate(gate_string_used, xy=arrow_to_prior_gate, xycoords="data", #xy positions bbox relative to axises
                    xytext=gate_pos, textcoords="data",
                    arrowprops=dict(arrowstyle="-",color='lightblue'),
                    va="center", ha="right",        #left
                    bbox=dict(boxstyle='square',edgecolor='w', fc=gate_color_used)
                    )  
    #gate_depth_by_qubit_idx[qubit_idx]+=1
    gate_depth_by_qubit_idx[qubit_idx] = new_depth
    my_gate.set_fontsize(gate_fontsize) 
    gate_string['new_gate']=0
    #return   gate_depth+1
#ax.axis(False)                    # hide axies
################# TEST gate_update()

def test_gate_update(prog):
    
    #‘data’	use the coordinate system of the object being annotated (default)
    plt.plot([0],[0],'o')
    plt.plot([1],[1],'o')
    idx_frm_label = {}
    

    wire_gap_my = 0.07
    for idx,q in enumerate(prog.get_qubits()):
        idx_frm_label.update({q:idx})
        gate_kwargs ={'gate_type':'q','gate_param':0,'gate_cntrl':0,'qubit_label':q,'modifiers':''}
        gate_update(**gate_kwargs,wire_gap=wire_gap_my)
        
    gate_update(gate_type='H',qubit_label=3,wire_gap=wire_gap_my)
    gate_update(gate_type='H',qubit_label=3,wire_gap=wire_gap_my)
    gate_update(gate_type='H',qubit_label=3,wire_gap=wire_gap_my)
    gate_update(gate_type='H',qubit_label=3,wire_gap=wire_gap_my)
    gate_update(gate_type='H',qubit_label=0,wire_gap=wire_gap_my)
    gate_update(gate_type='H',qubit_label=0,wire_gap=wire_gap_my)
    gate_update(gate_type='CNOT',qubit_label=0,gate_cntrl=2,wire_gap=wire_gap_my)
    gate_update(gate_type='H',qubit_label=0,wire_gap=wire_gap_my)


def instruction_generator(prog):
    """
    Yield, as kwargs for update(**kwargs), the sequential gate instructions that make up the program. Used for subsequent drawing of each gate.
    """
    for instr in prog:  
            #parameter, eg angle
            param = 0
            cntrl = 0
            mod=''
            if len(instr.modifiers) >0: mod='C'# different if it is double controlled
            if len(instr.params)>0:
                param = round(instr.params[0],1)            
            # cntrl  eg for  a CNOT or CZ gate
            for n, q in enumerate(instr.qubits):
                if n==0 and len(instr.qubits)>1:
                    #print(q,'=cntl',end=',')
                    cntrl = int(str(q))
            # qubit label number
            index = int(str(q))   #derive the number of a)what is controlled by the qubit,b)the qubit, and c) the gate description
            #print(index, ' index')  #debug
            yield{'gate_type':instr.name,'gate_param':param,'gate_cntrl':cntrl,'qubit_label':index,'modifiers':mod} #instr.modifiers[0]='CONTROLLED'
            #yield [instr.name,param_or_cntrl,int(str(q))]
def is_spanned(gate_kwargs_compare,gate_CNOT):
    """
    Return True if the latest gate cntrl to label has wires in common with either
    a) a CNOT gate to be compared to OR
    b) a simple gate to be compared
    """
    if type(gate_kwargs_compare) != dict: return False
    if gate_kwargs_compare['gate_type'][0] == 'C':
        min_compare = min(idx_frm_label[gate_kwargs_compare['qubit_label'] ], idx_frm_label[gate_kwargs_compare['gate_cntrl']] )
        max_compare = max(idx_frm_label[gate_kwargs_compare['qubit_label'] ], idx_frm_label[gate_kwargs_compare['gate_cntrl']] )
    else:
        min_compare = idx_frm_label[gate_kwargs_compare['qubit_label'] ]
        max_compare = min_compare
    min_latest = min(idx_frm_label[gate_CNOT['qubit_label'] ], idx_frm_label[gate_CNOT['gate_cntrl']] )
    max_latest = max(idx_frm_label[gate_CNOT['qubit_label'] ], idx_frm_label[gate_CNOT['gate_cntrl']] )
    #print(min_compare,max_compare)
    #print(min_latest,max_latest)
    for q_compare in range(min_compare,max_compare+1):
        for q_latest in range(min_latest,max_latest+1):
            if q_compare==q_latest: return True
    return False

def position_CNOT(circuit_matrix,gate_kwargs_latest):
    """
    Return the instruction or depth? of the CNOT already psotioned in the circuit_matrix that is both 
        a) SPANNED by the latest CNOT, and 
        b) is the DEEPEST example of this.
    """
    deepest_spanned_CNOT_depth = 0
    for wire in range(len(circuit_matrix)):   
        #idx_frm_label[n]             is the qubit number and n is the rank of that qubit
        len_wire = len(circuit_matrix[wire]) 
        for idx,gate_kwargs_in_matrix in enumerate(reversed(circuit_matrix[wire])):
            #if wire==3: print(depth,gate_kwargs['gate_type'])
            depth = len_wire-1-idx
            #print(depth_max, wire, ' depth_max, n')
            # select CNOT on current wire, deepest first. Test whether it overlaps with the CNOT gate to be position, the latest
            if type( gate_kwargs_in_matrix) == dict:
                if gate_kwargs_in_matrix['gate_type'][0]=='C':
                    if is_spanned(gate_kwargs_in_matrix,gate_kwargs_latest) :
                        deepest_spanned_CNOT_depth = max(deepest_spanned_CNOT_depth,depth)
                        #print(deepest_spanned_CNOT_depth)
                        #print(gate_kwargs_in_matrix,gate_kwargs_latest, depth, wire, ' = n')
    return deepest_spanned_CNOT_depth

def add_depth_CNOT_cntrl_latest(deepest_spanned_CNOT_depth,circuit_matrix,gate_kwargs_latest):
    # position of qubit_label
    wire_for_cnot= idx_frm_label[gate_kwargs_latest['qubit_label'] ] 
    depth_curr_latest = len(circuit_matrix[wire_for_cnot])
    # position gate_cntrl
    wire_for_cntrl = idx_frm_label[gate_kwargs_latest['gate_cntrl'] ] 
    depth_curr_cntrl = len(circuit_matrix[wire_for_cntrl])
    # propose position for latest CNOT gate
    depth_propose = max(depth_curr_cntrl, depth_curr_latest,deepest_spanned_CNOT_depth+1) 
    # ADD to position of qubit_label
    for idx in range(depth_curr_latest,depth_propose):
            circuit_matrix[wire_for_cnot].append( 0)
    circuit_matrix[wire_for_cnot].append( gate_kwargs_latest)
    # position gate_cntrl
    for idx in range(depth_curr_cntrl,depth_propose):
        circuit_matrix[wire_for_cntrl].append( 0)
    circuit_matrix[wire_for_cntrl].append( 0)

def reposition_old_gates(circuit_matrix,gate_kwargs_latest):
    """
    When positioning a CNOT gate, test if there are gates spanned by the cntrl to qubit of the CNOT and shuffle them to the Left on their wire
    """
    qubit_idx= idx_frm_label[gate_kwargs_latest['qubit_label']]
    depth_latest = len (circuit_matrix[qubit_idx])-1

    for q in range(len(circuit_matrix)):# range(min_wire,max_wire):
        if  q !=qubit_idx \
            and len(circuit_matrix[q]) > depth_latest:       # 1>2-1?
            #print(q ,'=q is deep enough',len(circuit_matrix[q]), '= len(circuit_matrix[q])', depth_latest, '=depth_latest')       #0,2
            gate_kwargs_compare = circuit_matrix[q][depth_latest]
            if type(gate_kwargs_compare) == dict:
                if is_spanned(gate_kwargs_compare,gate_kwargs_latest):
                    #print(q ,'=q is Spanned. with gate of ',gate_kwargs_compare['gate_type'] ,len(circuit_matrix[q]), '= len(circuit_matrix[q])', depth_latest, '=depth_latest')     
                    circuit_matrix[q].insert(depth_latest,0)
        #else:
         #   print(q ,'=q is NOT deep enough',len(circuit_matrix[q]), '= len(circuit_matrix[q])', depth_latest, '=depth_latest')      
          #  if q==1: print(q ,'=q',circuit_matrix[q])
def is_proposed_scanned_by_CNOT(gate_kwargs_latest):
    """
    Return Is an existing CNOT gate, spanning this proposed gate ?
    type bool

    param gate_kwargs_latest is the proposed gate, to be added to the circuit_matrix
    """
    
    existing_wire = circuit_matrix[idx_frm_label[gate_kwargs_latest['qubit_label'] ]  ]
    depth_to_search=len(existing_wire)
    for q in range(len(circuit_matrix)):
        if len(circuit_matrix[q]) > depth_to_search:                
            if type(circuit_matrix[q][depth_to_search]) == dict:
                existing_gate = circuit_matrix[q][depth_to_search]
                if  existing_gate['gate_type'][0] =='C':                #is there a 'C' gate, ie a control gate eg CZ or CNOT, at this depth? 
                    #print(existing_gate, ' = existing_gate')
                    is_spanned(gate_kwargs_latest,existing_gate)        #does the CNOT span gate_kwargs_latest
                    return True
    return False



#### Test Draw

def test_Draw():
    if 1:
        prog = Program()
        if 0:
            prog += Z(0)
        prog += Z(0)
        prog +=CNOT(1,2)
        prog += Z(1)
        
        
        #prog += CNOT(0,4)
        #prog += Z(0)
        #prog += Z(4)
        
        #prog += Y(0)
        
        #prog += Z(1)
        #prog += Z(1)
        #prog += Z(2)
    
        draw(prog,max_gates_to_draw=20)
    
        plt.show()