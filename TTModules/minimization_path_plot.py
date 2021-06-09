#import 2 key modules
import numpy as np
import matplotlib.pyplot as plt

def OpenListTxt_TT1(filename, decimal_places=3):
    """
    Filename is presumed to be headed by a Dict followed by a list containing numbers and lists
    The difference from OpenListTxt_TT is that this function does not return the contents of the Dict in the first line of the programme
    """
    #import ast # this module not used as it does not work with multiple Dicts in the first line of a file
    f = open(filename, "r+")        #does not erase file just reads it.
    DataOut = []
    sub_list = []
    liststr = ''
    string = f.readlines()              #read all the lines in the file
    #print('string[0]', string[0])                               # line by line  and create a list using the space delimiter
    dict_str = string[0]
    #print('dict_str', dict_str)
    #convertedDict = ast.literal_eval(dict_str)
    #print(convertedDict)
    string1 = string[1:]
    #print('string1 =', string1)
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
                    sub_list.append( float(('%3.'+str(decimal_places)+'f')% float(y)) )
                    if x.rfind(']') >=0:                #end of a sub_list
                        listfnd -= 1        # SUB number of ]
                        list_of_new_line.append(sub_list)
                        sub_list = []
                else:
                    list_of_new_line.append( float(('%3.'+str(decimal_places)+'f')% float(x)) )
            else:
                #if n < 3:
                    #print(n, 'list_of_new_line', list_of_new_line)
                DataOut.append( list_of_new_line)
    f.close()
    return  DataOut # convertedDict

def angle_plot_TT(filename, energy_break):
        """
        For an input file assumed to have a Dict as the first line, for a p = 1 type QAOA results file with the end and start angles assumed to be located at
        index position [4] and [5] respectively e.g. [2.45, 3.23], [6.234, 0.874], and the EV located at index position [0] in the results file:
        plots 3 graphs, showing
        1) all angle run trajectories colour coded according to the list energy breaks
        2) only the lowest band of energy trajectories, ie the best performing band
        3) only the second lowest band of energy trajectories, ie the second best performing band
        Filename is the name of the .txt results file
        energy_break is a list of 4 numbers specifying the breaks in the data to be used, e.g. [-440, -400, -300, -200]
        This results in there being five bands which are colour coded on the main plot
        Axis ranges are chosen by matplotlib rather than specified by this function
        """
    
        #Read in the data to a list DataOutList
        DataOutList = OpenListTxt_TT1(filename, decimal_places=5)
        
        #Sort the data by EV, assumed to be the first entry on each line of results data
        def takeFirst(elem):
            return elem[0]
        if 1:
            DataOutList.sort(key=takeFirst)
        
        # Calculate the line trajectories and store the data
        X_many = []
        Y_many = []
        Value_many = []
        End_point = []
        Start_point = []
        pi = np.pi
        for i in range(len(DataOutList)):
            temp_x = [DataOutList[i][5][0]/pi, DataOutList[i][4][0]/pi]
            temp_y = [DataOutList[i][5][1]/pi, DataOutList[i][4][1]/pi]
            temp_end = [DataOutList[i][4][0]/pi, DataOutList[i][4][1]/pi]
            temp_start = [DataOutList[i][5][0]/pi, DataOutList[i][5][1]/pi]
            X_many.append(temp_x)
            Y_many.append(temp_y)
            End_point.append(temp_end)
            Start_point.append(temp_start)
            Value_many.append(DataOutList[i][0])

        # Specify the plots

        ######### First plot
        fig, ax = plt.subplots()

        
        # Colours used for the lines are specified by the variable z
        z = ['r-', 'b-', 'g-', 'y-', 'm-']  ## yo  means yellow circle - Not used here

        # Energy breaks used by the colour coding and the detail graphs 


        for i in range(len(X_many)):
            if Value_many[i] < energy_break[0]:
                plot_colour = z[0]
            elif Value_many[i] <energy_break[1]:
                plot_colour = z[1]
            elif Value_many[i] < energy_break[2]:
                plot_colour = z[2]
            elif Value_many[i] < energy_break[3]:
                plot_colour = z[3]
            else:
                plot_colour = z[4]
            plt.plot(X_many[i], Y_many[i], plot_colour)
            plt.plot(End_point[i][0], End_point[i][1], 'b^') # note end point has a triangle to distinguish the end points from the start points
            #plt.plot(Start_point[i][0], Start_point[i][1], 'b*')


        plt.suptitle("Start and end angles with p = 1 for multiple random angle runs ")
        ax.set_title('Key red <{0}, blue < {1}, green < {2}, yellow < {3}, maroon >= {3}'.format(energy_break[0], energy_break[1], energy_break[2], energy_break[3]))
        ax.set_ylabel("Gamma/pi")
        ax.set_xlabel("Beta/pi")
        plt.grid(True)
        fig.set_figwidth(14.0)
        fig.set_figheight(6.0)


        plt.show()


        ######### Second plot - sub set of runs with best EV results

        fig, ax = plt.subplots()
        for i in range(len(X_many)):
            if Value_many[i] <energy_break[0]:
                plot_colour = z[0]
                plt.plot(X_many[i], Y_many[i], plot_colour)
                plt.plot(End_point[i][0], End_point[i][1], 'b^')
            
        ax.set_ylabel("Gamma/pi")
        ax.set_xlabel("Beta/pi")
        plt.suptitle("Start and end angles with p = 1 for best runs")
        ax.set_title('Runs with an EV below {0}'.format(energy_break[0]))
        plt.grid(True)
        fig.set_figwidth(14.0)
        fig.set_figheight(6.0)
        plt.show()

        ######### Third plot - sub set of runs with second best set of EV results

        fig, ax = plt.subplots()
        for i in range(len(X_many)):
            if Value_many[i] > energy_break[0] and  Value_many[i] < energy_break[1]:
                plot_colour = z[1]
                plt.plot(X_many[i], Y_many[i], plot_colour)
                plt.plot(End_point[i][0], End_point[i][1], 'r^')
        plt.suptitle("Start and end angles with p = 1 for runs in second best band")
        ax.set_title('Runs with an EV between {0} and {1}'.format(energy_break[0], energy_break[1])) # Note this format required to include data in the title
        ax.set_ylabel("Gamma/pi")
        ax.set_xlabel("Beta/pi")
        plt.grid(True)
        fig.set_figwidth(14.0)
        fig.set_figheight(6.0)
        plt.show()
