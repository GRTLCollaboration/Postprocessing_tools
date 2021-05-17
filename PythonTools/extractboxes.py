import csv
import matplotlib.pyplot as plt

'''Run processing.sh with pout.0 - outputs level, time, boxes on proc 0, total boxes'''
'''Only works on output with verbosity 0'''

i=0

data = open('pout.0_parsed')
with data as csvfile:
    filereader = csv.reader(csvfile, delimiter=',')
    timelist = [] #col 1
    totalboxlist = [] #col 3

    '''Hardcoded max level to plot'''
    while i < 10:
        for row in filereader:
            if int(row[0]) == i:
                '''Can input max time'''
#                if float(row[1]) < 300:
                    timelist.append(float(row[1]))
                    totalboxlist.append(int(row[3]))

        plt.plot(timelist, totalboxlist, label='%s' % i)
        plt.legend()
        plt.savefig('boxplot_level_%s.png' % i)

        timelist.clear()
        totalboxlist.clear()

        '''Resets file loop'''
        data.seek(0)
        
        i+=1

    print('%s' % i)


