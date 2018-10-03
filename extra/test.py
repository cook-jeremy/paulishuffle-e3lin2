import matplotlib.pyplot as plt
import sys
import csv

x = []
y = []

with open('test','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))

    plt.plot(x,y)
    plt.xlabel('gamma')
    plt.ylabel('<C>')
    plt.show()
