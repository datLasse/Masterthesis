import csv
import numpy as np


datatime = np.zeros(59)
datapoints = np.zeros(59)

with open('data.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        datatime[i] = row[0]
        datapoints[i] = row[1]
        i += 1

print(datatime,datapoints)
