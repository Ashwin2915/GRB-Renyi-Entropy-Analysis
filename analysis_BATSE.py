import pandas as pd 
import numpy as np 
import csv
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import random
import math

def distance(ra1, dec1, dist1, ra2, dec2, dist2):
    t1 = (np.sin(np.deg2rad(dec1)) * np.sin(np.deg2rad(dec2)) * np.cos(np.deg2rad(ra1)) * np.cos(np.deg2rad(ra2)))
    t2 = (np.sin(np.deg2rad(dec1)) * np.sin(np.deg2rad(dec2)) * np.sin(np.deg2rad(ra1)) * np.sin(np.deg2rad(ra2)))
    t3 = (np.cos(np.deg2rad(dec1)) * np.cos(np.deg2rad(dec2)))
    if (dist1**2) + (dist2**2) - ((2*dist1*dist2)*(t1 + t2 + t3))>0:
        dist = ((dist1**2) + (dist2**2) - ((2*dist1*dist2)*(t1 + t2 + t3)))**0.5
    elif ((dist1**2) + (dist2**2) - ((2*dist1*dist2)*(t1 + t2 + t3)))<0:
        if ((ra1==ra2) & (dec1==dec2) & (dist1 == dist2)):
            dist = 0
    else:
        dist = 0
    return dist 

#filename = 'ashwin_data2.csv'
filename = 'BATSE_data.csv'
fields = []
rows = []

with open(filename, 'r') as csvfile:
    # creating a csv reader object
    csvreader = csv.reader(csvfile)
 
    # extracting field names through first row
    fields = next(csvreader)
 
    # extracting each data row one by one
    for row in csvreader:
        rows.append(row)
 
    # get total number of rows
    print("Total no. of rows: %d" % (csvreader.line_num))


# Data that will be used 
redshift = []
Tt90 = []
Tt90_err = []
ra_r = []
dec_r = []

for i in rows:
    redshift.append(float(i[1]))
    Tt90.append(float(i[2]))
    Tt90_err.append(np.e**float(i[8]))

filename = 'ashwin_data.csv'
fields = []
rows = []

with open(filename, 'r') as csvfile:
    # creating a csv reader object
    csvreader = csv.reader(csvfile)
 
    # extracting field names through first row
    fields = next(csvreader)
 
    # extracting each data row one by one
    for row in csvreader:
        rows.append(row)
 
    # get total number of rows
    print("Total no. of rows: %d" % (csvreader.line_num))

for i in rows:
    redshift.append(float(i[1]))
    Tt90.append(float(i[2]))
    Tt90_err.append(np.e**float(i[8]))
    

id2 = []
T90 = []

file1 = open("batse_grb.txt","r")
id1 = []
ra = []
dec = []
std_dev = []
stat = []

raf = []
decf = []
std_devf = []
statf = []
id1f = []  

for i in range(2702):
    a = file1.readline()
    ra.append(float(a.split()[5]))
    dec.append(float(a.split()[6]))
    std_dev.append(float(a.split()[9]))
    stat.append(a.split()[12])
    id1.append(float(a.split()[0]))
file1.close()

file2 = open("duration_table.txt","r")
for i in range(2041):
    a = file2.readline().split()
    id2.append(float(a[0]))
    T90.append(float(a[4]))

for i in range(len(std_dev)):
    if std_dev[i] <= 6:
        raf.append(ra[i])
        decf.append(dec[i])
        std_devf.append(std_dev[i])
        statf.append(stat[i])
        id1f.append(id1[i])


for i in range(len(statf)):
    if statf[i] == 'Y':
        raf.pop(i)
        decf.pop(i)
        std_devf.pop(i)
        id1f.pop(i)

RA = []
DEC = []
t90 = []
err = []

for i in range(len(id1f)):
    for j in range(len(id2)):
        if id1f[i] == id2[j]:
            id1.append(id1f[i])
            t90.append(T90[j])
            RA.append(raf[i])
            DEC.append(decf[i])
            err.append(std_devf[i])


RA_final = []
DEC_final = []
for i in range(len(Tt90)):
    count = 0
    max = 100
    for j in range(len(t90)):
        if ((Tt90[i]+Tt90_err[i]>t90[j]) and (Tt90[i]-Tt90_err[i]<t90[j])):
            if np.abs(Tt90[i]-t90[j])<max:
                max = np.abs(Tt90[i]-t90[j])
                temp_RA = RA[j]
                temp_DEC = DEC[j]
    RA_final.append(temp_RA)
    DEC_final.append(temp_DEC)
    count+=1

    if count == 0:    
        RA_final.append('na')
        DEC_final.append('na')
    if count>1:
        print('problem2')


ra_r = RA_final
dec_r = DEC_final
redhshift = redshift


cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
dist = []
c = 0
for i in redshift:
    dist.append(cosmo.luminosity_distance(i).value)

dist = np.array(dist)

count = 0
for i in dist:
    count+=1

def mock(n, dis):
    ra_mock = []
    dec_mock = []
    dist_mock = []

    for i in range(n):
        u = random.random()
        v = random.random()
        z = random.random()
        d = math.asin(2*u - 1) #changed
        dd = 90 - d*(180/math.pi)
        r = 360*v
        R = dis*(z**(1/3))
        ra_mock.append(r)
        dec_mock.append(dd)
        dist_mock.append(R)
    result = []
    result.append(ra_mock)
    result.append(dec_mock)
    result.append(dist_mock)
    return result 

def analysis(order, radius, distt, ra, dec):

    i = 1000000
    a = mock(len(distt), i)
    print('number of points in database : ', len(ra))
    print('number of points in artificial distribution : ', i)
    ra_mock = a[0]
    dec_mock = a[1]
    dist_mock = a[2]

    result = []
    for rad in radius:
        filtered_dist = []
        filtered_ra = []
        filtered_dec = []
        for i in range(len(distt)):
            if distt[i]<np.max(distt)-rad:
                filtered_dist.append(distt[i])
                filtered_ra.append(ra[i])
                filtered_dec.append(dec[i])

        for i in range(len(filtered_dec)):
            filtered_dec[i] = (filtered_dec[i] - 90) * (-1)
        
        #calculating rho and f
        rho = []
        f = [] 
        for i in range(len(filtered_dist)):
            n = 0
            for j in range(len(filtered_dist)):
                if (distance(filtered_ra[i], filtered_dec[i], filtered_dist[i], filtered_ra[j], filtered_dec[j], filtered_dist[j]) < rad):
                    n+=1
        
            rho.append((3*n)/(4*np.pi*(rad**3)))
        for i in rho:
            f.append(i/np.sum(rho))

        g = []
        for i in range(len(filtered_dist)):
            n = 1
            for j in range(len(dist_mock)):
                if (distance(filtered_ra[i], filtered_dec[i], filtered_dist[i], ra_mock[j], dec_mock[j], dist_mock[j]) < rad):
                    n+=1
            g.append(n/len(dist_mock))

        sum = 0   
        for i in range(len(f)):
            sum = sum + (((f[i])**order)/((g[i])**(order-1)))
        Dq = (1/(order-1))*(np.log(sum))
        result.append(Dq)
        # r = 0
        # l = 0
        # if order == 1:
        #     for i in f:
        #         if i == 0:
        #             r = r+0
        #         else:
        #             r = r+(i*(np.log(i)))
        #     result.append(-r)


        # else:
        #     for i in f:
        #         l = l + i**order    
        #     result.append((1/(1-order))*(np.log10(l)))
    return(result)


l = []
i = 100
while i<50000:
    l.append(i)
    i+=100


i = 1.0000000001


markers = ['>', '+', '.', ',', 'o', 'v', 'x', 'X', 'D', '|']
orders = [float(i),2,3,4,5,6,7,8,9,10]
matrix = []
for i in orders:
    result = (analysis(i,l, dist, ra_r, dec_r))
    result = np.array(result)
    result_f = []
    y_f = []
    j = 9
    while j<len(l):
        result_f.append(result[j])
        y_f.append(int(100+j*100))
        j+=10
    matrix.append(result)
    plt.plot(y_f,result_f, label = 'order = ' + str(int(i)), linestyle = '--', marker = markers[int(i)-1])
for i in range(len(matrix[0])):
    if np.abs(matrix[0][i] - matrix[9][i]) < 0.05:
        plt.axvline(x=(i*100 + 100), color='red', linestyle=':', label='convergence line')
        break
plt.xlabel('r (MPc)')
#plt.ylabel('Sq(r)/Sq(r)max')
plt.ylabel('Dq(r)') 
plt.legend()
plt.title('BATSE GRB catalog')
plt.grid()
plt.savefig('BATSE_Dq(r)_.png', dpi=300)
plt.show()