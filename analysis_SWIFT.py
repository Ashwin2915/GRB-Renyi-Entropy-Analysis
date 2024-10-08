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

def ra(ra):
    l = ra.split()
    for i in range(len(l)):
        l[i] = float(l[i])
    degree = 0
    degree = degree + (l[0]*360/24) + (l[1]*360/(24*60))+ (l[2]*360/(24*60*60))
    return degree 

def dec(dec):
    l = dec.split()
    degree = 0
    if float(l[0])<0:
        degree = degree + float(l[0]) - (float(l[1])/60) - (float(l[2])/(60*60)) 
    else:
        degree = degree + float(l[0]) + (float(l[1])/60) + (float(l[2])/(60*60)) 
    return degree




filename = 'SWIFT_Redshift.csv'
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
ra_r = []
dec_r = []

for i in rows:
    if i[9][len(i[9])-1] == ' ':
        redshift.append(None)
    else:
        redshift.append(float(i[9]))

count = 0
indices = []
for i in range(len(redshift)):
    if redshift[i] != None:
        count+=1
        indices.append(i)
#print(len(fields))

for i in rows:
    ra_r.append(ra(i[3]))
    dec_r.append(dec(i[4]))
    if dec(i[4]) < -180:
        print(i)



cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
dist = []
c = 0
for i in indices:
    dist.append(cosmo.luminosity_distance(redshift[i]).value)

#dist = np.array(dist)
# print(count)
# print(len(ra_r))
# print(len(indices))
# print(len(dec_r))
# print(min(dec_r))
# print(max(dec_r))
# print(len(dist))
print(max(dist))
print(min(dist))
print(np.average(dist))

count = 0
for i in dist:
    count+=1
print(count)

l = []
i = 1000
while i<20000:
    l.append(i)
    i+=1000


# filtered_dist = []
# filtered_ra = []
# filtered_dec = []

# for i in range(len(dist)):
#     if dist[i]<18700-radius:
#         filtered_dist.append(dist[i])
#         filtered_ra.append(ra_r[i])
#         filtered_dec.append(dec_r[i])

# for i in range(len(filtered_dec)):
#     filtered_dec[i] = (filtered_dec[i] - 90) * (-1)

# rho = []
# f = []

# for i in range(len(filtered_dist)):
#     n = 0
#     for j in range(len(filtered_dist)):
#         if i!=j:
#             if (distance(filtered_ra[i], filtered_dec[i], filtered_dist[i], filtered_ra[j], filtered_dec[j], filtered_dist[j]) < radius):
#                 n+=1
#     rho.append((3*n)/(4*np.pi*(radius**3)))

# for i in rho:
#     f.append(i/np.sum(rho))

#q = 1

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
            if distt[i]<max(distt)-rad:
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
                if i!=j:
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

def find_max(l):
    max = 0
    for i in range(len(l)):
        if ((l[i]  > max) and (l[i] != np.inf)):
            max = l[i] 
    return max

dr_result = []

# markers = ['>', '+', '.', ',', 'o', 'v', 'x', 'X', 'D', '|']
# orders = [1,2,3,4,5,6,7,8,9,10]
# for i in orders:
#     result = (renyi(i,l, dist, ra_r, dec_r))
#     result = np.array(result)
#     d_result = [1]*len(result)
#     d_result = np.array(d_result)
#     result = result/find_max(result)
#     dr_result = d_result - result
#     #print(result)
#     plt.plot(l,result, label = 'order = ' + str(i),linestyle = '--', marker = markers[i-2])
# plt.xlabel('r (MPc)')
# #plt.ylabel('Sq(r)/Sq(r)max')
# plt.ylabel('Rq(r)') 
# plt.legend()
# plt.title('Swift GRB catalog')
# plt.grid()
# plt.show()


############### HOMOGENEOUS PART 
# homo_dist = []
# homo_ra = []
# homo_dec = []
# #R = max(dist)
# R = 20000
# N = 0
# while N<500:
#     r = np.random.uniform(low = 0.0, high = R) 
#     p = 3 * (r**2)/(R**3)
#     p_rand = np.random.uniform(low = 0.0, high = 3/R)
#     if p_rand >= p:
#         homo_dist.append(r)
#         u = random.random()
#         v = random.random()
#         d = math.acos(1-2*u)
#         dd = 90 - d*(180/math.pi)
#         r = 360*v
#         homo_ra.append(r)
#         homo_dec.append(dd)
#         N+=1
# print(len(homo_dist))
# print(len(homo_ra))
# print(len(homo_dec))

markers = ['>', '+', '.', ',', 'o', 'v', 'x', 'X', 'D', '|']

l = []
i = 200
while i<50000:
    l.append(i)
    i+=200

i = 1.00000001

orders = [float(i),2,3,4,5,6,7,8,9,10]
result_rr = [0]*19
result_rr = np.array(result_rr)
matrix = []
for i in orders:
    # result_r = [[0]*19]*10
    # result_r = np.array(result_r)
    result = (analysis(i,l, dist, ra_r, dec_r))
    result = np.array(result)
    result_f = []
    y_f = []
    j = 9
    while j<len(result):
        result_f.append(result[j])
        y_f.append(200+(j*200))
        j+=10
    matrix.append(result)
    plt.plot(y_f,result_f, label = 'order = ' + str(int(i)), linestyle = '--', marker = markers[int(i)-1])
    print('still going onnnnnn')
i = 10
while i< (len(matrix[0])):
    if np.abs(matrix[0][i] - matrix[9][i]) < 0.05:
        plt.axvline(x=(i*200 + 200), color='red', linestyle=':', label='convergence line')
        print('i is ')
        print(i)
        break
    i = i+1
l = []
i = 0
while i<50000:
    l.append(i)
    i+=5000

plt.xticks(l)
plt.xlabel('r (MPc)')
plt.ylabel('Dq(r)')
plt.legend()
plt.title('SWIFT GRB distribution')
plt.grid()
plt.savefig('SWIFT_Dq(r).png', dpi=300)
plt.show()
