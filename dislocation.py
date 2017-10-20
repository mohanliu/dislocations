#! /usr/bin/python
# Created by Mohan, Oct 4th, 2017

import numpy as np

# unit cell length
a = 1.5

# Function to create the original crystal structure
def create_mesh(N=7):
    points = np.mgrid[-N:N,-N:N,-N:N]
    p = points.reshape(3,-1).T
    return p

# Function to apply screw dislocation
def screw_dislocation(data,b=-a):
    new_data = []
    for d in data:
        if d[1] > -0/5*a:
            z = d[2]-b/np.pi*np.arctan((d[1]+0.5*a)/(d[0]+0.5*a))
            new_data.append(np.array([d[0],d[1],z]))
        else:
            new_data.append(d)
    return np.array(new_data)

# Function to write .xyz file
def write_xyz(data,filename='result.xyz'):
    f = open(filename,'w')
    f.write('%d\n' %len(data))
    f.write('Structure with screw dislocation\n')
    
    for v in data:
        f.write('Po %.3f   %.3f   %.3f\n'%(v[0],v[1],v[2]))
    f.close()

# Function to write .txt file
def write_txt(data,filename='result.txt'):
    f = open(filename,'w')
    for v in data:
        f.write('%.3f   %.3f   %.3f\n'%(v[0],v[1],v[2]))
    f.close()



if __name__ == "__main__":
    data = create_mesh()*a
    write_xyz(data,'before_dislocation.xyz')
    new_data = screw_dislocation(data)
    write_xyz(new_data,'after_dislocation.xyz')


