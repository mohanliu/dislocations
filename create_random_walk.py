#! /usr/bin/python
# Created by Mohan Liu, Oct 3rd, 2017

import numpy as np
import random
import sys

# Global variable: C-C single bond length and bond angle
a = 1.54 # Angstrom
angle = 109.5 # degree

def perpen_direction(vec):
    """
    Input:
        vec: A 3D numpy array 
    Return:
        A random unit vector that is perpendicular to the input vector
    """
    if vec[0] == vec[1] == vec[2] == 0:
        raise ValueError('zero-vector')

    # p1 is one vector perpendicular to input vector
    if vec[0] == 0:
        p1 = np.array([1,0,0])
    elif vec[1] == 0:
        p1 = np.array([0,1,0])
    elif vec[2] == 0:
        p1 = np.array([0,0,1])
    else:
        p1 = np.array([1,1,-1.0*(vec[0]+vec[1])/vec[2]])

    # p2 is another vector perpendicular to both p1 and input vector
    p2 = np.cross(vec,p1)

    # normalize p1 and p2
    p1_norm = p1/np.linalg.norm(p1)
    p2_norm = p2/np.linalg.norm(p2)

    # randomize a theta value 
    random_theta = random.uniform(-np.pi,np.pi)

    return p1_norm*np.cos(random_theta)+p2_norm*np.sin(random_theta)

def check_overlap(chain,vec,tol_fac=1.3):
    """
    Input:
        chain [list of numpy arrays]: an existed carbon chain         
        vec [numpy array]: postion of additional carbon 
        tol_fac [float]: tolerance factor (the shortest distance between                                   any two unbonded carbons should be a * tol_fac)
    Return:
        True: no overlaps in the carbon chain 
        False: otherwise
    """
    for v in chain[:-1]:
        if not np.linalg.norm(v-vec) > a*tol_fac:
            return False
    return True

def add_another_carbon(chain):
    """
    Input:
        chain [list of numpy arrays]: an existed carbon chain         
    Function:
        Append a new carbon atom to the existed chain
    """
    last1 = chain[-1]
    last2 = chain[-2]
    v_diff = last1-last2

    # Get new carbon position (make sure new atom doesn't overlap)
    phi = (180.0-angle)/180*np.pi
    v1 = v_diff*np.cos(phi)
    v2 = perpen_direction(v1)*a*np.sin(phi)
    new = last1+v1+v2
    while not check_overlap(chain,new):
        v2 = perpen_direction(v1)*a*np.sin(phi)
        new = last1+v1+v2

    chain.append(new)

def write_xyz(chain,filename='result.xyz'):
    """
    Input:
        chain [list of numpy arrays]: an existed carbon chain         
        filename [string]
    Function:
        Write the atomic postions of the carbon chain into a .xyz file
    """
    f = open(filename,'w')
    f.write('%d\n' %len(chain))
    f.write('This is a carbon chain\n')
    
    for v in chain:
        f.write('C %.3f   %.3f   %.3f\n'%(v[0],v[1],v[2]))
    f.close()

def write_txt(chain,filename='result.txt'):
    """
    Input:
        chain [list of numpy arrays]: an existed carbon chain         
        filename [string]
    Function:
        Write the atomic postions of the carbon chain into a .txt file
    """
    f = open(filename,'w')
    for v in chain:
        f.write('%.3f   %.3f   %.3f\n'%(v[0],v[1],v[2]))
    f.close()

if __name__ == "__main__":
    # initialize the carbon chain (start with two carbons)
    carbon_chain = [np.array([0,0,0]),np.array([0,0,a])]

    # number of carbons atoms in the chain (default is 50)
    try:
        N = int(sys.argv[1])
    except:
        N = 50

    print "Creating a chain with %d carbons" %N

    for i in range(N-2):
        add_another_carbon(carbon_chain)

    write_xyz(carbon_chain,'C'+str(N)+'.xyz')
    write_txt(carbon_chain,'C'+str(N)+'.txt')
    print "Done!"
