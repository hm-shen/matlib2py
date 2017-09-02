import numpy as np

from filtering import *
from utils import *
import matplotlib.pyplot as plt



if __name__ == '__main__':

    ''' Test filtering package '''

    dT = 30
    N = 7 * 24 * 60 / dT
    t = np.linspace(0,N-1,N) * dT
    pnoise = 0.3
    T1 = 12.4 * 60
    T2 = 24 * 60
    T3 = 15 * 24 * 60
    Tc = 10 * 60
    xn = 5 + 3 * np.cos(2 * np.pi * t / T1) + 2 * np.cos(2 * np.pi * t / T2) +\
            np.cos(2 * np.pi * t / T3)
    xn += pnoise * np.max(xn - np.mean(xn)) * (0.5 - np.random.rand(xn.size))
    [xs, c, h, Cx, f] = lanczos_filter(xn, spl_intvl=dT, cut_off_freq=1.0/Tc,\
            num_of_coeffs=100,mode='low')

    plt.figure(1)
    plt.plot(t,xn)
    plt.plot(t,xs)

    plt.figure(2)
    plt.plot(f,h)
    plt.plot(f,np.absolute(Cx) / np.max(np.absolute(Cx)))

    plt.show()



