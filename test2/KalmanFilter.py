"""

Extended kalman filter (EKF) localization sample

author: Atsushi Sakai (@Atsushi_twi)

"""

import numpy as np
import math
import matplotlib.pyplot as plt

class Observe (object):
    
    def __init__(self):
        self.Q    = np.diag([0.1, 0.1, np.deg2rad(1.0), 1])**2
        self.R    = np.diag([1.0, np.deg2rad(40.0)])**2
        self.Qsim = np.diag([0.5, 0.5])**2
        self.Rsim = np.diag([1.0, math.radians(30.0)])**2
        self.dt = 0.1
        self.mat_f = np.identity(4)
        self.mat_f[3,3] = 0.0
        self.Jacob = np.matrix ([
            [1, 0, 0, 0],
            [0, 1, 0, 0]
        ])
    
    def Model (self, x=np.ones(4), u=np.ones(2)):
        mat_b = np.zeros ((4,2))
        mat_b[0, 1] = self.dt * np.cos(x[2,0])
        mat_b[1, 1] = self.dt * np.sin(x[2,0])
        mat_b[2, 1] = self.dt
        mat_b[3, 0] = 0.0
        return self.mat_f * x + mat_b * u
    
    def JacobF (self, x=np.ones(4), u=np.ones(2)):
        """
        Jacobian of Motion Model
        motion model
        x_{t+1} = x_t+v*dt*cos(yaw)
        y_{t+1} = y_t+v*dt*sin(yaw)
        yaw_{t+1} = yaw_t+omega*dt
        v_{t+1} = v{t}
        so
        dx/dyaw = -v*dt*sin(yaw)
        dx/dv = dt*cos(yaw)
        dy/dyaw = v*dt*cos(yaw)
        dy/dv = dt*sin(yaw)
        """
        jF = np.identity (4)
        jF[0, 2] =-self.dt * np.sin(x[2,0])
        jF[0, 3] = self.dt * np.cos(x[2,0])
        jF[1, 2] = self.dt * np.cos(x[2,0])
        jF[1, 3] = self.dt * np.sin(x[2,0])
        return jF

class EKF_Calc (Observe):

    def __init__(self):
        super(EKF_Calc, self).__init__()
        self.Init()
    
    def Init (self):
        # State Vector [x y yaw v]'
        self.xEst  = np.matrix(np.zeros((4, 1)))
        self.xTrue = np.matrix(np.zeros((4, 1)))
        self.PEst  = np.eye(4)
        self.xDR   = np.matrix(np.zeros((4, 1)))  # Dead reckoning

        # history
        self.hxEst  = self.xEst
        self.hxTrue = self.xTrue
        self.hxDR   = self.xTrue
        self.hz     = np.zeros((1, 2))
    
    def Observe_model (self, u):
        self.xTrue = self.Model (self.xTrue, u)
        
        zx = self.xTrue[0, 0] + np.random.randn() * self.Qsim[0, 0]
        zy = self.xTrue[1, 0] + np.random.randn() * self.Qsim[1, 1]
        z = np.matrix([zx, zy])

        ud1 = u[0, 0] + np.random.randn() * self.Rsim[0, 0]
        ud2 = u[1, 0] + np.random.randn() * self.Rsim[1, 1]
        ud = np.matrix([ud1, ud2]).T
        
        self.xDR = self.Model (self.xDR, ud)
        return z, ud

    def Estmiate (self, z, u):
        #  Predict
        xPred = self.Model (self.xEst, u)
        jF = self.JacobF (xPred, u)
        PPred = jF * self.PEst * jF.T + self.Q

        #  Update
        jH = self.Jacob
        zPred = self.Jacob * xPred
        y = z.T - zPred
        S = jH * PPred * jH.T + self.R
        K = PPred * jH.T * np.linalg.inv(S)
        self.xEst = xPred + K * y
        self.PEst = (np.eye(len(self.xEst)) - K * jH) * PPred

        self.hxEst  = np.hstack((self.hxEst , self.xEst))
        self.hxTrue = np.hstack((self.hxTrue, self.xTrue))
        self.hxDR   = np.hstack((self.hxDR  , self.xDR))
        self.hz     = np.vstack((self.hz    , z))

def plot_covariance_ellipse(xEst, PEst):
    Pxy = PEst[0:2, 0:2]
    eigval, eigvec = np.linalg.eig(Pxy)

    if eigval[0] >= eigval[1]:
        bigind = 0
        smallind = 1
    else:
        bigind = 1
        smallind = 0

    t = np.arange(0, 2 * math.pi + 0.1, 0.1)
    a = math.sqrt(eigval[bigind])
    b = math.sqrt(eigval[smallind])
    x = [a * math.cos(it) for it in t]
    y = [b * math.sin(it) for it in t]
    angle = math.atan2(eigvec[bigind, 1], eigvec[bigind, 0])
    R = np.matrix([[math.cos(angle), math.sin(angle)],
                   [-math.sin(angle), math.cos(angle)]])
    fx = R * np.matrix([x, y])
    px = np.array(fx[0, :] + xEst[0, 0]).flatten()
    py = np.array(fx[1, :] + xEst[1, 0]).flatten()
    plt.plot(px, py, "--r")

if __name__ == '__main__':
    time = 0.0
    obj = EKF_Calc ()

    while 50.0 >= time:
        time += obj.dt
        u = np.matrix([1.0, 0.1]).T

        z, ud = obj.Observe_model(u)
        obj.Estmiate(z, ud)

        plt.cla()
        plt.plot(obj.hz[:, 0], obj.hz[:, 1], ".g")
        plt.plot(np.array(obj.hxTrue[0, :]).flatten(),
                 np.array(obj.hxTrue[1, :]).flatten(), "-b")
        plt.plot(np.array(obj.hxDR  [0, :]).flatten(),
                 np.array(obj.hxDR  [1, :]).flatten(), "-k")
        plt.plot(np.array(obj.hxEst [0, :]).flatten(),
                 np.array(obj.hxEst [1, :]).flatten(), "-r")
        plot_covariance_ellipse(obj.xEst, obj.PEst)
        plt.axis("equal")
        plt.grid(True)
        plt.pause(0.001)
    
    plt.savefig ("name.png")