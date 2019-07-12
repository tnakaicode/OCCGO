import numpy as np
from abc import abstractmethod
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import random
import math
import datetime

class Particle():
    def __init__(self, pos = [0.0, 0.0, 0.0], vel = [0.0, 0.0, 0.0], force = [0.0, 0.0, 0.0], mass = 1.0, type = -1):
        self.pos = np.array(pos)
        self.vel = np.array(vel)
        self.force = np.array(force)
        self.mass = mass
        self.type = type

    # レナードジョーンズポテンシャル下で粒子間に働く力
    @staticmethod
    def force(pos1, pos2, eps = 1.0, sigma = 1.0):
        #return np.array([0, 0, 0])
        r2 = ((pos2 - pos1)**2).sum()
        if abs(r2) < 0.00001:  
            return np.array([0.0, 0.0, 0.0])
        sig6_div_r6 = (sigma**2 / r2)**3
        #print("force : ", 24 * eps * (1/r2) * sig6_div_r6 * (1 - 2 * sig6_div_r6) * (pos2 - pos1))
        return 24 * eps * (1/r2) * sig6_div_r6 * (1 - 2 * sig6_div_r6) * (pos2 - pos1)

class Calculater():
    # 粒子のdt後の位置(pos)と速度(vel)をルンゲクッタ法で更新する
    @staticmethod
    def runge_kutta(particle, adjuscent_particles, ps, dt, t, eps, sigma):
        k1 = dt * particle.vel
        l1 = 0
        for p in adjuscent_particles:
            l1 += dt * Particle.force(pos1 = particle.pos, pos2 = p.pos, eps = eps, sigma = sigma) / particle.mass
        l1 += dt * ps.force(particle.pos, particle.vel, particle, t)

        k2 = dt * (particle.vel + k1 / 2)
        l2 = 0
        for p in adjuscent_particles:
            l2 = dt * Particle.force(particle.pos + l1/2, p.pos, eps, sigma) / particle.mass
        l2 += dt * ps.force(particle.pos + l1/2, particle.vel + k1/2, particle, t + dt/2)

        k3 = dt * (particle.vel + k2 / 2)
        l3 = 0
        for p in adjuscent_particles:
            l3 = dt * Particle.force(particle.pos + l2/2, p.pos, eps, sigma) / particle.mass
        l3 += dt * ps.force(particle.pos + l2/2, particle.vel + k2/2, particle, t + dt/2)

        k4 = dt * (particle.vel + k3)
        l4 = 0
        for p in adjuscent_particles:
            l4 = dt * Particle.force(particle.pos + l3, p.pos, eps, sigma) / particle.mass
        l4 += dt * ps.force(particle.pos + l3, particle.vel + k3, particle, t + dt)

        particle.pos += (k1 + 2*k2 + 2*k3 + k4)/6
        particle.vel += (l1 + 2*l2 + 2*l3 + l4)/6

    # 粒子のdt後の位置(pos)と速度(vel)をオイラー法で更新する
    @staticmethod
    def euler(particle, adjuscent_particles, ps, dt, t, eps, sigma):
        f = 0
        for p in adjuscent_particles:
            f += dt * Particle.force(pos1 = particle.pos, pos2 = p.pos, eps = eps, sigma = sigma) / particle.mass
        f += dt * ps.force(particle.pos, particle.vel, particle, t)
        particle.pos += dt * particle.vel
        particle.vel += dt * f

class ParticleSystem():
    # 
    def __init__(self, cutoff_r, NX, NY, NZ, e = 0.9, mass = 1, eps = 1, sigma = 1):
        self.cutoff_r = cutoff_r
        self.NX = NX
        self.NY = NY
        self.NZ = NZ
        self.X_MAX = cutoff_r * NX
        self.Y_MAX = cutoff_r * NY
        self.Z_MAX = cutoff_r * NZ
        self.e = e
        self.mass = mass
        self.eps = eps
        self.sigma = sigma
        self.time = 0.0
        self.prev_timestamp = 0.0
        self.fps = 0.1
        self.particles = []
        self.snapshots = []
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.gif_output_mode = False
        self.init_particles()
        self.update_cells()
        self.setup_particle_type_dict()

    # 全粒子に働く力とdt後の位置・速度を求める。
    def update(self, dt, calc_method):
        self.time += dt
        #update_forces_on_particle()
        self.update_pos_and_vel(dt, calc_method)
        self.update_cells()

    '''
    # 全粒子に働く力を計算する。
    # ルンゲクッタ法では使えそうにないため廃止。
    def update_forces_on_particle():
        for particle in self.particles:
            particle.force = 0
            # particle上に働く粒子間力を計算
            for p in get_particles_adjuscent_cell(particle):
                f = Particle.force(particle.pos, p.pos)
                particle.force += f
            particle.force += self.force(particle)
    '''

    # 全粒子のdt後の位置と時間をcalc_methodで計算する。
    def update_pos_and_vel(self, dt, calc_method):
        for particle in self.particles:
            calc_method(particle = particle, adjuscent_particles = self.get_particles_adjuscent_cell(particle), ps = self, dt = dt, t = self.time, eps = self.eps, sigma = self.sigma)
            # スペースの両端で跳ね返るようにする。(跳ね返り係数e)
            if particle.pos[0] < 0:
                #particle.pos[0] += 2*(0 - particle.pos[0])
                particle.pos[0] = 0.0
                particle.vel[0] = -self.e * particle.vel[0]
            if particle.pos[0] > self.X_MAX:
                #particle.pos[0] -= 2*(particle.pos[0] - self.X_MAX)
                particle.pos[0] = self.X_MAX - 0.0001
                particle.vel[0] = -self.e * particle.vel[0]
            if particle.pos[1] < 0:
                #particle.pos[1] += 2*(0 - particle.pos[1])
                particle.pos[1] = 0.0
                particle.vel[1] = -self.e * particle.vel[1]
            if particle.pos[1] > self.Y_MAX:
                #particle.pos[1] -= 2*(particle.pos[1] - self.Y_MAX)
                particle.pos[1] = self.Y_MAX - 0.0001
                particle.vel[1] = -self.e * particle.vel[1]
            if particle.pos[2] < 0:
                #particle.pos[2] += 2*(0 - particle.pos[2])
                particle.pos[2] = 0.0
                particle.vel[2] = -self.e * particle.vel[2]
            if particle.pos[2] > self.Z_MAX:
                particle.pos[2] = self.Z_MAX - 0.0001
                #particle.pos[2] -= 2*(particle.pos[2] - self.Z_MAX)
                particle.vel[2] = -self.e * particle.vel[2]
        if self.gif_output_mode is True:
            self.take_snapshot(save_to_snapshots = True)

    # particleセル含むセルと隣接26セル内に含まれる粒子リストを取得
    def get_particles_adjuscent_cell(self, particle):
        particles = []
        index = self.get_cellindex_from_pos(particle.pos)
        for i in range(index[0] - 1, index[0] + 1):
            for j in range(index[1] - 1, index[1] + 1):
                for k in range(index[2] - 1, index[2] + 1):
                    try:
                        for p in self.cell[i][j][k] if self.cell[i][j][k] is not None else []:
                            if p is not particle: 
                                particles.append(p)
                    except IndexError:
                        continue
        return particles

    # 粒子の位置からセルのインデックスを算出する。
    def get_cellindex_from_pos(self, pos):
        return [int(pos[0] / self.cutoff_r), int(pos[1] / self.cutoff_r), int(pos[2] / self.cutoff_r)]

    # 各粒子がどのセルに属しているか更新
    def update_cells(self):
        self.cell = [[[ [] for i in range(self.NX)] for j in range(self.NY)] for k in range(self.NZ)]
        for p in self.particles:
            try:
                index = self.get_cellindex_from_pos(p.pos)
                self.cell[index[0]][index[1]][index[2]].append(p)
            except Exception as e:
                print("Exception : ", e)
                print("pos : ", p.pos)
                #input()

    # 粒子リストを取得
    def get_particles(self):
        return self.particles

    # 時間を取得
    def get_time(self):
        return self.time

    # 粒子の状態のスナップショットをとる。表示は行わない。
    def take_snapshot(self, save_to_snapshots = False):
        if self.time - self.prev_timestamp > self.fps:
            self.ax.set_title("Time : {}".format(self.time))
            self.ax.set_xlim(0, self.X_MAX)
            self.ax.set_ylim(0, self.Y_MAX)
            self.ax.set_zlim(0, self.Z_MAX)
            scats = []
            for type, particles in self.particle_dict.items():
                x_list = []
                y_list = []
                z_list = []
                for p in particles[0]:
                    x_list.append(p.pos[0])
                    y_list.append(p.pos[1])
                    z_list.append(p.pos[2])
                scats.append(self.ax.scatter(x_list, y_list, z_list, c=particles[1]))
            if save_to_snapshots is True:
                self.snapshots.append(scats)
                print(len(self.snapshots))
            self.prev_timestamp = self.time

    # 粒子の状態を表示
    def show_snapshot(self):
        if self.time - self.prev_timestamp > 0.1:
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111, projection='3d')
            self.take_snapshot()
            #plt.savefig('box_gravity_time-{}.png'.format(self.time))
            plt.show()
            self.prev_timestamp = self.time

    # gifファイルへ書き込みモードオン
    def start_output_gif(self, fps = 0.1):
        self.gif_output_mode = True
        self.snapshots = []
        self.take_snapshot(save_to_snapshots = True)

    # gifファイル書き込みモードオフ。ファイルへ書き込み
    def stop_output_gif(self, filename = "hoge.gif"):
        print("stop_output_gif : ", len(self.snapshots))
        self.gif_output_mode = False
        ani = animation.ArtistAnimation(self.fig, self.snapshots)
        ani.save(filename, writer='imagemagick')
        self.snapshots = []

    # 粒子のタイプ毎のディクショナリを作っておく。
    def setup_particle_type_dict(self):
        self.particle_dict = {}
        for p in self.particles:
            if str(p.type) not in self.particle_dict:
                # 粒子のタイプごとにランダムな色を割り当てる
                self.particle_dict[str(p.type)] = [[], tuple([random.random() for _ in range(3)])]
            self.particle_dict[str(p.type)][0].append(p)

    # 全粒子の運動エネルギーを計算
    def get_kinetic_energy(self):
        return  sum([(p.mass*np.array(p.vel)**2).sum()/2 for p in self.particles])

    # 粒子の初期化
    @abstractmethod    
    def init_particles(self):
        raise NotImplementedError()

    # 系として働く力 (重力など)
    @abstractmethod
    def force(self, pos, vel, particle, t):
        raise NotImplementedError()

# 重力下にあるボックス状の粒子軍の系
class BoxGravitySystem(ParticleSystem):
    def __init__(self, cutoff_r, NX, NY, NZ, e, mass, eps, sigma):
        super().__init__(cutoff_r, NX = NX, NY = NY, NZ = NZ, e = e, mass = mass, eps = eps, sigma = sigma)

    # 下端にボックス状の10×10×10の粒子軍を生成。初速度は斜め上に進むように設定。
    def init_particles(self):
        for x in np.linspace(0.1*self.X_MAX, 0.2*self.X_MAX, 10):
            for y in np.linspace(0.1*self.Y_MAX, 0.2*self.Y_MAX, 10):
                for z in np.linspace(0.1*self.Z_MAX, 0.2*self.Z_MAX, 10):
                    self.particles.append(Particle(pos = [x, y, z], vel=[0.1*self.X_MAX, 0.05*self.Y_MAX, 0.5*self.Z_MAX], mass = self.mass))

    # 重力
    def force(self, pos, vel, particle, t):
        return np.array([0.0, 0.0, -particle.mass * 9.8])

# 壁に弾丸を打ち込む系
class BulletWallSystem(ParticleSystem):
    def __init__(self, cutoff_r, NX, NY, NZ, e, mass, eps, sigma):
        super().__init__(cutoff_r, NX = NX, NY = NY, NZ = NZ, e = e, mass = mass, eps = eps, sigma = sigma)

    # 初速度を持った球(弾丸)と固定された壁を設置
    def init_particles(self):
        # 球(弾丸)生成
        bullet_center = [0.35*self.X_MAX, 0.5*self.Y_MAX, 0.5*self.Z_MAX]
        for i in range(200):
            r = 0.05 * self.X_MAX * (random.random() - 0.5) * 2
            phi = 2 * np.pi * random.random()
            theta = np.pi * random.random()
            self.particles.append(Particle(pos = [bullet_center[0] + r*np.sin(theta)*np.cos(phi), bullet_center[1] + r*np.sin(theta)*np.sin(phi), bullet_center[2] + r*np.cos(theta)], vel=[0.2*self.X_MAX, 0, 0], mass = self.mass, type = 1))
        # 壁生成
        for x in np.linspace(0.49*self.X_MAX, 0.50*self.X_MAX, 2):
            for y in np.linspace(0.3*self.Y_MAX, 0.7*self.Y_MAX, 30):
                for z in np.linspace(0.3*self.Z_MAX, 0.7*self.Z_MAX, 30):
                    self.particles.append(Particle(pos = [x, y, z], vel=[0.0, 0.0, 0.0], mass = self.mass, type = 2))

    # 外部から受ける力はなし
    def force(self, pos, vel, particle, t):
        return np.array([0.0, 0.0, 0])

# 竜巻のように回転するような力を受ける系
class TornadeSystem(ParticleSystem):
    def __init__(self, cutoff_r, NX, NY, NZ, e, mass, eps, sigma):
        super().__init__(cutoff_r, NX = NX, NY = NY, NZ = NZ, e = e, mass = mass, eps = eps, sigma = sigma)

    # ランダムな位置に3000個粒子を作成
    def init_particles(self):
        for _ in range(3000):
            x = self.X_MAX * random.random()
            y = self.Y_MAX * random.random()
            z = self.Z_MAX * random.random()
            self.particles.append(Particle(pos = [x, y, z], vel=[0.0, 0.0, 0.0], mass = self.mass, type = 1))

    # z軸を中心に回転するような力(rot(r)∝[0,0,1]となるような力) × z
    def force(self, pos, vel, particle, t):
        return 0.05 * pos[2] * np.array([-(pos[1] - self.Y_MAX/2), pos[0] - self.X_MAX/2 , 0.0])# - 0.5 * np.array(vel)

if __name__ == '__main__':
    cutoff_r = 1.0
    e = 0.1
    mass = 1.0
    eps = 1
    sigma = 0.5
    dt = 0.01
    system = BulletWallSystem(cutoff_r = cutoff_r, NX = 100, NY = 100, NZ = 100, e = e, mass = mass, eps = eps, sigma = sigma)
    time = 0
    system.start_output_gif(fps = 0.1)
    while time <= 2.5:
        system.update(dt = dt, calc_method = Calculater.runge_kutta)
        time = system.get_time()
        print("Time : ", time)
        #system.show_snapshot()
    system.stop_output_gif(filename = "BoxGravitySystem_cutoffr-{}_mass-{}_eps-{}_sigma-{}_dt-{}.gif".format(cutoff_r, mass, eps, sigma, dt))
