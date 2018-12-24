#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 09:04:30 2017

@author: Yaojia Huang
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from random import randint
import scipy.odr 


class particle(object):

      def __init__(self, mass, R, pos, vel, clr = 'r'):
          self.__m = mass
          self.__R = R
          self.__pos = np.array(pos, dtype = 'float')
          self.__vel = np.array(vel, dtype = 'float')
          if type(self) is container:
             self.__patch = plt.Circle((0,0), -self.__R, 
                            ec = 'b', fill = False, ls = 'solid')
          else:
             self.__patch = plt.Circle(self.__pos, self.__R, fc = clr, ec = 'black')


      def __repr__(self):
          return "Particle (Mass: %r, Radius: %r)" %(self.__m, self.__R)


      def pos(self):
          return self.__pos


      def vel(self):
          return self.__vel


      def KE(self):
          return 0.5 * self.__m * np.dot(self.__vel, self.__vel)


      def momentum(self):
          return self.__m * self.__vel
    

      def move(self, dt):
          '''Let the ball move for dt'''
          self.__pos = self.__pos + self.__vel * dt
          
          '''Update the patch of the ball with the new position'''
          self.__patch.center = self.__pos
          

      def get_patch(self):
          return self.__patch
      
        
      def momentum_change(self, new_vel):
          vel_change = new_vel - self.vel()
          momentum_c = self.__m * np.sqrt(np.dot(vel_change, vel_change))
          return momentum_c


      def time_to_collision(self, other):
          '''
          Solve |pos1+vel1*t-pos2-vel2*t| = R1 + R2
          Use relative velocity and position to calculate the
          coefficients of the quadratic equation: at^2+bt+c=0
          '''
          vel_r = self.__vel - other.__vel
          pos_r = self.__pos - other.__pos
          a = np.dot(vel_r, vel_r)
          if a == 0:
              return np.Inf
          b = 2 * np.dot(pos_r, vel_r)
          c = np.dot(pos_r, pos_r) - (other.__R + self.__R)**2
                        
          determinant = b**2-4*a*c
          
          '''Find the minimal real positive solution of the equation'''
          if determinant <= 0:
              return np.Inf
          t = [(-b+np.sqrt(determinant))/(2*a), (-b-np.sqrt(determinant))/(2*a)]
          if t[1] >= 0:
              return t[1]
          elif t[0] >= 0:
              return t[0]
          else:
              return np.Inf


      def collide(self, other):     
          '''Normal direction of the contacting surface'''
          normal = self.__pos - other.__pos

          '''Find velocity component parallel to the normal direction'''
          vel1_para = normal * float(np.dot(self.__vel, normal)) / np.dot(normal, normal)
          vel2_para = normal * float(np.dot(other.__vel, normal)) / np.dot(normal, normal)

          '''Parallel component of centre of mass velocity'''
          vel_cm = (vel1_para * self.__m + vel2_para * other.__m) / (self.__m + other.__m)

          '''Find new normal velocities after the collision'''
          vel1_new_para = 2 * vel_cm - vel1_para
          vel2_new_para = 2 * vel_cm - vel2_para

          '''Combine perpendicular conponents with the new normal components'''
          self.__vel = ( self.__vel - vel1_para ) + vel1_new_para
          other.__vel = ( other.__vel - vel2_para ) + vel2_new_para

          return self, other







class container(particle):

    def __init__(self, mass, R, pos = [0,0], vel = [0,0]):
        particle.__init__(self, mass, -R, pos, vel)
        self.__R = R
    

    def __repr__(self):
        return "Container (Mass: %r, Radius: %r)" %(self.__m, -self.__R)
    
    
    def radius(self):
        return self.__R

    
    def perimeter(self):
        return 2 * np.pi * self.__R


    def surface_area(self):
        return 4 * np.pi * self.__R**2
        






class gas(container):
      '''Record the soonest time when a collision will happen'''
      min_time = 0
      
      '''Find which two particles are involved in the next collision'''
      colliding_particles = []
      
      '''The total momentum change of particles after they collied with the container'''
      accumulated_momentum_change = 0
      
      '''Record the speeds of particles'''
      spd = []
      
      '''Record the pressure in the container'''
      pressure = 0
      
      '''Record the total momentum of all objects at a moment'''
      momentum = 0
      
      '''Record the total kinetic energy of all objects at a moment'''
      ke = 0
      


      def __init__(self, particle_number, particle_m, particle_r, average_ke, container, dt = 0.001,
                   dimension = '2d', spd_distribution = False, record_stats = False):
          self.__number_of_particles = particle_number
          self.__particle_m = particle_m
          self.__particle_r = particle_r
          self.__ke_average = average_ke
          self.__container1 = container
          self.__dt = dt
          
          '''Use to record velocity distribution between frames'''
          self.__spd_distribution = spd_distribution
          self.__record_stats = record_stats
          self.__text0 = None
          
          '''a list that stores all particle objects'''
          self.__particle_objects = []
          
          '''Check input dimension'''
          if dimension in ['2d', '3d', '2D', '3D', 2, 3]:
              self.__d = int(str(dimension)[0])
          else:
              raise ValueError("Wrong dimension!")
                       
          self.__Random_Particles()
          self.__all_objects = [*self.__particle_objects, self.__container1]
         
            
          
      '''This method is used to randomly generate positions of particles.'''
      def __Random_Particles(self):
          
          '''Calculate how many particles can be aligned along the radius of the container.
             The coefficient 1.1 makes space between particles.'''
          particle_diameter = 2 * self.__particle_r
          p = int(self.__container1.radius() / (1.1 * particle_diameter))
          
          '''Prevent too many particles being put in the container.'''
          if (4/3) * np.pi * p**self.__d < self.__number_of_particles:
              raise Exception("Too many particles in the container!")
          
          '''a list that stores occupied positions.'''
          location_occupied = [] 
          for i in range(self.__number_of_particles):
              loop = 0
              while True:
                  pos = np.array([randint(-p,p) * particle_diameter * 1.1 for i in range(self.__d)])
                  
                  '''Set particle postions with respect to the position of the container.'''
                  pos = list(pos + self.__container1.pos())
                  
                  '''Distance away from the center of the container'''
                  distance = np.sqrt(np.dot(pos - self.__container1.pos(), pos - self.__container1.pos()))
                  
                  if pos not in location_occupied and distance < self.__container1.radius() - self.__particle_r:
                      location_occupied.append(pos)
                      break
                  loop += 1
                  '''Prevent dead or extremely long loop'''
                  if loop == 10 * self.__number_of_particles:
                     raise Exception("Too many particles in the container!")
                      
              
              '''Set each particles' kinetic energy equal to the average kinetic energy'''
              average_vel = np.sqrt(2 * self.__ke_average / self.__particle_m)
              vel = [(-1 + 2 * randint(0,1)) * average_vel] + [0 for i in range(self.__d - 1)]
              
              '''Put the newly generated particle in the list'''
              self.__particle_objects.append(particle(self.__particle_m, self.__particle_r, pos, vel))
              
          return None



      '''Find the soonest time at which a collision will happen and which two particles are involved.
         The results will be assigned to class variables.'''
      def __next_collision(self):
          objects = self.__all_objects
          collision = []

          '''Find the next collision time of every pair of objects.'''
          for x, i in enumerate(objects):
              for j in objects[x+1::]:
                  collision.append([i, j, i.time_to_collision(j)])

          '''Assign the shortest time interval to gas.min_time and
             the corresponding pair of objects to gas.colliding_particles'''
          gas.min_time = min(list(zip(*collision))[2])
          colliding_particles = [i[0:2] for i in collision if i[2] == gas.min_time][0]
          gas.colliding_particles = colliding_particles
          
          return None



      def init_figure(self):
        '''Initialise the text that shows the framenumber.'''
        self.__text0 = ax.text(-19.5, 18, "f={:4d}".format(0, fontsize = 12))
        patches = [self.__text0]
        
        '''Add the patches of the particles to the plot.'''
        for b in self.__all_objects:
            pch = b.get_patch()
            ax.add_patch(pch)
            patches.append(pch)
        
        return patches



      def next_frame(self, framenumber):
          try:
              '''This line does not work if animation method is not called.'''
              self.__text0.set_text("f={:4d}".format(framenumber))
          except:
              pass
          
          patches = [self.__text0]
    
          '''The time step for each frame'''
          dt = self.__dt
    
          '''Record for how long the particles will be moving in the following while-loop in one frame.''' 
          time_in_while_loop = 0
    
          '''Record the next collision time and colliding particles'''
          self.__next_collision()
    
          '''Continuously check whether there will be a collision in less than dt.'''
          while gas.min_time < self.__dt:
                  
              '''Although the actual velocity of each particle is constant when no collisions
                 are taking place, the change of time step among frames makes VISUAL velocities 
                 vary in the animation. In order to keep the time step constant, particles can't be 
                 allowed to move for a time interval longer than self.__dt in this while-loop. 
                 After the while-loop is broken, the particles move for (self.__dt - time_in_while_loop),
                 so the total moving time in each frame is always constant.'''
              each_time = 0.99 * gas.min_time
              time_in_while_loop += each_time
              if time_in_while_loop < self.__dt:
                 dt = self.__dt - time_in_while_loop
              else:
                 break
    
              '''Move every particle to the position corresponding to the time at which 
                 the next collision will take place.
                 The coefficient (0.99) in "each_time" makes sure that two particles are not 
                 progressing too much (due to the rounding of gas.min_time) to overlap.'''
              for b in self.__all_objects:
                  b.move(each_time)
    
              '''Get new velocities of the particles(or a particle and a container) that collide.
                 Before new velocities are assigned, the previous velocity of colliding_particles[0], which 
                 must be a particle object, is recorded and the momentum change is calculated if 
                 colliding_particles[1] is a container. It is used to calculated the pressure in the container.'''
              vel_before_collision = gas.colliding_particles[0].vel()
              gas.colliding_particles[0].collide(gas.colliding_particles[1]) # change velocities
              if gas.colliding_particles[1] == self.__container1:
                 gas.accumulated_momentum_change += gas.colliding_particles[0].momentum_change(vel_before_collision)
    
              '''Find the next collision and update gas.colliding_particles and gas.min_time.'''
              self.__next_collision()
    
          '''Move the particles again and update patches of the particles.
             Note that dt = self.__dt - time_in_while_loop, so the total moving time
             in this frame is self.__dt.'''
          for b in self.__all_objects:
              b.move(dt)
              patches.append(b.get_patch())
              
               
          '''The next collision time interval is reduced by self.__dt at the end of each frame.'''
          gas.min_time -= self.__dt
          
          

          '''The velocities of particles within the chosen frames will be recorded.'''       
          if self.__spd_distribution != False and self.__spd_distribution[0] <= framenumber < self.__spd_distribution[1]:
              for b in self.__all_objects[0:-1]:
                  gas.spd.append(np.sqrt(np.dot(b.vel(), b.vel())))
          
             
          if self.__record_stats == True:
              '''Calculate the pressure in the container.'''
              x = [None, None, self.__container1.perimeter(), self.__container1.surface_area()]
              gas.pressure = gas.accumulated_momentum_change / (self.__dt * (framenumber + 1) * x[self.__d])
          
          
              '''Calculate the total momentum and kinetic energy in each frame.'''
              gas.momentum = [0 for i in range(self.__d)] # initialisation
              gas.ke = 0 # initialisation
              for x in self.__all_objects:
                  gas.momentum += x.momentum()
                  gas.ke += x.KE() 
              
          return patches
      
        
        
      '''Plot the histogram of speed distribution and its best-fit function'''
      def histogram(self, bin_number, guess = 300):
          plt.figure()
          
          counts, bins, patch = plt.hist(gas.spd, bins = bin_number, label = "Simulation Result", normed = True, color = "skyblue", ec = "skyblue")
          step = (bins[1] - bins[0]) / 2
          bins = [ i + step for i in bins[:-1] ]
          k = 1.38e-23
          m = self.__particle_m 
          
          '''This is normalised 2D Maxwell-Boltzmann distribution equation.'''
          def boltzmann_2d(p, x):
              T = p
              return (m/T/k)*x*np.exp((-m*x**2.)/2./k/T)
          
          '''This is unnormalised 3D Maxwell-Boltzmann distribution equation.'''
          def boltzmann_3d(p, x):
              T = p
              return ((m/(2.*np.pi*k*T))**1.5)*4*np.pi*(x**2)*np.exp((-m*x**2.)/2./k/T)
          
          func_list = [None, None, boltzmann_2d, boltzmann_3d]
          
          '''Use optimal distance regression method to fit the histogram.
             The estimation must be close to the real value and cannot be too large.'''
          init_guess = [guess]
          model = scipy.odr.Model(func_list[self.__d])
          data = scipy.odr.RealData(bins, counts)
          odr = scipy.odr.ODR(data, model, beta0 = init_guess)
          out = odr.run()
          out.pprint()
            
          x_fit = np.linspace(0, bins[-1], 2000)
          y_fit = func_list[self.__d](out.beta, x_fit)
          fit_label = [None, None, "Unnormalised 2D Maxwell-\nBoltzmann Distribution Fit", "Unnormalised 3D Maxwell-\nBoltzmann Distribution Fit"]
          plt.plot(x_fit, y_fit, c = 'b', label = fit_label[self.__d])
          plt.title("Probability Density   v.s.   Particle Speed")
          plt.xlabel("Particle Speed (m/s)")
          plt.ylabel("Probability Density")
          plt.scatter(bins, counts, c ='k', s = 4)
          plt.legend(loc='best')
          plt.show()
          
          return None
      
        
        
      def Pressure(self):
          return gas.pressure
      
        
      def Momentum(self):
          return gas.momentum
      
        
      def KE(self):
          return gas.ke
          



if __name__ == "__main__":
    '''gas(particle_number, particle_m, particle_r, average_ke, container, dt, dimension, spd_distribution, record_stats)
       Position and velocity of the container and the dimension parameter can be set to 3D.'''
    movie = gas(10, 0.028/(6.02e23), 1, 6.213e-21, container(9999999999999, 15, pos = [0,0], vel = [0,0]),
                dt = 0.001, dimension = 2, spd_distribution = False, record_stats = False)
    
    
    '''The following is the animation method.'''
    fig = plt.figure()
    ax = plt.axes(xlim = (-20, 20), ylim = (-20, 20))
    ax.axes.set_aspect('equal')
    
    
    anim = animation.FuncAnimation( fig,
                                    movie.next_frame,  
                                    interval = 20,
                                    init_func = movie.init_figure, 
                                    blit = True)
    plt.show()
    

    '''
    # The animation part can be commented and the following code can be used to get statistics.
    # To get speed distribution, set spd_distribution to be [start frame, end frame].
    # To get pressure, total kinetic energy and momentum, set record_stats to be True.
    for i in range(5000):
        movie.next_frame(i)
        
        if i % 500 == 0:
           print ("The pressure is " + str(movie.Pressure()))
           print ("The total momentum is " + str(movie.Momentum()))
           print ("The total KE is " + str(movie.KE()))   
        
        
    movie.histogram(bin_number = 30, guess = 300)
   '''
                
