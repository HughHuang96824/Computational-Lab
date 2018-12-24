The program is a simulation of gas molecule in a container.

The gas-container system is an instance of the "gas" class. An instance needs the following parameters:
1. Number of particles
2. Mass of a single particle
3. Radius of a particle
4. Mean kinetic energy of particles
5. Container object (e.g. container(mass, radius, position, velocity))
6. dt, which specifies the time lapse in between two frames.
7. Dimension (must be either 2 or 3)
8. Speed distribution Range. If a histogram of speed distribution is needed, set spd_distribution to be a list of two elements, where the first and second elements denote the range of frames within which velocities of particles will be recorded. Otherwise set it to be False.
9.  record_stats. User will be able to get the pressure, total momentum and kinetic energy of the system if this parameter is set to be True

e.g.
movie = gas(30, 1, 1, 100000, container(9999999999999, 10, pos = [0,0,0], vel = [0,0,0]),
                dt = 0.001, dimension = 3, spd_distribution = False, record_stats = True)

Things to notice: the position and velocity of the container must have the same dimension as indicated by the dimension parameter.



Run the program


I. Use animation

1.) If only the animation method is called, the pressure, total momentum and kinetic energy of the
     system and the speed-distribution histogram will not be showed.
2.) 3D gas can also be animated, but motions in the third dimension are not observable, which makes
     the animation look weird.


II. Get statistics

1.) Animation method should be commented.
2.) Use a for-loop to loop over the next_frame method. Print the statistics the user wants to see. 
     e.g. print ("The pressure is " + str(instance.Pressure()))
            print ("The total momentum is " + str(instance.Momentum()))
            print ("The total KE is " + str(instance.KE())) 
3.) use instance.histogram(bin_number, guess) to show the histogram. The number of bins can be
     customised. The histogram will come out fitted with the Maxwell-Boltzmann distribution function.
     The fitting function has temperature as the undetermined parameter. The user must estimate the
     value of the temperature, which cannot be too far off or too large (e.g. 5e29) or the program will
     not be able to fit the histogram.
