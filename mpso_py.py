import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import random

def schaffer_f6(particle):
	# particle - 2 dim numpy vect
	xy2 = particle[0]**2 + particle[1]**2
	return 0.5 + (np.sin(np.sqrt(xy2))**2 - 0.5) / ((1 + 0.001 * xy2)**2)

def plot_func():
	N = 201
	x = np.linspace(-100, 100, N)
	y = np.linspace(-100, 100, N)
	X, Y = np.meshgrid(x, y)

	Z = np.zeros((N, N))

	for i in range(len(x)):
		for j in range(len(y)):
			Z[i, j] = schaffer_f6([x[i], y[j]])

	#fig = plt.figure()
	#ax = fig.add_subplot()
	#ax.contour(X, Y, Z)
	#plt.grid()
	#plt.show()
	
	#plt.plot(Z[:, N//2])
	#plt.show()		

	#fig = plt.figure()
	#ax = plt.axes(projection='3d')
	#ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
        #        cmap='winter', edgecolor='none')
	#plt.show()


def run_mpso(X_init, V_init, w, verbose=1):
	
	X = X_init.copy()
	V = V_init.copy()
	found_min = False
	n = 0
	while(n < max_iters):
		
		if n==0:
			X_best = np.copy(X)
			fitness_best = []
			for particle in X_best:
				fitness_best.append(schaffer_f6(particle))
			fitness_best = np.array(fitness_best)
		else:
			fitness = []
			for particle in X:
				fitness.append(schaffer_f6(particle))
			fitness = np.array(fitness)
			mask = fitness < fitness_best
			X_best[mask] = X[mask]
			fitness_best[mask] = fitness[mask]

		global_best_idx = np.argmin(fitness_best)
		# check if we have a minimum
		if min(fitness_best) < tol:
			if verbose:
				print("Minimum found in {0} iterations \n".format(n))
				print("Particle : ", X[global_best_idx])
				print("Fitness value : ", min(fitness_best))
			found_min = True
			break
		
		# update particles according to local best and global best		
		for k in range(pop_size):
			V[k, :] = w * V[k, :] + c1 * np.random.uniform() * (X_best[k, :] - X[k, :]) + c2 * np.random.uniform() * (X_best[global_best_idx, :] - X[k, :])
			# check for boundries
			if V[k, 0] > 2:
				V[k, 0] = 2
			if V[k, 1] > 2:
				V[k, 1] = 2
			# update positon			
			X[k, :] = X[k, :] + V[k, :]
			# check for boundries
			if X[k, 0] > 100:
				X[k, 0] = 100
			if X[k, 0] < -100:
				X[k, 0] = -100
			if X[k, 1] > 100:
				X[k, 1] = 100
			if X[k, 1] < -100:
				X[k, 1] = -100
			

		n += 1
	
	if not found_min:
		if verbose:
			print("\nMinimul nu a fost gasit.\n")
		return 0, -1
	else:
		return 1, n


if __name__ == "__main__":
	
	pop_size = 20
	max_velocity = 2
	max_iters = 4000
	c1 = 2
	c2 = 2
	tol = 0.0001
	ALL_ITERS = []	
	
	for exp in range(30):
		
		print("Experiment : {0}/30".format(exp+1))
		# init particle's position
		X_init = np.random.uniform(low=-100, high=100, size=(pop_size, 2))
		# init particle's velocity
		V_init = np.random.uniform(low=0, high=2, size=(pop_size, 2))

		total_iters = []
      
		for w in [1.4, 1.2, 1.1, 1.05, 1, 0.95, 0.9, 0.85, 0.8, 0]:
			found, n_iters = run_mpso(X_init, V_init, w, verbose=0)
			total_iters.append(n_iters)

		
		print(total_iters)
		ALL_ITERS.append(total_iters)

	ALL_ITERS = np.array(ALL_ITERS)
	np.savetxt("all_experiments", ALL_ITERS, delimiter=' ', fmt="%.2f")
	
    """
	exps = np.loadtxt("all_experiments", delimiter=' ')
	print("Failures : \n")	
	print(np.sum(exps == -1, axis=0))
    """






	
