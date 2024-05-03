# -*- coding: utf-8 -*-

######### PROGRAM SYNTHETIC PHENOLOGY MAX ENTROPY ###########
# This programs generates synthetic phenological configurations under certain restrictions
# We solve a multivariate system of inequations to determine the phase space of phenological configurations compatible with the matrix of interactions 
# We solve the system using linear programming and a objective function to minimize/maximize
# In this case we maximize the entropy associated to the fraction of species sharing the same starting date

import sys
import numpy as np
from scipy.optimize import minimize, Bounds, LinearConstraint, differential_evolution, shgo
import matplotlib.pyplot as plt
from scipy import stats
import warnings 
import datetime

###############################################
############### FUNCTIONS #####################
###############################################

###########  ENTROPY OF THE MIDDLE DATES ##################
## Calculates the negative Shannon entropy of the number of species with initial date = ti
def entropy(x, n_total, lbound, ubound, periods):
      
      #Calculating the vector of middle dates
      tm = np.empty(shape=(n_total))
      for i in range(n_total):
            tmiddle = x[i] + (float(periods[i])-1)/2.0
            tm[i] = tmiddle
      

      #Create a sequence of the bins (1 day bin width)
      seq = np.linspace(lbound,ubound,(ubound-lbound)+1)

      #Count how many species share the same starting date
      values, bins = np.histogram(tm,seq)

      #Compute the Shannon entropy
      ent = 0
      for xi in values:
            if xi > 0:
                  #Calculate fraction
                  pi=float(xi)/float(n_total)
                  #Sum to the entropy
                  ent = ent - pi*np.log(pi)

      #We take the negative, so that the minimization is actually a maximization
      ent = -ent

      return ent

def print_convergence(x,convergence): #definition for the minimizers 'minimize' and 'differential_evolution'

      global lbound, ubound, periods

      #Calculating the vector of middle dates
      tm = np.empty(shape=(n_total))
      for i in range(n_total):
            tmiddle = x[i] + (float(periods[i])-1)/2.0
            tm[i] = tmiddle

      #Create a sequence of the bins (1 day bin width)
      seq = np.linspace(lbound,ubound,(ubound-lbound)+1)

      #Count how many species share the same starting date
      values, bins = np.histogram(tm,seq)

      #Compute the Shannon entropy
      ent = 0
      for xi in values:
            if xi > 0:
                  #Calculate fraction
                  pi=float(xi)/float(n_total)
                  #Sum to the entropy
                  ent = ent - pi*np.log(pi)

      ent = -ent
      #We take the negative, so that the minimization is actually a maximization
      print ("entropy=",ent)
      print ("converged fraction=", convergence)



###############################################
######## INITIALIZATION OF VARIABLES ##########
###############################################


#reading input argument: index of the network to study. 2 for the Illinois datase, 6 for the Kantsa dataset 1st year and 7 for the Kantsa dataset 2nd year
idoc = sys.argv[1]

#opening files
doc_dimension = open ("general"+idoc+".dat","r")
doc_matrix = open ("matrix"+idoc+".dat","r")

#reading dimension (number of rows, number of columns)
n_rows = int(doc_dimension.readline())
n_cols = int(doc_dimension.readline())
n_total = n_rows + n_cols
doc_dimension.close()

#reading interaction matrix
matrix = []
matrix = [[int(num) for num in line.split('\t')] for line in doc_matrix]
doc_matrix.close()

#calculating the number of links L
L=0
for i in range(n_rows):
      for j in range(n_cols):
            L = L + matrix[i][j]



###############################################
###### GENERATION OF SYNTHETIC PERIODS ########
###############################################

#Reading periods of activity from file 
doc_periods = open ("periods"+idoc+".dat","r")
periods = []
periods_rows = []
for i in range(n_rows):
      p = int(doc_periods.readline())
      periods_rows.append(p)
      periods.append(p)

periods_cols = []
for j in range(n_cols):
      p = int(doc_periods.readline())
      periods_cols.append(p)
      periods.append(p)




###############################################
####### GENERATION OF STARTING DATES ##########
###############################################


###### Formating the constraints using the LienarConstrains object
#### Inequalities must be written in the standard form A*x <= b. We have L inequalities, one per link
A = np.zeros(shape=(L,n_total)) #2-D array of coefficients for x. Size: LxN
lb = np.zeros(shape=(L)) #vector of lower bounds, size L
ub = np.zeros(shape=(L)) #vector of upper bounds, size L

#### Construction of the inequalities 
iiter = 0
for i in range(n_rows):
      for j in range(n_cols):                                 
            #If there's a mutualistic link, construct the two corresponding inequalities
            if (matrix[i][j] == 1):
                  #Coefficients for the starting dates
                  A[iiter][i] = 1
                  A[iiter][j+n_rows] = -1
                  #Lower bound
                  lb[iiter] = -periods_rows[i] + 1
                  #Upper bound
                  ub[iiter] = periods_cols[j] - 1
                  iiter = iiter + 1



####Formating the constraints using a dictionary
cons_dict = [] #sequence of dictionaries
for i in range(n_rows):
      for j in range(n_cols):                                 
            #If there's a mutualistic link, construct the two corresponding inequalities
            if (matrix[i][j] == 1):
                  #Lower bound
                  cons_dict.append({'type': 'ineq', 'fun': lambda x : -x[j+n_rows]+periods_rows[i]-1+x[i]})
                  #Upper bound
                  cons_dict.append({'type': 'ineq', 'fun': lambda x : x[j+n_rows]+periods_cols[j]-1-x[i]})


#### Setting the lower and upper bounds (limits of the season). We constrain the starting dates to coincide within one month.
lbound = 100
ubound = 200
veclb = np.empty(shape=(n_total))
vecub = np.empty(shape=(n_total))
for i in range(n_rows):
      veclb[i] = lbound 
      vecub[i] = ubound - periods_rows[i]+1
      
for j in range(n_cols):
      veclb[j+n_rows] = lbound 
      vecub[j+n_rows] = ubound - periods_cols[j]+1
season_bounds = Bounds(veclb,vecub) #Valid for the methods 'minimize' and 'differential_evolution'
seq_bounds = np.column_stack((np.array(veclb),np.array(vecub)))
seq_bounds = seq_bounds.tolist()


#### We solve the system of inequalities by linear programming


########Local+global minimization using the functions 'minimize' and 'differential_evolution'
Nsample = 15
#We start constructing a population by solving iteratively the problem with a local minimizer, then use each of these solutions as a initializer
init_pop = []
for i in range(Nsample):
      #Draw an initial point
      t0 = np.random.uniform(lbound,ubound,size=(n_total)) 
      #Maximize the entropy with the function 'minimize
      res = minimize(entropy, t0, args = (n_total,lbound,ubound,periods), method='trust-constr', bounds=season_bounds, constraints=(LinearConstraint(A, lb, ub, keep_feasible=False)), tol=0.1,callback=None, options={'disp':False})
      #Save the initial times
      ti = res.x
      init_pop.append(ti)

init_pop = np.reshape(init_pop,(Nsample,n_total))

print(init_pop)

#Maximize the entropy using the global optimizer 'differential_evolution' and the local solutions as input
res=differential_evolution(entropy, bounds=season_bounds, args=(n_total, lbound, ubound,periods), strategy='best1exp', maxiter=1000000000000, tol=0.1,disp=False, polish=True, init=init_pop,constraints=(LinearConstraint(A, lb, ub,keep_feasible=False)))

#Print results
ti = res.x
print ("Solution of starting dates=", ti)
print("Results=",res)

#Initial times for rows
ti_rows = []
for i in range(n_rows):
      ti_rows.append(int(round(ti[i])))      

#Initial times for columns
ti_cols = []
for i in range(n_cols):
      ti_cols.append(int(round(ti[i+n_rows])))



###############################################
###### RESULTS: Plots and writing files########
###############################################

####Construct matrices of presence/absence
#Minimum initial date
min_ti_rows = min(ti_rows)
min_ti_cols = min(ti_cols)
min_ti = min(min_ti_rows, min_ti_cols)

#Ending date
tf_rows = []
for i in range(n_rows):
      tf = ti_rows[i] + periods_rows[i]
      tf_rows.append(tf)

tf_cols = []
for i in range(n_cols):
      tf = ti_cols[i] + periods_cols[i]
      tf_cols.append(tf)

#Maximum final date
max_tf_rows = max(tf_rows)
max_tf_cols = max(tf_cols)
max_tf = max(max_tf_rows, max_tf_cols)

#Matrix of presence for rows
matrix_rows = []
for i in range(n_rows):
      matrix_rows.append([])
      ti = ti_rows[i]
      tf = tf_rows[i]
      for t in range (min_ti, max_tf):
            #within the period of activity, append a 1
            if (ti <= t < tf):
                  matrix_rows[i].append(1)
            #outside the period of activity, append a 0
            else:
                  matrix_rows[i].append(0)

#Matrix of presence for columns
matrix_cols = []
for i in range(n_cols):
      matrix_cols.append([])
      ti = ti_cols[i]
      tf = tf_cols[i]
      for t in range (min_ti, max_tf):
            #within the period of activity, append a 1
            if (ti <= t < tf):
                  matrix_cols[i].append(1)
            #outside the period of activity, append a 0
            else:
                  matrix_cols[i].append(0)


######## Writing the periods and starting dates on a file

#Opening files
doc_periods = open ("maxent_periods"+idoc+".dat","w+")
doc_ti = open ("maxent_ti"+idoc+".dat","w+")

#Writing periods
for i in range(n_total):
      doc_periods.write(str(periods[i])+"\n")

#Writing initial times
for i in range(n_rows):
      doc_ti.write(str(ti_rows[i])+"\n")
for i in range(n_cols):
      doc_ti.write(str(ti_cols[i])+"\n")


#####Drawing figure: we plot the timeline of periods, distinguishing between plants (rows) and pollinators (columns)
fig, axs = plt.subplots(2)
axs[0].imshow(matrix_rows,cmap='Greens', aspect='auto')
axs[0].set(ylabel='Plant species')
axs[0].xaxis.set_ticks_position('none') 
axs[0].set_xticks([]) 
plt.xlabel('Day of the year')
plt.ylabel('Pollinator species')
axs[1].imshow(matrix_cols,cmap='BuPu', aspect='auto')
plt.draw()
plt.show()
fig.savefig('maxent_periods'+idoc+'.png')
plt.clf()


