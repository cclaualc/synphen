# -*- coding: utf-8 -*-

######### PROGRAM SYNTHETIC PHENOLOGY SYNCHRONIZED ######
# This programs generates synthetic phenological configurations under certain restrictions
# In this case we perturb the starting dates, starting from a centered situation but keeping the phenology compatible with the matrix of interactions 


###############################################
############### FUNCTIONS #####################
###############################################

######### RANDOMIZATION OF TI FOR ROWS ##################
#This function randomizes the starting dates by extracting a new date, from a uniform distribution, 
#set within the range (tinf,tsup) which warrants that interacting species have non-null overlap
def randomization_ti_rows(ti_rows,ti_cols,n_rows,n_cols,periods_rows,periods_cols,matrix):
      ti_rows = []
      for indx in range (n_rows):
            tmax = 0.0 #trial maximum date
            tmin = 10000000.0 #trial minimum date
            #we look for the maximum starting date and the minimum final date among the columns partner's
            for i in range (n_cols):
                  if (matrix[indx][i] == 1):
                        #initial time for row i
                        t1 = ti_cols[i]
                        p = periods_cols[i]
                        #final time for row i
                        t2 = t1 + p -1
                        #find the maximum initial time (latest species to start)
                        if (t1 > tmax):
                              tmax = t1
                        #find the minimum final time (first species to finish)
                        if (t2 < tmin):
                              tmin = t2
            #Upper boundary (tsup) = minimum final time
            tsup = tmin 
            #Lower boundary (tinf) =  maximum initial time - p_cols + 1
            tinf = tmax - periods_rows[indx] + 1
            #Now that we have the boundaries, we may draw the starting date for columns 'indx'
            if (tinf < tsup):
                  ti = np.random.random_integers(tinf, tsup, None)
                  ti_rows.append(ti)
            elif (tinf == tsup):
                  ti = tinf
                  ti_rows.append(ti)
            else:
                  print("uep! index=",indx)

      return ti_rows


######### RANDOMIZATION OF TI FOR COLUMNS ##################
#This function randomizes the starting dates by extracting a new date, from a uniform distribution, 
#set within the range (tinf,tsup) which warrants that interacting species have non-null overlap
def randomization_ti_cols(ti_rows,ti_cols,n_rows,n_cols,periods_rows,periods_cols,matrix):
      ti_cols = []
      for indx in range (n_cols):
            tmax = 0.0 #trial maximum date
            tmin = 10000000.0 #trial minimum date
            #we look for the maximum starting date and the minimum final date among the rows partner's
            for i in range (n_rows):
                  if (matrix[i][indx] == 1):
                        #initial time for row i
                        t1 = ti_rows[i]
                        p = periods_rows[i]
                        #final time for row i
                        t2 = t1 + p -1
                        #find the maximum initial time (latest species to start)
                        if (t1 > tmax):
                              tmax = t1
                        #find the minimum final time (first species to finish)
                        if (t2 < tmin):
                              tmin = t2
            #Upper boundary (tsup) = minimum final time
            tsup = tmin 
            #Lower boundary (tinf) =  maximum initial time - p_cols + 1
            tinf = tmax - periods_cols[indx] + 1
            #Now that we have the boundaries, we may draw the starting date for columns 'indx'
            if (tinf < tsup):
                  ti = np.random.random_integers(tinf, tsup, None)
                  ti_cols.append(ti)
            elif (tinf == tsup):
                  ti = tinf
                  ti_cols.append(ti)
            else:
                  print("uep! index=",indx)
      return ti_cols


###############################################
######## INITIALIZATION OF VARIABLES ##########
###############################################

#libraries
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
from scipy import stats


#reading input argument: index of the network to study
idoc = sys.argv[1]

#opening files
doc_dimension = open ("general"+idoc+".dat","r")
doc_matrix = open ("matrix"+idoc+".dat","r")

#reading dimension (number of rows, number of columns)
n_rows = int(doc_dimension.readline())
n_cols = int(doc_dimension.readline())
doc_dimension.close()

#reading interaction matrix
matrix = []
matrix = [[int(num) for num in line.split('\t')] for line in doc_matrix]
doc_matrix.close()

#calculating degree sequences
degree_rows = []
degree_cols = []
for i in range (n_rows):
      degree_rows.append(0)
      for j in range(n_cols):
            degree_rows[i] = degree_rows[i] + matrix[i][j]

for j in range (n_cols):
      degree_cols.append(0)
      for i in range(n_rows):
            degree_cols[j] = degree_cols[j] + matrix[i][j]

#parameters
if idoc == "2":
      distr_plants = "beta"
      distr_pollinators = "beta"
      beta_param_plants = [4.36,3.69,1.25,64.62] #parameters for Burkle
      beta_param_pollinators = [0.92,1.04,2.0,79.22] #parameters for Burkle
if idoc== "6":
      distr_plants = "lognormal"
      distr_pollinators = "exponential"
      lognorm_param_plants = [0.28,-37.40,69.85] #parameters for Kantsa 6
      expon_param_pollinators = [1.0,24.4] #parameters for Kantsa 6
if idoc=="7":
      distr_plants = "lognormal"
      distr_pollinators = "exponential"
      lognorm_param_plants = [0.42,-19.52,47.70] #parameters for Kantsa 7
      expon_param_pollinators = [1.0,32.56] #parameters for Kantsa 7



###############################################
###### MAIN: GENERATION OF PHENOLOGIES ########
###############################################

#Reading periods of activity from file 
doc_periods = open ("periods"+idoc+".dat","r")

periods_rows = []
for i in range(n_rows):
      p = int(doc_periods.readline())
      periods_rows.append(p)

periods_cols = []
for j in range(n_cols):
      p = int(doc_periods.readline())
      periods_cols.append(p)



##### Generation of starting dates ######

#vectors with starting dates
ti_rows = []
ti_cols = []

#free guild: we set randomly the starting dates for columns (pollinators)
mean = 120 #mean of the normal distribution
sigma = 0.5  #standard deviation of the normal distribution
for i in range (n_cols):
      #middle time
      tm = int(np.random.normal(mean, sigma, None))
#      tm = mean
      #print(tm)
      #starting date
      p = periods_cols[i]
      ti = int(tm - (p-1)/2.0)
      ti_cols.append(ti)

#we make 100 re-randomizations of both guilds
for niter in range(100):
      #conditioned guild: we set the starting dates for rows (plants), pseudo-randomly
      #by imposing that the species with a connection in the matrix must have non-zero overlap
      ti_rows = randomization_ti_rows(ti_rows,ti_cols,n_rows,n_cols,periods_rows,periods_cols,matrix)
      #we re-randomize now the originally free guild of columns (pollinators)
      ti_cols = randomization_ti_cols(ti_rows,ti_cols,n_rows,n_cols,periods_rows,periods_cols,matrix)


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
doc_periods = open ("cen_periods"+idoc+".dat","w+")
doc_ti = open ("cen_ti"+idoc+".dat","w+")

#Writing periods
for i in range(n_rows):
      doc_periods.write(str(periods_rows[i])+"\n")
for i in range(n_cols):
      doc_periods.write(str(periods_cols[i])+"\n")

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
#plt.draw()
#plt.show()
fig.savefig('cen_synthetic_periods'+idoc+'.png')
plt.clf()


