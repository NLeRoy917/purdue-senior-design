"""
This python script serves to model and calculate the powder flow characteistics of our
product-excipient formulation after they have been mixed together.

The equations and theory for this script come from the following journal article by
Sklubalova et. al.

Šklubalová, Z., & Zatloukal, Z. (2013). Flow rate and flow equation of pharmaceutical free-flowable powder excipients. 
Pharmaceutical Development and Technology, 18(1), 106–111. https://doi.org/10.3109/10837450.2011.640686
DOI: 10.3109/10837450.2011.640686

Script designed and written by: Nathan LeRoy
"""


#import necessary modules
import numpy as np
import sys
import matplotlib.pyplot as plt


#define constants
g = 981 #cm/s/s
pi = np.pi #pi
rho = 1.25 # g/cm^3
n = 0.4

#Initiate lists
flow_rates = np.logspace(1,1000,100) #logarithmically spaced vector of flow rates


#calculate orififace sizes
oriface_sizes = [((4/pi)*(Q/(rho*np.sqrt(g))))**(1/n) for Q in flow_rates]

