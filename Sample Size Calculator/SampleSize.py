# -*- coding: utf-8 -*-
"""
Created on Sat Aug 28 10:11:39 2021

@author: jddrumheller
"""

## create a class to compute sample sizes for four different tests

## IMPORT LIBRARIES

import numpy as np
from scipy.stats import norm, t

class SampleSize:
    
    '''
    A class to hold/compute sample size information. 
    
    Attributes:
        mu1: mean/proportion 1
        mu2: mean/proportion 2
        sigma: common variance
    '''
    
    def __init__(self, mu1, mu2, sigma):
        self.mu1 = mu1
        self.mu2 = mu2
        self.sigma = sigma

        
        
    ## CLASS METHODS ## 
    
    ## z-test sample size
    
    def zSampleSize(self, test = 'equal', numSim = 2000, maxN = 500, alpha = 0.05, power = 0.8):
        
        """ Compute the sample size for a test to a known mean parameter
        muHyp: the hypothesized mean parameter -user must supply
        sigma: a known population variance -user must supply
        
        test: type of test on parameter eg. one-sided or equal
        mu: the population mean parameter with default 0
        numSim: the number of simulations to be used with default 2000
        maxN: maxium number of samples with default 500
        alpha: the prescribed test size with default value 0.05
        power: the prescribed power of the test with default value 0.8
        
        """
        
        ## get critical value based on test type
        
        if test == 'equal':
            criticalValue = np.abs(norm.ppf(alpha/2))
        else:
            criticalValue = np.abs(norm.ppf(1-alpha))
                    
        ## loop over sample size
        
        for i in range(2, numSim):
            
            ## set a counter to see how many times we reject the null
            rejectedCounter = 0
            
            for j in range(numSim):
                
                ## simulate normal random variable and compute z-statistic
                simulated = np.random.normal(self.mu1, self.sigma, i)
                zStatistic = (simulated.mean() - self.mu2) / (self.sigma / np.sqrt(i))
                
                ## check if the statistic exceeds the specified critical value and increment counter
                if test == 'equal':
                    if np.abs(zStatistic) > criticalValue:
                        rejectedCounter += 1
                else:
                    if np.sign(self.mu1 - self.mu2) > 0:
                        if zStatistic > criticalValue:
                            rejectedCounter += 1
                    else:
                        if zStatistic < -criticalValue:
                            rejectedCounter += 1
                            
            ## check if we achived sufficient power
            if rejectedCounter / numSim >= power:
                return(i)
                break
            
                           
                            
    ## Independent Samples t-test sample-sie ##                
    
    def tSampleSize(self, test = 'equal', numSim = 2000, maxN = 500, alpha = 0.05, power = 0.8):
        
        """ Compute the sample size for a indpendent sample t-test
        mu1: the hypothesized mean parameter for group 1 -user must supply
        mu2: the hypothesized mean parameter for group 2 -user must supply
        sigma: a known population variance -user must supply
        
        test: type of test on parameter eg. one-sided or equal
        mu: the population mean parameter with default 0
        numSim: the number of simulations to be used with default 2000
        maxN: maxium number of samples with default 500
        alpha: the prescribed test size with default value 0.05
        power: the prescribed power of the test with default value 0.8 
        """
        
        ## loop over sample size
        for i in range(3, maxN):
            
            ## determine number of times we reject
            rejectedCounter = 0
            
            ## get critical value based on test type             
            if test == 'equal':
                criticalValue = np.abs(t.ppf(alpha/2, 2*i -2))
            else:
                criticalValue = np.abs(t.ppf(alpha, 2*i-2))
                
            ## create simulations of data and test-statistic
            
            for j in range(numSim):
                
                ## simulate from two normal distributions
                simulated1 = np.random.normal(self.mu1, self.sigma, i)
                simulated2 = np.random.normal(self.mu2, self.sigma, i)
                
                ## calculate test statistic
                muHat1 = simulated1.mean()
                muHat2 = simulated2.mean()
                S1 = simulated1.var()
                S2 = simulated2.var()
                varPooled = (S1 + S2) / 2
                tStatistic = (muHat1 - muHat2) / np.sqrt(varPooled * (2/i))
                
                ## check if t-statistic exceed critical value an increment counter
                
                if test == 'equal':
                    if np.abs(tStatistic) > criticalValue:
                        rejectedCounter += 1
                else:
                    if np.sign(self.mu1 - self.mu2) > 0:
                        if tStatistic > criticalValue:
                            rejectedCounter += 1
                        else:
                            if tStatistic < -criticalValue:
                                rejectedCounter += 1
                
                ## check if we've achived sufficient power
                if rejectedCounter / numSim > power:
                    return(i)
                    break
        
    ## Test of a single proportion sample size
    
    def PropSampleSize(self, test = 'equal', numSim = 2000, maxN = 500, alpha = 0.05, power = 0.8):
         
        """ Compute the sample size for a test to a known proportion
        pHyp: the hypothesized proportion value -user must supply
        pKnw: the population proportion parameter - user must supply
        
        test: type of test on parameter eg. one-sided or equal
        numSim: the number of simulations to be used with default 2000
        maxN: maxium number of samples with default 500
        alpha: the prescribed test size with default value 0.05
        power: the prescribed power of the test with default value 0.8
        
        """
        
        ## check type of test and get critical value
        
        if test == 'equal':
            criticalValue = np.abs(norm.ppf(alpha/2))
        else:
            criticalValue = np.abs(norm.ppf(alpha))
            
        ## loop over the sample size
        
        for i in range(2, maxN): 
            
            ## counter to determine number of times test statistic exceed the critical value
            rejectedCounter = 0
            
            ## create simulations based on user input
            for j in range(numSim):
                
                ## simulated a Bernoulli random variable and compute z-statistc
                simulated = np.random.binomial(1, self.mu1, i)
                phat = simulated.mean()
                
                ## get value for continuity correction
                if self.mu2 > phat: 
                    c = 1 / (2*i)
                elif self.mu2 < phat:
                    c = -1 / (2*i)
                elif np.abs(self.mu2 - phat) < 1 / (2*i):
                    c = 0
                    
                zStatistic = ((phat - self.mu2) + c) / np.sqrt((self.mu2 * (1 - self.mu2)) / i)
                
                ## check if test statistic exceeds the critical value
                
                if test == 'equal':
                    if np.abs(zStatistic) > criticalValue:
                        rejectedCounter += 1
                else:
                    if np.sign(phat - self.mu2) > 0:
                        if zStatistic > criticalValue:
                            rejectedCounter += 1
                    else:
                        if zStatistic < -criticalValue:
                            rejectedCounter += 1
                
                ## check if we've achived suffcient power
                
                if rejectedCounter / numSim >= power:
                    return(i)
                    break
                
                
    ## Test between two sample proportions sample size
    
    def TwoPropSampleSize(self, test = 'equal', numSim = 2000, maxN = 500, alpha = 0.05, power = 0.8):
        
        ## check for type of test and get critical value
        
        if test == 'equal':
            criticalValue = np.abs(norm.ppf(alpha/2))
        else:
            criticalValue = np.abs(norm.ppf(alpha))
            
        ## loop over sample size
        
        for i in range(4, maxN):
            
            ## counter to determine number of times the null is rejected
            rejectedCounter = 0
            
            ## create simulations based on input
            
            for j in range(numSim):
                
                ## simulate two Bernoulli random variables and compute z-statistic
                simulated1 = np.random.binomial(1, self.mu1, i)
                simulated2 = np.random.binomial(1, self.mu2, i)
                phat1 = simulated1.mean()
                phat2 = simulated2.mean()
                
                ## helper terms
                pnet = np.concatenate([simulated1, simulated2]).mean()                
                denom = np.sqrt(pnet * (1 - pnet) * 2 / i)
                
                ## check for too small denominator and add some noise
                if denom < 0.0001:
                    denom = denom + 0.0001
                
                zStatistic = (phat1 - phat2) / denom
                
                ## check if statistic exceeds the critical value
                
                if test == 'equal':
                    if np.abs(zStatistic) > criticalValue:
                        rejectedCounter += 1
                else:
                    if np.sign(self.mu1, self.mu2) > 0:
                        if zStatistic > criticalValue:
                            rejectedCounter += 1
                            
                    else:
                        if zStatistic < -criticalValue:
                            rejectedCounter += 1
                
                ## check if we've achieved sufficient power
                
                if rejectedCounter / numSim >= power:
                    return(i)
                    break
                
            
##############################################################################   
 ## main function for testing
 

def functionParameters():
    
    ''' Function to get the user input for power calculation specification
    
    no inputs
    returns: the type of test (one-sided/two-sided), alpha/test size, and power
    as specified by the user
    
    '''
    
    test_type = input('One sided [1] or two sided [2] test: ')
    alpha_input = float(input('Alpha level: '))
    power_input = float(input('Power level: '))
            
    if test_type == '2':
        test_input = 'equal'
    else:
        test_input = 'one-sided'
        
    return(test_input, alpha_input, power_input)

def main():
    
    quitFlag = True
    
    while quitFlag:
    
        print("Welcome to Sample Size by Simulation. Please make a selection: \n \
           [1]: Compare a mean to a known value (z-test) \n \
           [2]: Compare two indpendent sample means (t-test) \n \
           [3]: Compare a proportion to a known value (z-test) \n \
           [4]: Compare two indepent proportions (z-test) \n \
           [q]: To Quit" )
              
        testType = input('Power Caclulation Type: ')
        
        
        if testType == "1":
                
            mu1_input = float(input('Hypothesized Sample Mean: '))
            mu2_input = float(input('Population Mean: '))
            sigma_input = float(input('Population standard deviation: ')) 
            test_input, alpha_input, power_input = functionParameters()
            
            SampCalc = SampleSize(mu1_input, mu2_input, sigma_input)
            SampSize = SampCalc.zSampleSize(test = test_input, alpha = alpha_input, power = power_input)

            print('You need {} observations to achieve sufficient power.'.format(SampSize))
            print('***************************************************************** \n \n')
            
            
            
        elif testType == "2":
            mu1_input = float(input('Sample Mean of Group 1: '))
            mu2_input = float(input('Sample Mean of Group 2: '))
            sigma_input = float(input('Common Standard Deviation: ')) 
            test_input, alpha_input, power_input = functionParameters()
            
            SampCalc = SampleSize(mu1_input, mu2_input, sigma_input)
            SampSize = SampCalc.tSampleSize(test = test_input, alpha = alpha_input, power = power_input)
            
            print('You need {} observations to achieve sufficient power.'.format(SampSize))
            print('************************************************************* \n \n')
            
        elif testType == "3":
            mu1_input = float(input('Hypothesized Sample Proportion: '))
            mu2_input = float(input('Population Proportion: '))
            test_input, alpha_input, power_input = functionParameters()
            sigma_input = None 
            
            SampCalc = SampleSize(mu1_input, mu2_input, sigma_input)
            SampSize = SampCalc.PropSampleSize(test = test_input, alpha = alpha_input, power = power_input)
            
            print('You need {} observations to achieve sufficient power.'.format(SampSize))
            print('***************************************************************** \n \n')
            
        elif testType == "4":
            mu1_input = float(input('Sample proportion of group 1: '))
            mu2_input = float(input('Sample proportion of group 2: '))
            test_input, alpha_input, power_input = functionParameters()
            sigma_input = None    
            SampCalc = SampleSize(mu1_input, mu2_input, sigma_input)
            SampSize = SampCalc.TwoPropSampleSize(test = test_input, alpha = alpha_input, power = power_input)
            
            print('You need {} observations to achieve sufficient power.'.format(SampSize))
            print('***************************************************************** \n \n')
            
        elif testType == "q":
            quitFlag = False
            
        else:
            print("Please make a valid selection.")

 #   mySampSize = SampleSize(mu1 = 0.5, mu2 = 0.75, sigma = None)
 #   print(mySampSize.TwoPropSampleSize(test = 'equal'))

if __name__ == "__main__":
    main()
     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
