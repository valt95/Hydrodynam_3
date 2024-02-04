import math
import numpy as np
def rrb(x) :
    def g(y):
       if y<1 :
         a=0
       else :
         a=math.pi/4
       i=1
       while i<11:
           a=a-(a+(1/2)*math.sin(2*a)+(math.pi/2)*(1-y*(math.cos(a))**(9/4)))/(1+math.cos(2*a)+(9*math.pi/8)*y*(math.cos(a))**(5/4)*math.sin(a))
           i=i+1
       return a
    if x<0.25 :
          aa=6.087667*(0.25)**(9/4)*x**(-9/4)
    else :
          tt=g(x**(-9/8))
          aa=(9*(tt+math.sin(2*tt)/2+math.pi/2)**2)/(2*(math.cos(tt))**4*(8*(math.cos(tt))**2++9*math.tan(tt)*(tt+math.sin(2*tt)/2+math.pi/2)))
    return aa

def rrbbert(x) :
    def g(y):
       if y<1 :
         a=0
       else :
         a=math.pi/4
       i=1
       while i<11:
           a=a-(a+(1/2)*math.sin(2*a)+(math.pi/2)*(1-y*(math.cos(a))**(9/4)))/(1+math.cos(2*a)+(9*math.pi/8)*y*(math.cos(a))**(5/4)*math.sin(a))
           i=i+1
       return a
    tt=g(x**(-9/8))
    aa=(9*(tt+math.sin(2*tt)/2+math.pi/2)**2)/(2*(math.cos(tt))**4*(8*(math.cos(tt))**2++9*math.tan(tt)*(tt+math.sin(2*tt)/2+math.pi/2)))
    return aa

def mdbertsh(x) :
    def g(y):
       if y<1 :
         a=0
       else :
         a=math.pi/4
       i=1
       while i<11:
           a=a-(a+(1/2)*math.sin(2*a)+(math.pi/2)*(1-y*(math.cos(a))**(9/4)))/(1+math.cos(2*a)+(9*math.pi/8)*y*(math.cos(a))**(5/4)*math.sin(a))
           i=i+1
       return a
    bb=x**(3/4)*(math.cos(g(x**(-9/8))))**(-3/2)*(3*math.pi/4)**2
    return bb

# Функция переделанная под работу с векторами
def mdmt(x) :
    def g(y):
       # if y<1 :
       #   a=0
       # else :
       #   a=math.pi/4

       # Переделка функции под работу с векторами
       a = np.pi/4 * np.where(y>=1,1,0)

       i=1
       while i<11:
           a=a-(a+(1/2)*np.sin(2*a)+(np.pi/2)*(1-y*(np.cos(a))**(9/4)))/(1+np.cos(2*a)+(9*np.pi/8)*y*(np.cos(a))**(5/4)*np.sin(a))
           i=i+1
       return a
    # if x<0.25 :
    #       aa=1.07615765*x**(3/4)
    # else :
    #       aa=0.3804792+x**(3/4)*(np.cos(g(x**(-9/8))))**(-3/2)*(3*np.pi/4)**2-3.67345

    # Переделка функции под работу с векторами
    aa = 1.07615765*x**(3/4)*np.where(x<0.25,1,0) + 0.3804792+x**(3/4)*(np.cos(g(x**(-9/8))))**(-3/2)*(3*np.pi/4)**2-3.67345 * np.where(x>=0.25,0,1)
    return aa

#Оригинальная функция
# def mdmt(x) :
#     def g(y):
#        if y<1 :
#          a=0
#        else :
#          a=math.pi/4
#        i=1
#        while i<11:
#            a=a-(a+(1/2)*math.sin(2*a)+(math.pi/2)*(1-y*(math.cos(a))**(9/4)))/(1+math.cos(2*a)+(9*math.pi/8)*y*(math.cos(a))**(5/4)*math.sin(a))
#            i=i+1
#        return a
#     if x<0.25 :
#           aa=1.07615765*x**(3/4)
#     else :
#           aa=0.3804792+x**(3/4)*(math.cos(g(x**(-9/8))))**(-3/2)*(3*math.pi/4)**2-3.67345
#     return aa


#############################################
# далее мои дополнения

def pressure(Rho, Eps):
    gamma = 5/3
    P = (gamma - 1)*Rho*Eps
    return P


