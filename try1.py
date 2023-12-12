# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import matplotlib.pyplot as plt

# sigma  angstrom
# epsilon  kJ/mol
#Boltzmann constant  kJ/mol⋅K

Ar={'sigma':3.401,'epsilon':0.978638}
He={'sigma':2.556,'epsilon':0.08368}
bolt=8.314462618e-3  

numbsit=[269]
lab=['Ar']
box=np.array([[215.4434],[215.4434],[215.4434]])
rad=1.88
b=8.314462618e-3*273
step=10000
sigma=[3.401]
epsilon=[0.978638]


class model:
    def __init__(self,size,num,rad,b,config,lab,sigma,epsilon):
        self.size=size
        self.num=num
        self.rad=rad
        self.b=b
        self.config=config
        self.tot=0
        for i in num:
            self.tot+=i
        self.lab=lab
        self.perio='off'
        self.sigma=sigma
        self.epsilon=epsilon
        
    def dist(self,i,j,k,z):
        return round(np.linalg.norm(self.config[i].T[j]-self.config[k].T[z]),5)
    
    def initconfig(self):
        self.config=[]
        con=np.random.rand(len(self.size),self.tot)*self.size
        ac=0
        for i in self.num:
            self.config.append(con[:,ac:ac+i])
            ac=ac+i
        while self.checkconfig()==False:
            self.config=[]
            con=np.random.rand(len(self.size),self.tot)*self.size
            ac=0
            for i in self.num:
                self.config.append(con[:,ac:ac+i])
                ac=ac+i
            


    def picture(self,step):
        if len(self.size)==3:
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            plt.title('step : {}'.format(step))
            for i in range(len(self.config)):
                ax.scatter(self.config[i][0],self.config[i][1],self.config[i][2],label=self.lab[i])
            ax.legend()        
            ax.set_xticklabels([])
            ax.set_zticklabels([])
            ax.set_yticklabels([])
            plt.show()
        if len(self.size)==2:
            plt.figure()
            plt.title('step : {}'.format(step))
            for i in range(len(self.config)):
                plt.scatter(self.config[i][0],self.config[i][1],label=self.lab[i])
            plt.legend()
            plt.xlim(0,self.size[0])
            plt.ylim(0,self.size[1])
            
            
    def checkconfig(self):
        chek=0
        for i in range(len(self.config)):
            for j in range(len(self.config[i].T)):
                for k in range(len(self.config)):
                    for z in range(len(self.config[k].T)):
                        if i==k and j==z:
                            pass
                        else:
                            if self.dist(i,j,k,z)<self.rad:
                                chek=chek+1
                            else:
                                pass
        for i in range(len(self.config)):
            for j in range(len(self.config[i].T)):
                if (self.config[i].T[j]>self.size.T-self.rad).any():
                    chek=chek+1
                elif (self.config[i].T[j]<self.rad).any():
                    chek=chek+1
                else:
                    pass
        if chek==0:
            return True
        else:
            return False
        
    def configE(self,index):
        if self.perio=='off':
            if index==1:
                utot=0
                for i in range(len(self.config)):
                    for j in range(len(self.config[i].T)):
                        for k in range(len(self.config)):
                            for z in range(len(self.config[k].T)):
                                if i==k and j==z:
                                    pass
                                else:
                                    utot=utot+4*((epsilon[i]*epsilon[k])**0.5)*((((self.sigma[i]+self.sigma[k])/2)/self.dist(i,j,k,z))**12-(((self.sigma[i]+self.sigma[k])/2)/self.dist(i,j,k,z))**6)
                return utot/2
        else:
            return 0
    #計算config能量
    # 1 : Lennard-Jones
    
    def ranchoiceandmove(self,steplength):
        while True:
            cho0=np.random.choice([1,-1])
            cho1=np.random.randint(0,len(self.config))
            cho2=np.random.randint(0,len(self.config[cho1].T))
            cho3=np.random.randint(0,len(self.size))
            cho4=steplength*float(np.random.rand(),5)
            oldE=self.configE(1)
            self.config[cho1][cho3,cho2]=self.config[cho1][cho3,cho2]+(cho0*cho4*self.size[cho3])/100
            if self.checkconfig():
                if np.exp(-self.b*(self.configE(1)-oldE))>=1:
                    break
                else:
                    if np.random.random()<np.exp(-self.b*(self.configE(1)-oldE)):
                        break
                    else:
                        self.config[cho1][cho3,cho2]=self.config[cho1][cho3,cho2]-(cho0*cho4*self.size[cho3])/100
            else:
                self.config[cho1][cho3,cho2]=self.config[cho1][cho3,cho2]-(cho0*cho4*self.size[cho3])/100
    #隨機取樣與移動
    def rdf(self,lim):
        
        step=0.1
        for i in range(len(self.lab)):
            for k in range(len(self.lab)):
                data=np.zeros(10*lim)
                for z in range(len(self.config[i].T)):
                    dist=[]
                    for j in range(len(self.config[k].T)):
                        dist.append(self.dist(i,z,k,j))
                    for a in range(len(data)):
                        for b in dist:
                            if b>a*step and b<a*step+step:
                                data[a]=data[a]+1
                for c in range(len(data)):
                    data[c]=data[c]/((4*np.pi*(((c+1)*step))**3-(c*step)**3)/3)
                data=data/len(self.config[i].T)
                data[0]=0
                plt.figure()
                plt.plot(np.arange(0,10,0.1),data)
                
                
            
    
mod=model(box,numbsit,rad,b,[],lab,sigma,epsilon)
mod.initconfig()
mod.picture(0)
eig=[]
eig.append(mod.configE(1))
for i in range(step):
    mod.ranchoiceandmove(10)
    eig.append(mod.configE(1))

plt.figure()
plt.plot(range(step+1),eig)
mod.picture(step)