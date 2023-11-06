# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import matplotlib.pyplot as plt

numbsit=[5,5]
lab=['Ar','He']
box=np.array([[100],[100]])
rad=0.5
b=1
step=20000
sigma=[10,4]
epsilon=[4,2]



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
        return np.linalg.norm(self.config[i].T[j]-self.config[k].T[z])
    
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
            oldE=self.configE(1)
            self.config[cho1][cho3,cho2]=self.config[cho1][cho3,cho2]+(cho0*steplength*self.size[cho3])/100
            if self.checkconfig():
                if np.exp(-self.b*(self.configE(1)-oldE))>=1:
                    break
                else:
                    if np.random.random()<np.exp(-self.b*(self.configE(1)-oldE)):
                        break
                    else:
                        self.config[cho1][cho3,cho2]=self.config[cho1][cho3,cho2]-(cho0*steplength*self.size[cho3])/100
            else:
                self.config[cho1][cho3,cho2]=self.config[cho1][cho3,cho2]-(cho0*steplength*self.size[cho3])/100
    #隨機取樣與移動
    
    
mod=model(box,numbsit,rad,b,[],lab,sigma,epsilon)
mod.initconfig()
mod.picture(0)
eig=[]
eig.append(mod.configE(1))
for i in range(step):
    mod.ranchoiceandmove(3)
    eig.append(mod.configE(1))

plt.figure()
plt.plot(range(step+1-50),eig[50:])
mod.picture(step)