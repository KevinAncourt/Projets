#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 13:32:14 2020

@author: kevinancourt
"""
import numpy as np
from random import *
from math import *
import time

def ReadingFile(input_filename):    
    file = open(input_filename, "r")   #Opens the file 
    line_0=file.readline() 
    line_0=line_0.split(" ")  
    nx=int(line_0[0]) # Number of lines 
    ny=int(line_0[1]) # Number of columns
    nd=int(line_0[2]) # Number of  vacuum cleaners
    list_line=[]
    robots=[]
    for i in range (nx): 
        line_i=file.readline()
        line_i=line_i.rstrip('\n')
        line=list(line_i)
        list_line.append(line)    
    for i in range(nd):
        line_i=file.readline()
        line_i=line_i.rstrip('\n')
        line_i=line_i.split(" ")# splits the line 
        robot=list(line_i)
        robots.append(robot)
    file.close()    
    return nx,ny,nd,list_line,robots
    
def InitGrid(nx,ny,nd,robots):
    grid=np.zeros((nx,ny)) #Creates a grid , 1 if we already visited this place , otherwise 0
    for i in range(0,nd):
        x=int(robots[i][1])
        y=int(robots[i][2])
        grid[x][y]=1 #Initializes the position of the vacuum cleaners
    return grid

def HexeToBin(n):
    return format(int(n,16),'#010b')[6:] #Converts to binary
   
def AuthorizedMove(list_line,x,y,robots,r,X,Y):
    robot_line=x 
    robot_column=y
    box_described=list_line[robot_line][robot_column] # the position of the walls around the vacuum cleaner
    move=HexeToBin(box_described);
    W=0
    S=0
    E=0
    N=0
    for i in range(len(robots)):
        if i != r:
            if (move[0]=='1') or ((Y[i]==y-1) and X[i]==x):
                W=1
            if (move[1]=='1') or ((X[i]==x+1) and Y[i]==y):
                S=1
            if (move[2]=='1') or ((Y[i]==y+1) and X[i]==x):
                E=1
            if (move[3]=='1') or ((X[i]==x-1) and Y[i]==y):
                N=1   
    return W,S,E,N        

def Direction_Chosen(W,S,E,N):
    move_disp=[W,S,E,N]  
    move_potential=[]
    for i in range(4):
        if move_disp[i]==0:
            move_potential.append(i)
    if move_potential!=[]:
        direction=choice(move_potential) # Chooses a direction randomly among those available
    else :
        direction=1 # if we are blocked, we take the first direction     
    return direction,move_disp
    

def Movement(x,y,grid,list_line,direction,move_disp):
    list_direction=['W','S','E','N']
    if direction==0 and y-1 >= 0 and move_disp[direction]==0:  #Checks if we can move on the direction chosen
        y=y-1
        grid[x][y]=1 #adds 1 to the case clean on the grid 
    if direction==1  and x+1<=len(list_line)-1 and move_disp[direction]==0:
        x=x+1
        grid[x][y]=1
    if direction==2  and y+1<=len(list_line[0])-1 and move_disp[direction]==0:
        y=y+1
        grid[x][y]=1
    if direction==3  and x-1>=0 and move_disp[direction]==0:
        x=x-1
        grid[x][y]=1   
    return x,y,grid  

def Movement_Loop(robots,r,x,y,grid,list_line,test,X,Y): # Loop on the Movement until the cleaner cross a wall
	W,S,E,N=AuthorizedMove(list_line,x,y,robots,r,X,Y)
	espglout = randint(1,10)
	if(espglout<7):
		W,S,E,N = Glout(W,S,E,N,r,grid,robots,X,Y,list_line)
	direction,move_disp= Direction_Chosen(W,S,E,N)
	x,y,grid=Movement(x,y,grid,list_line,direction,move_disp)
	list_direction=['W','S','E','N']
	test.append(robots[r][0]+list_direction[direction])
	while move_disp[direction] !=1:
		W,S,E,N=AuthorizedMove(list_line,x,y,robots,r,X,Y)
		move_disp=[W,S,E,N]
		x,y,grid=Movement(x,y,grid,list_line,direction,move_disp)
	return x,y,grid,test

def Glout(W,S,E,N,r,grid,robots,X,Y,list_line): #privileges the dirtiest direction 
    compt_W=0
    compt_S=0
    compt_E=0
    compt_N=0
    if W==0:
        x_prime=X[r]
        y_prime=Y[r]
        W_prime=W
        S_prime=S
        E_prime=E
        N_prime=N
        while W_prime!=1:
            y_prime=y_prime-1
            W_prime,S_prime,E_prime,N_prime=AuthorizedMove(list_line,x_prime,y_prime,robots,r,X,Y)
            if grid[x_prime][y_prime]==0:
                compt_W = compt_W +1
    if E==0:
        x_prime=X[r]
        y_prime=Y[r]
        W_prime=W
        S_prime=S
        E_prime=E
        N_prime=N
        while E_prime!=1:
            y_prime=y_prime+1
            W_prime,S_prime,E_prime,N_prime=AuthorizedMove(list_line,x_prime,y_prime,robots,r,X,Y)
            if grid[x_prime][y_prime]==0:
                compt_E = compt_E +1
    if S==0:
        x_prime=X[r]
        y_prime=Y[r]
        W_prime=W
        S_prime=S
        E_prime=E
        N_prime=N
        while S_prime!=1:
            x_prime=x_prime+1
            W_prime,S_prime,E_prime,N_prime=AuthorizedMove(list_line,x_prime,y_prime,robots,r,X,Y)
            if grid[x_prime][y_prime]==0:
                compt_S = compt_S +1    
    if N==0:
        x_prime=X[r]
        y_prime=Y[r]
        W_prime=W
        S_prime=S
        E_prime=E
        N_prime=N
        while N_prime!=1:
            x_prime=x_prime-1
            W_prime,S_prime,E_prime,N_prime=AuthorizedMove(list_line,x_prime,y_prime,robots,r,X,Y)
            if grid[x_prime][y_prime]==0:
                compt_N = compt_N +1
    if compt_W<max(compt_S,compt_N,compt_E): W=1
    if compt_S<max(compt_W,compt_N,compt_E): S=1
    if compt_E<max(compt_S,compt_N,compt_W): E=1
    if compt_N<max(compt_S,compt_W,compt_E): N=1
    
    return W,S,E,N
                
def CheckGrid(grid,nx,ny): # Checks if the all grid is clean 
    End=0
    for i in range(nx):
        for j in range(ny):
            if grid[i][j]==0:
                End=End+1
    return End
     
def Try(max_iter,robots,list_line,grid,nx,ny):   #Tries a lot of posibilities      
    iter=0
    End=1
    X=[]
    Y=[]
    for i in range(len(robots)):
        X.append(int(robots[i][1]))
        Y.append(int(robots[i][2]))
    test=[]
    while iter<max_iter and End!=0:
        r= randint(0,len(robots)-1)
        x=X[r]
        y=Y[r]
        x,y,grid,test=Movement_Loop(robots,r,x,y,grid,list_line,test,X,Y)
        End=CheckGrid(grid,nx,ny)
        X[r]=x
        Y[r]=y
        iter=iter+1
    it=len(test) #lenght of the deplacements
    return test,it,grid

def mainprog(filename,i_max=1000000):
    nx,ny,nd,list_line,robots=ReadingFile(filename)
    grid=InitGrid(nx,ny,nd,robots)
    way=[]
    gridf=[]
    max_iter=nx*ny
    i=0
    while i< i_max:
        grid=InitGrid(nx,ny,nd,robots)
        test,it,grid=Try(max_iter,robots,list_line,grid,nx,ny)
        i=i+1
        if it<max_iter:
            max_iter=it #if we find a good answer we don't look for a slower one
            way=test
            gridf=grid
    #print(max_iter)
    #print(way)
    #print(gridf)

    return max_iter,way

"""
start_time = time.time()
max_iter,way=mainprog("Case_Aspi_R_4.txt")
print(" %s secondes " % (time.time() - start_time))
print(max_iter)
print(way)
"""
