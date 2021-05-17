from AspiR import *
import time


def Test_Wall(filename):
    nx,ny,nd,list_line,robots=ReadingFile(filename)
    X=[]
    Y=[]
    for i in range(len(robots)):
        X.append(int(robots[i][1]))
        Y.append(int(robots[i][2]))
    
    print("The room has ", nx*ny , " boxes, and we want to find  ",12)
    print("\n")
    r=0
    x=X[r]
    y=Y[r]
    W,S,E,N=AuthorizedMove(list_line,x,y,robots,r,X,Y)
    print("W,S,E,N for the first robot : " ,W,S,E,N)
    print("the answer  must be 1 1 0 0 because a robot is on our left")
    print("\n")
    
    r=1
    x=X[r]
    y=Y[r]
    W,S,E,N=AuthorizedMove(list_line,x,y,robots,r,X,Y)
    print("W,S,E,N for the second robot : " ,W,S,E,N)
    print("the answer must be  1 0 0 1 a robot to the north and a wall to the west")
    print("\n")
    
    x=3
    y=2
    W,S,E,N=AuthorizedMove(list_line,x,y,robots,r,X,Y)
    print("W,S,E,N for the boxe 12 ", W,S,E,N)
    print("we want  0 1 1 0, unless you have a robot on an adjacent boxes")
    print("\n")
    
    X[0]=2
    Y[0]=2
    W,S,E,N=AuthorizedMove(list_line,x,y,robots,r,X,Y)
    print("W,S,E,N  for the boxe 12 ", W,S,E,N)
    print("we want 0 1 1 1,because we have a robot in the north of our actual position ")
    print("\n")
    



def Test_Move(filename):
    nx,ny,nd,list_line,robots=ReadingFile(filename)
    X=[]
    Y=[]
    for i in range(len(robots)):
        X.append(int(robots[i][1]))
        Y.append(int(robots[i][2]))
    r=0
    x=X[r]
    y=Y[r]
    print("We consider the robot ", robots[r][0])

    print("he's on  x,y = ",x,y )
    W,S,E,N=AuthorizedMove(list_line,x,y,robots,r,X,Y)
    grid=InitGrid(nx,ny,nd,robots)
    print(grid)
    move_disp=[W,S,E,N] 
    print("my available movements (null values) are: ", move_disp)
    print("\n")
    x,y,grid=Movement(x,y,grid,list_line,2,move_disp) # 0=W 1=S 2=E 3=N
    print("After going east I am in position x,y = ",x,y)
    print("\n")
    print(grid)
    
    print("\n")
    print("\n")
    print("\n")
    

    r=1
    print("I will now consider the robot", robots[r][0])
    print("he's on  x,y = ",x,y )
    W,S,E,N=AuthorizedMove(list_line,x,y,robots,r,X,Y)
    grid=InitGrid(nx,ny,nd,robots)
    move_disp=[W,S,E,N] 
    print("my available movements (null values) are ", move_disp)
    print("\n")
    print("I force the robot to make a forbidden path, if the grid or the position updates it is that there is an error ")
    print("\n")
    print(grid)
    x,y,grid=Movement(x,y,grid,list_line,0,move_disp) # 0=W 1=S 2=E 3=N
    print("After going east I am in position x,y = ",x,y)
    print("\n")
    print(grid)
    


def Test_AspiR():
    
    start_time = time.time()
    max_iter,way=mainprog("Case_Aspi_R_0.txt",2000)
    print("The optimal path ",max_iter, " and it's ", way)
    print("The optimal size's path is  6")
    print(" %s seconds " % (time.time() - start_time))
    print("\n")
    start_time = time.time()
    max_iter,way=mainprog("Case_Aspi_R_1.txt",2000)
    print("The optimal path ",max_iter, " and it's ", way)
    print("The optimal size's path is  6")
    print(" %s seconds " % (time.time() - start_time))
    print("\n")
    start_time = time.time()
    max_iter,way=mainprog("Case_Aspi_R_2.txt", 2000)
    print("The optimal path ",max_iter, "and it's ", way)
    print("The optimal size's path is 10")
    print(" %s seconds " % (time.time() - start_time))
    print("\n")
    
    start_time = time.time()
    max_iter,way=mainprog("Case_Aspi_R_3.txt",300000)
    print("The optimal path ",max_iter, " and it's ", way)
    print("The optimal size's path is 16")
    print(" %s seconds " % (time.time() - start_time))
    print("\n")
    start_time = time.time()
    max_iter,way=mainprog("Case_Aspi_R_4.txt")
    print("The optimal path ",max_iter, " and it's ", way)
    print("The optimal size's path is  12")
    print(" %s seconds " % (time.time() - start_time))
    print("\n")
    
    start_time = time.time()
    max_iter,way=mainprog("Case_Aspi_R_5.txt",50000)
    print("The optimal path ",max_iter, " and it's ", way)
    print("The optimal size's path is 12")
    print(" %s seconds " % (time.time() - start_time))
    print("\n")
    start_time = time.time()
    max_iter,way=mainprog("Case_Aspi_R_6.txt")
    print("The optimal path",max_iter, " and it's ", way)
    print("The optimal size's path is 14")
    print(" %s seconds " % (time.time() - start_time))
    print("\n")
    
    start_time = time.time()
    max_iter,way=mainprog("Case_Aspi_R_7.txt",100000)
    print("The optimal path ",max_iter, " and it's ", way)
    print("The optimal size's path is 14")
    print(" %s seconds " % (time.time() - start_time))
    

#Test_Wall("Case_Aspi_R_0.txt")
#Test_Move("Case_Aspi_R_0.txt")
Test_AspiR()

    
    
    

     
