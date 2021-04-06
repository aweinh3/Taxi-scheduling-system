# Implement the Extreme test cases (one vehicle with 100 requests)
import gurobipy as gp
import numpy as np
from gurobipy import *
import matplotlib.pyplot as plt
 

# Invoke Gurobi to solve the QP
def solve_ip():

    # Model setup
    gm = gp.Model("test_ip")

    # Create the data
    K=30
    R=90
    n=2*R+2
    alpha=0.5
    f=3
    g=2
    h=1
    M=2000


    #network data
    np.random.seed(903519110)
    D=np.random.randint(3,50,(n,n))
    for i in range(n):
        D[i,i]=0
        for j in range(n):
            D[i,j]=D[j,i]
    T={}
    for i in range(n):
        for j in range(n):
            T[i,j]=D[i,j]

    #data input
    q={}   
    #vertices loads
    q[0]=q[2*R+1]=0
    for i in range(1,R+1):
        q[i]=np.random.randint(1,4)
        q[i+R]=-q[i]

        
    #maximum waittime and assignment status of each request
    m={}
    P={}
    for i in range(1,R+1):
        m[i]=np.random.randint(5,15)
        P[i]=0
       
    
    #taxi capacity
    Q={}
    for k in range(K):
        Q[k]=4

      
    # variable setting
    x={}
    c={}
    t={}
    w={}
    
    for i in range(n):
        for j in range(n):
            for k in range(K):
                x[i,j,k]= gm.addVar(vtype=gp.GRB.BINARY,name='x'+ str(i + 1)+ str(j + 1) + str(k))
     
    for i in range(1,1+R):
        c[i]= gm.addVar(vtype=gp.GRB.BINARY,name='c'+ str(i))
         
    for i in range(n):
        for k in range(K):
            t[i,k] = gm.addVar(name='t'+ str(i + 1)+ str(k))
            w[i,k] = gm.addVar(vtype=gp.GRB.INTEGER,name='w'+ str(i + 1)+ str(k))
      
        
    #objective function        
    gm.setObjective(gp.quicksum(gp.quicksum(gp.quicksum((f+D[i,i+R]*g)*x[i,j,k] for k in range(K)) for j in range(n)) for i in range(1,1+R))
                        -gp.quicksum(gp.quicksum(gp.quicksum(x[i,j,k]*D[i,j]*h for j in range(n))for i in range(n))for k in range(K))
                        -alpha*(gp.quicksum(c[i]*(gp.quicksum(t[i,k] for k in range(K))-m[i])*(f+D[i,i+R]*g) for i in range(1,1+R))),gp.GRB.MAXIMIZE)
           
    #constraints
    #(2)
    for i in range(1,1+R):
        gm.addConstr(gp.quicksum(gp.quicksum(x[i,j,k] for k in range(K)) for j in range(n))<=1)
        for k in range(K):
             gm.addConstr(gp.quicksum(x[j,i,k] for j in range(n))==gp.quicksum(x[j,i+R,k] for j in range(n)))
    
    #(1)
    for k in range(K):
        gm.addConstr(gp.quicksum(x[0,i,k]*(1-P[i]) for i in range(1,R+1))<=1)

   #(3) #problematic
    for k in range(K):
        gm.addConstr(gp.quicksum(x[0,i,k] for i in range(n))==1)
        gm.addConstr(gp.quicksum(x[i,1+2*R,k] for i in range(n))==1)
        gm.addConstr(gp.quicksum(x[i,0,k] for i in range(n))==0)
        gm.addConstr(gp.quicksum(x[1+2*R,i,k] for i in range(n))==0)
    
    #(4)
    for i in range(1,2*R+1):
        for k in range(K):
            gm.addConstr(gp.quicksum(x[i,j,k] for j in range(n))==gp.quicksum(x[j,i,k] for j in range(n)))
               
    #(5),(6)
    for i in range(n):
        for j in range(n):
            for k in range(K):
                gm.addConstr(t[j,k]>=t[i,k]+T[i,j]+M*(x[i,j,k]-1))
                gm.addConstr(w[j,k]>=w[i,k]+q[i]+M*(x[i,j,k]-1))
    
    for k in range(K):
        gm.addConstr(t[0,k]==0)
        gm.addConstr(w[0,k]==0)
        
    #(8)           
    for i in range(n):
        for k in range(K):
             gm.addConstr(w[i,k]>=0)
             gm.addConstr(w[i,k]>=q[i])
             gm.addConstr(w[i,k]<=Q[k])
             gm.addConstr(w[i,k]<=Q[k]+q[i])
             
    #(9)
    for i in range(1,1+R):  
        gm.addConstr(-(gp.quicksum(t[i,k] for k in range(K))-m[i])+M*c[i]<=M)
        gm.addConstr(-(gp.quicksum(t[i,k] for k in range(K))-m[i])+M*c[i]>=0)
        

    # Solve the model
    gm.update()
    gm.Params.Method = 2.0
    gm.Params.MIPFocus = 1.0
    gm.optimize()
    return(gm)

mg=solve_ip()
print("Optimal Objective Value", mg.objVal)



