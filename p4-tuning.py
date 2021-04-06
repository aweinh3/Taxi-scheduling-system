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
    K=100
    R=500
    n=200
    alpha=0.5
    f=3
    g=2
    h=1
    M=2000
    T0=2

    #network data
    np.random.seed(903519110)
    D=np.random.randint(3,50,(n,n))
    for i in range(n):
        D[i,i]=0
    T={}
    for i in range(n):
        for j in range(n):
            T[i,j]=D[i,j]

    #trip data
    o=np.random.randint(1,n+1,R)
    d=np.random.randint(1,n+1,R)
    P=np.zeros(R)
    A=np.zeros(R)
    q=np.random.randint(1,4,R)
    t=np.random.randint(0,3,R)
    w=np.random.randint(5,15,R)

    #vehicle data

    u=np.random.randint(1,n+1,K)
    v=np.random.randint(1,n+1,K)
    Q=4

      
    # variable setting
    x={}
    c={}
    B={}
    m={}
    Pu={}
    Au={}
    tu={}
    y={}
    z={}
    for r in range(R):
        for k in range(K+1):
            x[r,k] = gm.addVar(vtype=gp.GRB.BINARY,name='x'+ str(r + 1)+ str(k))
    
    for r in range(R):
        c[r] = gm.addVar(vtype=gp.GRB.BINARY,name='c' + str(r + 1))
        Pu[r] = gm.addVar(vtype=gp.GRB.BINARY,name='Pu'+ str(r + 1))
        Au[r] = gm.addVar(vtype=gp.GRB.INTEGER,lb=0,name='Au'+ str(r + 1))
        tu[r] = gm.addVar(vtype=gp.GRB.INTEGER,lb=0,name='tu'+ str(r + 1))
        for k in range(1,K+1):
            B[r,k] = gm.addVar(vtype=gp.GRB.BINARY,name='B'+ str(r + 1)+ str(k))
            z[r,k] = gm.addVar(vtype=gp.GRB.BINARY,name='B'+ str(r + 1)+ str(k))

    for k in range(1,K+1):
        m[k] = gm.addVar(vtype=gp.GRB.INTEGER,lb=0,name='m'+ str(k))
    for i in range(n):
        for j in range(n):
            for k in range(1,K+1):
                y[i,j,k]= gm.addVar(vtype=gp.GRB.BINARY,name='y'+ str(i + 1)+ str(j + 1) + str(k))

    #objective function        
    gm.setObjective(gp.quicksum(gp.quicksum((f+D[o[r]-1,d[r]-1]*g)*x[r,k] for k in range(1,K+1)) for r in range(R))
                        -gp.quicksum(gp.quicksum(gp.quicksum(y[i,j,k]*D[i,j]*h for j in range(n))for i in range(n))for k in range(1,K+1))
                        -alpha*(gp.quicksum(c[r]*(t[r]-w[r])*(f+D[o[r]-1,d[r]-1]*g) for r in range(R))+gp.quicksum(T[u[k-1]-1,o[r]-1]*z[r,k] for k in range(1,K+1))),gp.GRB.MAXIMIZE)

    #constraints
    #(1)
    for k in range(1,K+1):
        gm.addConstr(gp.quicksum(x[r,k]*(1-P[r]) for r in range(R))<=1)
        gm.addConstr(gp.quicksum(x[r,k] for r in range(R))<=2)
    #(2)
    for r in range(R):
        gm.addConstr(gp.quicksum(x[r,k] for k in range(1,K+1))<=1)

    #(3),(16),(12),(13)
    for r in range(R):
        gm.addConstr(x[r,int(A[r])]*K>=int(A[r]))
        gm.addConstr(t[r]+sum(T[u[k-1]-1,o[r]-1]*x[r,k] for k in range(1,K+1))-w[r]<=M*c[r])    
        gm.addConstr(-(t[r]+sum(T[u[k-1]-1,o[r]-1]*x[r,k] for k in range(1,K+1))-w[r])+M*c[r]<=M)
        gm.addConstr(gp.quicksum(x[r,k]*P[r] for k in range(1,K+1))>=P[r])
        
    #(5),(9)
    for k in range(1,K+1):
         gm.addConstr(m[k]==gp.quicksum(q[r]*P[r]*x[r,k] for r in range(R)))
         gm.addConstr(gp.quicksum(y[u[k-1]-1,i,k] for i in range(n))==1)
    #(6),(7),(8)
    for r in range(R):
        for k in range(1,K+1):
             gm.addConstr(M*(1-(1-P[r])*x[r,k])>=m[k]+q[r]-Q)
             gm.addConstr(1-y[u[k-1]-1,v[k-1]-1,k]<=sum(x[r,k]*(1-P[r]) for r in range(R)))
             gm.addConstr(y[u[k-1]-1,o[r]-1,k]>=(1-P[r])*x[r,k])
    #(10),(11)
    for r in range(R):
        for k in range(1,K+1):
            gm.addConstr(v[k-1]-P[r]*o[r]*x[r,k]+M*B[r,k]>=sum(x[r,k]*(1-P[r]) for r in range(R)))
            gm.addConstr(-v[k-1]+P[r]*o[r]*x[r,k]+M*(1-B[r,k])>=sum(x[r,k]*(1-P[r]) for r in range(R)))
            gm.addConstr(v[k-1]-P[r]*o[r]*x[r,k]+M*B[r,k]>=0)
            gm.addConstr(v[k-1]-P[r]*o[r]*x[r,k]+M*B[r,k]<=M)
            gm.addConstr(z[r,k]<=c[r])
            gm.addConstr(z[r,k]<=x[r,k])
            gm.addConstr(z[r,k]>=c[r]+x[r,k]-1)

    
    #(14),(15)
    for r in range(R):  
        gm.addConstr(Pu[r]==sum(x[r,k] for k in range(1,K+1)))
        gm.addConstr(Au[r]==sum(x[r,k]*k for k in range(1,K+1)))
        gm.addConstr(tu[r]==t[r]+T0*(1-sum(x[r,k] for k in range(1,K+1))))
    
     # Solve the model
    gm.update()
    gm.optimize()
    return(gm)

mg=solve_ip()
print("Optimal Objective Value", mg.objVal)
mg.tune()
count=mg.tuneResultCount
mg.getTuneResult(0)
mg.Params.LogFile="tune"
mg.write('tune'+str(0)+'.prm')

