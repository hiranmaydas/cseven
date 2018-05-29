import math as m
import time
import random


#-----------------------------------------------------------------------
#declaring necessary constants:

TIME_STEPS=10**(-14)
Temperature=90 #K
n=4
N=n**3
M=37.5*0.001 #Kg per mol
N_A=6.023*(10**23) # per mol
m = M/N_A # amu
R=8.314 #Jule/mol.Kelvin
density=2800 #Kg/m^3
box=(n)*((m/density)**(1/3))
rc= (box/n)*3 #cutoff length
sigma=3.40*((10**(-10)))
epsilon = 997/N_A
No_of_loops=100000



#--------------------------------------------------------------------------------------------------------------------------------------------
#initialisation 

def initialisation(particles,Temperature,density):
    i=0
    s_f = (m/density)**(1/3) ## scaling factor  based on density


#placing the atoms in a lattice


    for x in range(n):
        for y in range(n):
            for z in range(n):
                particles[i][0]=[float(x)*s_f,float(y)*s_f,float(z)*s_f]
                i+=1


#giving random velocities to the particles
    rms_v=0.0
    vel_cm=[0.0,0.0,0.0]
    for i in range(i):
        vx=float(random.randint(-100,100))
        vy=float(random.randint(-100,100))
        vz=float(random.randint(-100,100))

        particles[i][1]=[vx,vy,vz]

        vel_cm[0]=vel_cm[0]+vx/N
        vel_cm[1]=vel_cm[1]+vy/N
        vel_cm[2]=vel_cm[2]+vz/N


        rms_v=rms_v+(vx*vx + vy*vy+ vz*vz)/N

    rms_v = rms_v**0.5

    rms_v_actual=((3*R*Temperature)/M)**0.5
    print("rms_v_obtained:",rms_v,"  rms_v_actual:",rms_v_actual,"    vel_of_cm:",vel_cm,"\n")
    correction_factor = rms_v_actual/rms_v


#correcting the velocity in order to match the temperature

    for i in range(N):
        r_x = particles[i][0][0]+(particles[i][1][0]-vel_cm[0])*correction_factor*TIME_STEPS
        r_y = particles[i][0][1]+(particles[i][1][1]-vel_cm[1])*correction_factor*TIME_STEPS
        r_z = particles[i][0][2]+(particles[i][1][2]-vel_cm[2])*correction_factor*TIME_STEPS

        particles[i][1] = [r_x,r_y,r_z] #updating the current position

    return particles



#----------------------------------------------------------------------------------------------------------------------------


def cal_force(particles,force):
    en=0
    for i in range(N):
        force[i]=[0.0,0.0,0.0]

    for i in range(N):
        for j in range(i+1,N):
            xr=particles[i][1][0]-particles[j][1][0] 

            yr=particles[i][1][1]-particles[j][1][1]

            zr=particles[i][1][2]-particles[j][1][2]

            r2= xr**2 + yr**2 + zr**2

            sigma6=sigma**6
            if r2<rc**2:
                r2i=1/r2
                r6i=r2i**3
                ff=r2i*r6i*(r6i*sigma6-0.5)

                force[i][0]=force[i][0] + ff*xr  #updating force imparted on ith particle
                force[i][1]=force[i][1] + ff*yr  #need to revisit if there is problem in simulation. sign may change.
                force[i][2]=force[i][2] + ff*zr


                force[j][0]=force[j][0] - ff*xr  # updating force imparted on jth particle
                force[j][1]=force[j][1] - ff*yr
                force[j][2]=force[j][2] - ff*zr

                en = en + (sigma6*r6i)**2 - (sigma6*r6i) 

    for i in range(N):

        force[i][0] = force[i][0]*48*epsilon*sigma6
        force[i][1] = force[i][1]*48*epsilon*sigma6
        force[i][2] = force[i][2]*48*epsilon*sigma6

    en = en*epsilon*4
    return [force,en]


#==============================================================================================================================================

def integration(particles,force):
    temp=[0,0,0]
    for i in range(N):
        tempx=particles[i][1][0]
        tempy=particles[i][1][1]
        tempz=particles[i][1][2]

        particles[i][1][0]= 2*tempx-particles[i][0][0]+(force[i][0]/m)*(TIME_STEPS**2) #   x
        particles[i][1][1]= 2*tempy-particles[i][0][1]+(force[i][1]/m)*(TIME_STEPS**2) #   y
        particles[i][1][2]= 2*tempz-particles[i][0][2]+(force[i][2]/m)*(TIME_STEPS**2) #   z

        particles[i][0]=[tempx,tempy,tempz]
        

# reflecting boundary condition 
        if particles[i][1][0]>box :                         #   Applying periodic boundary condition 
            particles[i][1][0]=2*box-particles[i][1][0]
            particles[i][0][0]=2*box-particles[i][0][0]


        elif(particles[i][1][0]<0):
            particles[i][1][0]= -particles[i][1][0]
            particles[i][0][0]= -particles[i][0][0]

        if(particles[i][1][1]>box): 
            particles[i][1][1]=2*box -particles[i][1][1] 
            particles[i][0][1]=2*box-particles[i][0][1]

        elif(particles[i][1][1]<0):
            particles[i][1][1]= - particles[i][1][1]
            particles[i][0][1]= -particles[i][0][1]

        if(particles[i][1][2]>box): 
            particles[i][1][2]= 2*box - particles[i][1][2]
            particles[i][0][2]=2*box-particles[i][0][2]

        elif(particles[i][1][2]<0):
            particles[i][1][2]= -particles[i][1][2]
            particles[i][0][2]= -particles[i][0][2]



    return particles   

def cal_vel(particles):
    avg_vel=0
    rms_vel=0

    avg_vx=0
    avg_vy=0
    avg_vz=0
 
    for i in range(N):
        vx=(particles[i][1][0]-particles[i][0][0])/TIME_STEPS
        vy=(particles[i][1][1]-particles[i][0][1])/TIME_STEPS
        vz=(particles[i][1][2]-particles[i][0][2])/TIME_STEPS
        
        v=(vx**2 +vy**2 +vz**2)**(0.5)

        avg_vel=avg_vel + v/N

        avg_vx = avg_vx + vx/N
        avg_vy = avg_vy + vy/N
        avg_vz = avg_vz + vz/N

        rms_vel = rms_vel+ (vx**2 +vy**2 +vz**2)/N
        
    
    rms_vel = rms_vel**0.5
    avg_v = (avg_vx**2 + avg_vy**2 + avg_vz**2)**0.5



    return [avg_vel,rms_vel,avg_vx,avg_vy,avg_vz,avg_v]


def draw_graph(particles,limit):
    
    graph=[]
    for i in range(20):
        graph.insert(i,0)
    
    for j in  range(20):
        ranges=[(limit/20)*(j-1),(limit/20)*(j+1)]
        for i in range(N):
            vx=(particles[i][1][0]-particles[i][0][0])/TIME_STEPS
            vy=(particles[i][1][1]-particles[i][0][1])/TIME_STEPS
            vz=(particles[i][1][2]-particles[i][0][2])/TIME_STEPS
        
            v=(vx**2 +vy**2 +vz**2)**(0.5)
            if(v>ranges[0] and v<=ranges[1]):
                graph[j]+=1
    for i in range(20):
       # for j in range(graph[i]):
       #     print(" ")
      #  print("[*]\n \n \n")
        print(graph[i])
        

#-----------------------------------------------------------------------------------------------------------------------------------------------
#testing the initialisation 
def check_net_force(forces):
    fx=0
    fy=0
    fz=0
    for i in range(N):
        fx = fx + forces[i][0]
        fy = fy + forces[i][1]
        fz = fz + forces[i][2]

    return [fx,fy,fz,forces[0][0]]

def main():

    particles=[]
    for i in range(N):
        particles.insert(i,[[0,0,0],[0,0,0]])


    particles= initialisation(particles,Temperature,density)  #initialisation happenes here
    
    force=[] 

    for i in range(N):
        force.insert(i,[0,0,0])
    energy=0
#----------------------------------------------------------------
    for i in range(No_of_loops):
        force_en = cal_force(particles,force)

        force=force_en[0]

        potential_energy=force_en[1]

        kinetic_energy= 0.5* N * m * cal_vel(particles)[1]**2

        energy = potential_energy + kinetic_energy

        integration(particles,force)

        if(i%1000==0):
            print(i,"\n")
            print("Velocities:",cal_vel(particles),"\n")
            print("Temperature: ",(cal_vel(particles)[1]**2)*(M/(3*R)),"\n")
            print("Total Energy: ",energy*(N_A/N),"\n")
            print("Kinetic Energy: ",kinetic_energy*(N_A/N),"\n")
            print("Potential Energy: ",potential_energy*(N_A/N),"\n")
            print("Net force on system: ", check_net_force(force),"\n")

    print(cal_vel(particles)) 
        

main()
