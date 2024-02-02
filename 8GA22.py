import copy
import random as ra
import math as ma
import matplotlib.pyplot as plt
## Initializes constants and variables  
crossover_rate = .9
mutation_rate = .04
accuracy_rate = .75
gen_count = 10
best_fitness = 0
best_chromosome = []
pop1 = []
popcount = 4
grid = {}
length = 5
breadth = 5
inputnum = 8
hlayerend = [0,0,0,0]
hlayer1  = [0 for i in range(length*breadth*2)]
database = []
for i in range(popcount+1):
    database.append([])
def sig(x):
    return 1/(1+ma.exp(-x))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
indiv_coord=[0,0]

chrom_input = {}
eneminput = {}

def limit_func(x,limit):
    if x+1 <= limit:
        return x+1
    else:
        return x -1


## Creates the blueprint for the whole population
def chrom_initiation():
    indiv=[[ra.uniform(-4,4) for i in range(length*breadth*inputnum)] for j in hlayer1], [[ra.uniform(-4,4) for y in hlayer1] for z in hlayerend]

    return list(indiv)
## Creates the population of individuals
def pop_initiation(pop_size):
    global pop1
    pop1 = [chrom_initiation() for x in range (1,pop_size) ]
    return pop1

## Creates the input map for feach individual's choice algorithm
for i in range(length+1):
    chrom_input[i]={}
    for j in range(breadth+1):
        chrom_input[i][j]=[0 for i in range(inputnum)]

for i in range(length+1):
    eneminput[i]={}
    for j in range(breadth+1):
        eneminput[i][j]=[0 for i in range(inputnum)]

## Creates grid where battle takes place
for i in range(length+1):
    grid[i] = {}
    for j in range(breadth+1):
        grid[i][j]=[1,0,0,0,[0,0]]

## Cross references points on grid to map
def crossreference(grid,chrom_input,control,enemcontrol):    
    for i in range(length+1):
        for j in range(breadth+1):
            for x in range(0,4):
                if grid[i][j][x] == 1:
                    chrom_input[i][j][x] =1
                else:
                    chrom_input[i][j][x] =0
            for y in ((1,0,1),(-1,0,2),(0,1,3),(0,-1,4)): 
                try:
                    if grid[abs(control[0]+y[0])][abs(control[1]+y[1])]==[]:
                        chrom_input[i][j][y[2]+2] = 0
                    else:
                        if grid[abs(control[0]+y[0])][abs(control[1]+y[1])][-1][0] >= grid[control[0]][control[1]][-1][0]:
                            chrom_input[i][j][y[2]+2] = 1
                        else:
                            chrom_input[i][j][y[2]+2] = 0
                except KeyError:
                    if grid[abs(control[0]-y[0])][abs(control[1]-y[1])]==[]:
                        chrom_input[i][j][y[2]+2] = 0
                    else:
                        if grid[abs(control[0]-y[0])][abs(control[1]-y[1])][-1][0] >= grid[control[0]][control[1]][-1][0]:
                            chrom_input[i][j][y[2]+2] = 1
                        else:
                            chrom_input[i][j][y[2]+2] = 0
    return chrom_input

def enem_crossreference(grid,eneminput,control,enemcontrol):
    grid[control[0]][control[1]][3], grid[enemcontrol[0]][enemcontrol[1]][4], grid[control[0]][control[1]][4], grid[enemcontrol[0]][enemcontrol[1]][3]= 0,0,1,1
    for i in range(length+1):
        for j in range(breadth+1):
            if grid[i][j][0] == 1:
                eneminput[i][j][0] = 1
            else:
                eneminput[i][j][0] = 0
            for f in range(1,2):
                if grid[i][j][f]==1:
                    eneminput[i][j][f]=0
                else:
                    eneminput[i][j][f]=1
            for x in range(3,4):
                if grid[i][j][x] == 1:
                    eneminput[i][j][x] =1
                elif grid[i][j][x] ==0:
                    eneminput[i][j][x] =0
            for y in ((1,0,1),(-1,0,2),(0,1,3),(0,-1,4)): 
                try:
                    if grid[abs(control[0]+y[0])][abs(control[1]+y[1])]==[]:
                        eneminput[i][j][y[2]+2] = 0
                    else:
                        if grid[abs(control[0]+y[0])][abs(control[1]+y[1])][-1][0] >= grid[control[0]][control[1]][-1][0]:
                            eneminput[i][j][y[2]+2] = 1
                        else:
                            eneminput[i][j][y[2]+2] = 0
                except KeyError:
                    if grid[abs(control[0]-y[0])][abs(control[1]-y[1])]==[]:
                        eneminput[i][j][y[2]+2] = 0
                    else:
                        if grid[abs(control[0]-y[0])][abs(control[1]-y[1])][-1][0] >= grid[control[0]][control[1]][-1][0]:
                            eneminput[i][j][y[2]+2] = 1
                        else:
                            eneminput[i][j][y[2]+2] = 0
    return eneminput


## Creates new members of the population based off of the previous gen. Either through duplication with some mutation, or crossover
def new_pop(pop1):
    new_pop = []
    while len(new_pop) != popcount:
        if len(new_pop) > popcount:
            while len(new_pop) > popcount:
                new_pop.remove(ra.choice(new_pop))
            break
        if len(new_pop) < popcount:
            pass
        if crossover_rate >= ra.uniform(0,1):
            new_chromo = crossover(pop1)
            new_pop.append(new_chromo[0])
            new_pop.append(new_chromo[1])
        else:
            new_chromo = mutation(roulette_selection(pop1)) 
            new_pop.append(new_chromo)
    return new_pop

## Selects a random individual in a weighted roulette wheel selection
def roulette_selection(pop1,chrome_input,eniminput):
    wheel = []
    for x in pop1:
        if database.index[pop1.index(x)] != []:
            fitnessvalue = database.index[pop1.index(x)]
        else:    
            fitnessvalue = war(x,chrom_input,eneminput,ra.choice(pop1),pop1)
            database[pop.index(x)] = fitnessvalue
        for i in range(fitnessvalue):
            wheel.append(x)
    choice = ra.choice(wheel)
    return choice

## Crossovers the individuals. Is the basis for the reproduction system
def crossover(pop1):
    a = roulette_selection(pop1)
    b = roulette_selection(pop1)
    p1 = ra.randint(1,len(hlayer1))
    try:
        p2 = ra.randint(1,len(hlayer1)-p1)
    except:
        p2=1
    individuals = two_point_crossover(a,b,p1,p2)
    c = individuals[0]
    d = individuals[1]
    return mutation(c), mutation(d)

## The details of how the crossover functions. Two point crossover.
def two_point_crossover(a, b, p1, p2):
    c = a[:p1] + b[p1:p2] + a[p2:] 
    d = b[:p1] + a[p1:p2] + b[p2:]
    for i in hlayer1:
        p3 = ra.randint(1,length*breadth*inputnum)
        num = length*breadth*inputnum-p3
        
        try:
            p4 = ra.randint(1,num)
        except:
            p4=1
        c[0][hlayer1.index(i)] = a[0][hlayer1.index(i)][:p3] + b[0][hlayer1.index(i)][p3:p4] + a[0][hlayer1.index(i)][p4:] 
        d[0][hlayer1.index(i)] = b[0][hlayer1.index(i)][:p3] + a[0][hlayer1.index(i)][p3:p4] + b[0][hlayer1.index(i)][p4:] 
    
    for i in hlayerend:
        p3 = ra.randint(1,len(hlayer1))
        try:
            p4 = ra.randint(1,len(hlayer1)-p3)
        except:
            p4=1
 
        c[1][hlayerend.index(i)] = a[1][hlayerend.index(i)][:p3] + b[1][hlayerend.index(i)][p3:p4] + a[1][hlayerend.index(i)][p4:] 
        d[1][hlayerend.index(i)] = b[1][hlayerend.index(i)][:p3] + a[1][hlayerend.index(i)][p3:p4] + b[1][hlayerend.index(i)][p4:] 
    return c, d

## Mutates the individual.
def mutation(indiv):
    if mutation_rate >= ra.uniform(0,1):
            choice1 = ra.choice(indiv)
            choice2 = ra.choice(indiv)
            a, b = indiv.index(choice1), indiv.index(choice2)
            indiv[b], indiv[a] = indiv[a], indiv[b]
    for i in indiv:
        for j in indiv[0][indiv.index(i)]:
            if mutation_rate >= ra.uniform(0,1):
                choice1 = ra.choice(indiv[0][indiv.index(i)])
                choice2 = ra.choice(indiv[0][indiv.index(i)])
                a, b = indiv[0][indiv.index(i)].index(choice1), indiv[0][indiv.index(i)].index(choice2)
                indiv[0][indiv.index(i)][b], indiv[0][indiv.index(i)][a] = indiv[0][indiv.index(i)][a], indiv[0][indiv.index(i)][b]
        for j in indiv[1][indiv.index(i)]:
            if mutation_rate >= ra.uniform(0,1):
                choice1 = ra.choice(indiv[1][indiv.index(i)])
                choice2 = ra.choice(indiv[1][indiv.index(i)])
                a, b = indiv[1][indiv.index(i)].index(choice1), indiv[1][indiv.index(i)].index(choice2)
                indiv[1][indiv.index(i)][b], indiv[1][indiv.index(i)][a] = indiv[1][indiv.index(i)][a], indiv[1][indiv.index(i)][b]
    return indiv

def lanchester_battle(indiv,enem):
    enem_pop = enem[-1][0]/2
    indiv_pop = indiv[-1][0]/2
    
    while indiv_pop and enem_pop >= 1:
        
        pop1 = enem_pop
        pop2 = indiv_pop
        indiv_pop -= (pop1*(enem[-1][1]*accuracy_rate))
        enem_pop -= (pop2*(indiv[-1][1]*accuracy_rate))
 
        if enem_pop or indiv_pop <= 0:
            if enem_pop >= indiv_pop:
                winner = ('enem',enem_pop)
            else:
                winner = ('indiv',indiv_pop)

            break
    return winner
                
def war(indiv,chrom_input,eneminput,enem,pop1):
    global database
    global grid
    gride = grid
    gride[0][0]=[0,0,1,1,0,[1500,.5]]
    gride[length][breadth]=[0,1,0,0,1,[1500,.5]]
    enemtemplate = 0,1,0,0,1
    indivtemplate = 0,0,1,1,0
    for i in range(length+1):
        for j in range(breadth+1):
                if gride[i][j][3] == 1:
                    troops_in_control = [i,j]
                if gride[i][j][4]==1:
                    enemcontrol = [i,j]
    for i in range(length+1):
        chrom_input = crossreference(gride,chrom_input,troops_in_control,enemcontrol)
        eneminput = enem_crossreference(gride,eneminput,troops_in_control,enemcontrol)
        dec = decision(chrom_input,indiv,gride)
        output = max(dec)
        grid = turn (gride,indiv,troops_in_control,dec,output,indivtemplate)
        dec2 = decision(eneminput,enem,gride)
        output2 = max(dec2)
        grid=turn(gride,enem,enemcontrol,dec2,output2,enemtemplate)
    squarecount = 0 
    enemcount = 0
    for i in range(length+1):
        for j in range(breadth+1):
            if gride[i][j][2] == 1:
                squarecount += 1
            if gride[i][j][1] == 1:
                enemcount += 1
    database[pop1.index(enem)] = enemcount
    return squarecount

def turn(grid,indiv,troops_in_control,dec,output,template):
    if dec.index(output) == 0:
        if grid[limit_func(troops_in_control[0],length)][troops_in_control[1]][1] == 1:
            simulation = lanchester_battle(grid[troops_in_control[0]][troops_in_control[1]],grid[limit_func(troops_in_control[0],length)][troops_in_control[1]])
            if simulation[0] == 'enem':
                grid[troops_in_control[0]][troops_in_control[1]][-1][0] = grid[troops_in_control[0]][troops_in_control[1]][-1][0]/2
                grid[limit_func(troops_in_control[0],length)][troops_in_control[1]][-1][0] = abs(simulation[1]-(grid[limit_func(troops_in_control[0],length)][troops_in_control[1]])[-1][0]/2)
                return grid
            else:
                grid[limit_func(troops_in_control[0],length)][troops_in_control[1]][-1][0] = grid[troops_in_control[0]][troops_in_control[1]][-1][0]/2
                grid[troops_in_control[0]][troops_in_control[1]][-1][0] = abs(simulation[1]-(grid[troops_in_control[0]][troops_in_control[-1]])[-1][0]/2)
                return grid
        grid[limit_func(troops_in_control[0], length)][troops_in_control[1]] = [[template[x-1] for x in range(1, len(template))], [grid[limit_func(troops_in_control[0], length)][troops_in_control[1]][-1][0] + grid[troops_in_control[0]][troops_in_control[1]][-1][0] / 2, .5]]

        half = grid[troops_in_control[0]][troops_in_control[1]][-1][0]/2
        grid[troops_in_control[0]][troops_in_control[1]] = [[template[x-1] for x in range(1,len(template))[half,.5]]]
        return grid
    if dec.index(output) == 1:
        if grid[abs(troops_in_control[0]-1)][troops_in_control[1]][1] == 1:
            simulation = lanchester_battle(grid[troops_in_control[0]][troops_in_control[1]],grid[abs(troops_in_control[0]-1)][troops_in_control[1]])
            if simulation[0] == 'enem':
                grid[troops_in_control[0]][troops_in_control[1]][-1][0] = grid[troops_in_control[0]][troops_in_control[1]][-1][0]/2
                grid[abs(troops_in_control[0]-1)][troops_in_control[1]][-1][0] = abs(simulation[1]-(grid[abs(troops_in_control[0]-1)][troops_in_control[1]])[-1][0]/2)
                return grid
            else:
                grid[abs(troops_in_control[0]-1)][troops_in_control[1]][-1][0] = grid[troops_in_control[0]][troops_in_control[1]][-1][0]/2
                grid[troops_in_control[0]][troops_in_control[1]][-1][0] = abs(simulation[1]-(grid[troops_in_control[0]][troops_in_control[1]])[-1][0]/2)
                return grid
        grid[abs(troops_in_control[0]-1)][troops_in_control[1]] = [template[x-1] for x in range(1,len(template))[grid[abs(troops_in_control[0]-1)][troops_in_control[1]][-1][0]+grid[troops_in_control[0]][troops_in_control[1]][-1][0]/2,.5]]
        half = grid[troops_in_control[0]][troops_in_control[1]][-1][0]/2
        grid[troops_in_control[0]][troops_in_control[1]] = [template[x-1] for x in range(1,len(template))[half,.5]]
        return grid
    if dec.index(output) == 2:
        if grid[troops_in_control[0]][limit_func(troops_in_control[1],breadth)][1] == 1:
            simulation = lanchester_battle(grid[troops_in_control[0]][troops_in_control[1]],grid[troops_in_control[0]][limit_func(troops_in_control[1],breadth)])
            if simulation[0] == 'enem':
                grid[troops_in_control[0]][troops_in_control[1]][-1][0] = grid[troops_in_control[0]][troops_in_control[1]][-1][0]/2
                grid[troops_in_control[0]][limit_func(troops_in_control[1],breadth)][-1][0] = abs(simulation[1]-(grid[troops_in_control[0]][limit_func(troops_in_control[1],breadth)])[-1][0]/2)
                return grid
            else:
                grid[troops_in_control[0]][limit_func(troops_in_control[1],breadth)][-1][0] = grid[troops_in_control[0]][troops_in_control[1]][-1][0]/2
                grid[troops_in_control[0]][troops_in_control[1]][-1][0] = abs(simulation[1]-(grid[troops_in_control[0]][troops_in_control[1]])[-1][0]/2)
                return grid
        grid[troops_in_control[0]][limit_func(troops_in_control[1],breadth)] = [template[x-1] for x in range(1,len(template))[grid[troops_in_control[0]][limit_func(troops_in_control[1],breadth)][-1][0]+grid[troops_in_control[0]][troops_in_control[1]][-1][0]/2,.5]]
        half = grid[troops_in_control[0]][troops_in_control[1]][-1][0]/2
        grid[troops_in_control[0]][troops_in_control[1]] = [template[x-1] for x in range(1,len(template))[half,.5]]
        return grid
    if dec.index(output) == 3:
        if grid[troops_in_control[0]][abs(troops_in_control[1]-1)][1] == 1:
            simulation = lanchester_battle(grid[troops_in_control[0]][troops_in_control[1]],grid[troops_in_control[0]][abs(troops_in_control[1]-1)])
            if simulation[0] == 'enem':
                grid[troops_in_control[0]][troops_in_control[1]][-1][0] = grid[troops_in_control[0]][troops_in_control[1]][-1][0]/2
                grid[troops_in_control[0]][abs(troops_in_control[1]-1)][-1][0] = abs(simulation[1]-(grid[troops_in_control[0]][abs(troops_in_control[1]-1)])[-1][0]/2)
                return grid
            else:
                grid[troops_in_control[0]][abs(troops_in_control[1]-1)][-1][0] = grid[troops_in_control[0]][troops_in_control[1]][-1][0]/2
                grid[troops_in_control[0]][troops_in_control[1]][-1][0] = abs(simulation[1]-(grid[troops_in_control[0]][troops_in_control[1]])[-1][0]/2)
                return grid
        grid[troops_in_control[0]][abs(troops_in_control[1]-1)] = [template[x-1] for x in range(1,len(template))[grid[troops_in_control[0]][abs(troops_in_control[1]-1)][-1][0]+grid[troops_in_control[0]][troops_in_control[1]][-1][0]/2,.5]]
        half = grid[troops_in_control[0]][troops_in_control[1]][-1][0]/2
        grid[troops_in_control[0]][troops_in_control[1]] = [template[x-1] for x in range(1,len(template))[half,.5]]    
        return grid

def decision(chrom_input,indiv,grid):
    for i in (chrom_input[x][y][j] for x in range(length+1) for y in range(breadth+1)for j in range(inputnum)):
        for j in hlayer1:
            count=0
            k=hlayer1.index(j)
            # if k < len(indiv[0]):
            #             print('true1')
            # else:
            #     print('false')
            try:
                hlayer1[hlayer1.index(j)]=sig(j+(i*indiv[0][count][list(chrom_input).index(i)]))
            except Exception as e:
                count+=1
                if count+1> len(indiv[0]):
                    count=0
                print(indiv)
                print(len(indiv[0]))
                print(count)
                print(j,i)
                print(hlayer1.index(j),list(chrom_input).index(i))

                if indiv[0][hlayer1.index(j)-1] and indiv[0][hlayer1.index(j)+1] != None:
                    print('true') 
                print(e)
                quit()
    for i in hlayer1:
        for j in hlayerend:
            hlayerend[hlayerend.index(j)]=sig(j+(i*indiv[1][hlayerend.index(j)][list(hlayer1).index(i)]))        
    outputlayer = hlayerend
    return outputlayer

def besfitness(pop1,popfitnesss):
    indices = popfitnesss.index(max(popfitnesss))
    return[pop1[indices],max(popfitnesss)]

def popfitnesses(pop1,chrome_input,eniminput):
 
    popfitness = []
    for indiv in pop1:
        if database[pop1.index(indiv)] != []:
            popfitness.append(database.index[pop1.index(indiv)])
        else:    
            fitnessvalue = war(indiv,chrom_input,eneminput,ra.choice(pop1),pop1)
            database[pop1.index(indiv)] = fitnessvalue 
            popfitness.append(war(indiv,chrom_input,eneminput,ra.choice(pop1),pop1))
    return popfitness

# dele = open("GAgraphinfo.txt","w")
# dele.close()
# text = open('GAgraphinfo.txt', 'a')
pop1 = pop_initiation(popcount+1)
y = []
y2 = []
x = []
gen = 0

for i in range(1,gen_count+1):
    database = []
    for i in range(popcount+1):
        database.append([])
    popfitness = [0]
    print('go')
    gen += 1
    popfitness = popfitnesses(pop1,chrom_input,eneminput)
    print ("Gen: "+str(gen))
    # print ("Best Solution: "+str(besfitness(pop,popfitness)[0]))
    print ("Best Fitness: "+str(besfitness(pop1,popfitness)[1]))
    print(len(pop1))
    for i in pop1:
        if type(i) == 'Tuple':
            print('yes')
    y.append(gen)
    x.append(besfitness(pop1,popfitness)[1])
    pop_copy = copy.deepcopy(pop1)
    newpop = new_pop(pop_copy)
    pop = newpop
 
print('end')

plt.plot(y,x)
plt.title('Best Fitness Over Generations')
plt.xlabel('Generations')
plt.ylabel('Fitness Value')
plt.show()

