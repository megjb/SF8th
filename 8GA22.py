import random as ra
import math as ma
import matplotlib.pyplot as plt
## Initializes constants and variables  
crossover_rate = .9
mutation_rate = .4
accuracy_rate = .75
gen_count = 200
best_fitness = 0
best_chromosome = []
pop = []
popcount = 750
grid = {}
length = 5
breadth = 5
hlayerend = [0,0,0,0]
hlayer1  = [0 for i in range(length*breadth*2)]
def sig(x):
    return 1/(1+ma.exp(-x))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
indiv_coord=[0,0]
inputnum = 4
chrom_input = {}
enem_coords = 5,5 

## Creates the blueprint for the whole population
def chrom_initiation():
    indiv=[[ra.uniform(-4,4) for i in range(length*breadth*inputnum)] for j in hlayer1], [[ra.uniform(-4,4) for y in hlayer1] for z in hlayerend] 
    return indiv

## Creates the population of individuals
def pop_initiation(pop_size):
    global pop
    pop = [chrom_initiation() for x in range (1,pop_size) ]
    return pop

## Creates the input map for feach individual's choice algorithm
for i in range(length):
    chrom_input[i]={}
    for j in range(breadth):
        chrom_input[i][j]=[0 for i in range(inputnum)]

## Creates grid where battle takes place
for i in range(length):
    grid[i] = {}
    for j in range(breadth):
        grid[i][j]=[1,0,[0,0]]

## Cross references points on grid to map
def crossreference(grid,chrom_input):    
    for i in range(length):
        for j in range(breadth):
            for z in grid[i][j]:
                if z == 1:
                    chrom_input[i][j][z] = 1
                elif z==0:
                    chrom_input[i][j][z] = 0    

## Creates new members of the population based off of the previous gen. Either through duplication with some mutation, or crossover
def new_pop(pop,popfitness):
    newpop = []
    while len(newpop) != len(pop):
        new_chromo = []
        if crossover_rate >= ra.uniform(0,1):
            new_chromo = crossover(pop,popfitness)
            if len(newpop)+2 > len(pop):
                while len(newpop)+2 > len(pop):
                    choi = ra.choice(new_chromo)
                    new_chromo.remove(choi)
                newpop.append(new_chromo)    
            newpop.append(new_chromo[0])
            newpop.append(new_chromo[1])
            return newpop
        else:
            new_chromo = mutation(roulette_selection(pop,popfitness)) 
            if len(newpop)+1 > len(pop):
                    return newpop
            newpop.append(new_chromo)
    return newpop

## Selects a random individual in a weighted roulette wheel selection
def roulette_selection(pop,popfitness):
    global count 
    wheel = []
    for x in pop:
        indices = pop.index(x)
        number = round(popfitness[indices]) 
        if number<= 0:
            if max(popfitness) <= 0:
                numbers = popfitness
                sorted_numbers = sorted(numbers)
                number = sorted_numbers.index(round(popfitness[indices]))
                if number == 0:
                    number =1
        for i in range(number):
            wheel.append(x)
    choice = ra.choice(wheel)
    return choice

## Crossovers the individuals. Is the basis for the reproduction system
def crossover(pop,popfitness):
    a = roulette_selection(pop,popfitness)
    b = roulette_selection(pop,popfitness)
    p1 = ra.randint(1,hlayer1)
    p2 = ra.randint(1,hlayer1-p1)
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
        p4 = ra.randint(1,(length*breadth*inputnum)-p3)
        c[0][hlayer1.index(i)] = a[0][hlayer1.index(i)][:p3] + b[0][hlayer1.index(i)][p3:p4] + a[0][hlayer1.index(i)][p4:] 
        d[0][hlayer1.index(i)] = b[0][hlayer1.index(i)][:p3] + a[0][hlayer1.index(i)][p3:p4] + b[0][hlayer1.index(i)][p4:] 
    
    for i in hlayerend:
        p3 = ra.randint(1,hlayer1)
        p4 = ra.randint(1,(hlayer1)-p3)
        c[1][hlayerend.index(i)] = a[1][hlayerend.index(i)][:p3] + b[1][hlayerend.index(i)][p3:p4] + a[1][hlayerend.index(i)][p4:] 
        d[1][hlayerend.index(i)] = b[1][hlayerend.index(i)][:p3] + a[1][hlayerend.index(i)][p3:p4] + b[1][hlayerend.index(i)][p4:] 
    return c, d

## Mutates the individual.
def mutation(indiv):
    for i in indiv:
        if mutation_rate >= ra.uniform(0,1):
            choice1 = ra.choice(indiv)
            choice2 = ra.choice(indiv)
            index1 = indiv.index(choice1)
            index2 = indiv.index(choice2)
            indiv.pop(index1)
            indiv.pop(index2)
            indiv.insert(index1,choice2)
            indiv.insert(index2,choice1)
        for j in indiv[0][i]:
            if mutation_rate >= ra.uniform(0,1):
                choice1 = ra.choice(indiv[0][i])
                choice2 = ra.choice(indiv[0][i])
                index1 = indiv[0][i].index(choice1)
                index2 = indiv[0][i].index(choice2)
                indiv[0][i].pop(index1)
                indiv[0][i].pop(index2)
                indiv[0][i].insert(index1,choice2)
                indiv[0][i].insert(index2,choice1)
        for j in indiv[1][i]:
            if mutation_rate >= ra.uniform(0,1):
                choice1 = ra.choice(indiv[1][i])
                choice2 = ra.choice(indiv[1][i])
                index1 = indiv[1][i].index(choice1)
                index2 = indiv[1][i].index(choice2)
                indiv[1][i].pop(index1)
                indiv[1][i].pop(index2)
                indiv[1][i].insert(index1,choice2)
                indiv[1][i].insert(index2,choice1)
    return indiv

def lanchester_battle(indiv,enem):
    enem_pop = enem[4][0]/2
    indiv_pop = indiv[4][0]/2
    
    while indiv_pop and enem_pop >= 1:
        
        pop1 = enem_pop
        pop2 = indiv_pop
        indiv_pop -= (pop1*(enem[4][1]*accuracy_rate))
        enem_pop -= (pop2*(indiv[4][1]*accuracy_rate))
 
        if enem_pop or indiv_pop <= 0:
            if enem_pop >= indiv_pop:
                winner = ('enem',enem_pop)
            else:
                winner = ('indiv',indiv_pop)

            break
    return winner
                
def war(indiv):
    grid[0][0]=[1,0,1,1[1500,.5]]
    for i in range(10):
        for i in range(length):
            for j in range(breadth):
                if grid[i][j][2] == 1:
                    troops_in_control = [i,j]
        chrom_input = crossreference(grid,chrom_input)
        dec = decision(chrom_input,indiv)
        output = max(dec)
        if dec.index(output) == 0:
            if grid[troops_in_control[0]+1][troops_in_control[1]][1] == 1:
                simulation = lanchester_battle(grid[troops_in_control[0]][troops_in_control[1]],grid[troops_in_control[0]+1][troops_in_control[1]])
                if simulation[0] == 'enem':
                    grid[troops_in_control[0]][troops_in_control[1]][4][0] = grid[troops_in_control[0]][troops_in_control[1]][4][0]/2
                    grid[troops_in_control[0]+1][troops_in_control[1]][4][0] = abs(simulation[1]-=(grid[troops_in_control[0]+1][troops_in_control[1]])[4][0]/2)
                    continue
                else:
                    grid[troops_in_control[0]+1][troops_in_control[1]][4][0] = grid[troops_in_control[0]][troops_in_control[1]][4][0]/2
                    grid[troops_in_control[0]][troops_in_control[1]][4][0] = abs(simulation[1]-=(grid[troops_in_control[0]][troops_in_control[1]])[4][0]/2)
                    continue
            grid[troops_in_control[0]+1][troops_in_control[1]] = [1,0,1,1[grid[troops_in_control[0]+1][troops_in_control[1]][4][0]+=grid[troops_in_control[0]][troops_in_control[1][4][0]/2],.5]]
            grid[troops_in_control[0]][troops_in_control[1]] = [1,0,1,0[grid[troops_in_control[0]][troops_in_control[1][4][0]/2],1]]
        if dec.index(output) == 1:
            if grid[troops_in_control[0]-1][troops_in_control[1]][1] == 1:
                simulation = lanchester_battle(grid[troops_in_control[0]][troops_in_control[1]],grid[troops_in_control[0]-1][troops_in_control[1]])
                if simulation[0] == 'enem':
                    grid[troops_in_control[0]][troops_in_control[1]][4][0] = grid[troops_in_control[0]][troops_in_control[1]][4][0]/2
                    grid[troops_in_control[0]-1][troops_in_control[1]][4][0] = abs(simulation[1]-=(grid[troops_in_control[0]-1][troops_in_control[1]])[4][0]/2)
                    continue
                else:
                    grid[troops_in_control[0]-1][troops_in_control[1]][4][0] = grid[troops_in_control[0]][troops_in_control[1]][4][0]/2
                    grid[troops_in_control[0]][troops_in_control[1]][4][0] = abs(simulation[1]-=(grid[troops_in_control[0]][troops_in_control[1]])[4][0]/2)
                    continue
            grid[troops_in_control[0]-1][troops_in_control[1]] = [1,0,1,1[grid[troops_in_control[0]-1][troops_in_control[1]][4][0]+=grid[troops_in_control[0]][troops_in_control[1][4][0]/2],.5]]
            grid[troops_in_control[0]][troops_in_control[1]] = [1,0,1,0[grid[troops_in_control[0]][troops_in_control[1][4][0]/2],1]]
        if dec.index(output) == 2:
            if grid[troops_in_control[0]][troops_in_control[1]+1][1] == 1:
                simulation = lanchester_battle(grid[troops_in_control[0]][troops_in_control[1]],grid[troops_in_control[0]][troops_in_control[1]+1])
                if simulation[0] == 'enem':
                    grid[troops_in_control[0]][troops_in_control[1]][4][0] = grid[troops_in_control[0]][troops_in_control[1]][4][0]/2
                    grid[troops_in_control[0]][troops_in_control[1]+1][4][0] = abs(simulation[1]-=(grid[troops_in_control[0]][troops_in_control[1]+1])[4][0]/2)
                    continue
                else:
                    grid[troops_in_control[0]][troops_in_control[1]+1][4][0] = grid[troops_in_control[0]][troops_in_control[1]][4][0]/2
                    grid[troops_in_control[0]][troops_in_control[1]][4][0] = abs(simulation[1]-=(grid[troops_in_control[0]][troops_in_control[1]])[4][0]/2)
                    continue
            grid[troops_in_control[0]][troops_in_control[1]+1] = [1,0,1,1[grid[troops_in_control[0]][troops_in_control[1]+1][4][0]+=grid[troops_in_control[0]][troops_in_control[1][4][0]/2],.5]]
            grid[troops_in_control[0]][troops_in_control[1]] = [1,0,1,0[grid[troops_in_control[0]][troops_in_control[1][4][0]/2],1]]
        if dec.index(output) == 3:
            if grid[troops_in_control[0]][troops_in_control[1]-1][1] == 1:
                simulation = lanchester_battle(grid[troops_in_control[0]][troops_in_control[1]],grid[troops_in_control[0]][troops_in_control[1]-1])
                if simulation[0] == 'enem':
                    grid[troops_in_control[0]][troops_in_control[1]][4][0] = grid[troops_in_control[0]][troops_in_control[1]][4][0]/2
                    grid[troops_in_control[0]][troops_in_control[1]-1][4][0] = abs(simulation[1]-=(grid[troops_in_control[0]][troops_in_control[1]-1])[4][0]/2)
                    continue
                else:
                    grid[troops_in_control[0]][troops_in_control[1]-1][4][0] = grid[troops_in_control[0]][troops_in_control[1]][4][0]/2
                    grid[troops_in_control[0]][troops_in_control[1]][4][0] = abs(simulation[1]-=(grid[troops_in_control[0]][troops_in_control[1]])[4][0]/2)
                    continue
            grid[troops_in_control[0]][troops_in_control[1]-1] = [1,0,1,1[grid[troops_in_control[0]][troops_in_control[1]-1][4][0]+=grid[troops_in_control[0]][troops_in_control[1][4][0]/2],.5]]
            grid[troops_in_control[0]][troops_in_control[1]] = [1,0,1,0[grid[troops_in_control[0]][troops_in_control[1][4][0]/2],1]]



def decision(chrom_input,indiv):
    for i in (chrom_input[x][y][j] for x in range(length) for y in range(breadth)for j in range(3)):
        for j in hlayer1:
            hlayer1[hlayer1.index(j)]=sig(j+(i*indiv[0][hlayer1.index(j)][list(grid).index(i)]))
    for i in hlayer1:
        for j in hlayerend:
            hlayerend[hlayerend.index(j)]=sig(j+(i*indiv[1][hlayerend.index(j)][list(hlayer1).index(i)]))        
    outputlayer = hlayerend
    return outputlayer
print(outputlayer)    
 