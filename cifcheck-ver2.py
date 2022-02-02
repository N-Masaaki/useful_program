from pymatgen.io.cif import CifParser
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import re
import BVSparameter

# A='cifファイル名', B=Cation原子番号,C=Anion原子番号
# D = 対象Cationの価数, E = Cation.価数.Anion.価数.reference (例：'Pb2O-2bs')
def cifcheck(A,B,C,D,E):
    r0,b = BVSparameter.bvsdict(E)
    plt.rcParams["font.size"] = 16
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.direction'] = 'in'
    name1 = "{}"
    name = name1.format(A)  #組成名
    parser = CifParser(name,occupancy_tolerance=1.0,site_tolerance=0.0001)
    structure = parser.get_structures()[0]
    structure.to(filename="poscar")
    num = structure.atomic_numbers
    data = np.genfromtxt(A,delimiter='\n',dtype=str)
    label1 = "{}-{}"
    label1 = label1.format(str(structure.species[num.index(B)]),str(structure.species[num.index(C)]))
    n = 0
    while n < len(data):
        if re.match(r'_chemical_formula_sum.*',data[n]):
            a = data[n]
        n += 1
    a0 = str(a)
    a1 = re.split('\s', a0)
    a1 = list(filter(lambda str:str != '', a1))
    del a1[0]
    name2 = ' '.join(a1).replace('\'', '')
    
    plt.figure(figsize=(6,4.5))
    for k in range(num.index(B),num.index(B)+num.count(B)):                         
        dis = np.array(structure.get_neighbors(structure[k],6))
        Bond = []
        for i in range(len(dis)):
            if '[PeriodicSite: %s' %str(structure.species[num.index(C)]) in str(dis[i]):
                Bond.append(dis[:,1][i])

        Bond.sort()
        BVS = []
        for j in range(len(Bond)):
            BVS.append(np.exp((r0-Bond[j])/b))
    
        plt.bar(Bond,BVS,width=0.01,color='b')
        for n in range(len(Bond)):
            if n < len(Bond)-1:
                BVS[n+1] = BVS[n] + BVS[n+1]
        
        plt.plot(Bond,BVS,'o')
        x = np.arange(int(Bond[0]),int(Bond[0])+2,0.01)
        y = np.exp((r0-x)/b)
        plt.plot(x,y,label=label1,color='b')
        plt.xlabel('Cation-O distance ($\mathrm{\AA}$)')
        plt.ylabel('Bond valence (-)')
        plt.title(name2)
        x1 = np.arange(int(Bond[0]),int(Bond[0])+2.1, 0.5)
        y1 = np.arange(0.0, D+1.1, 1.0)
        plt.xticks(x1)
        plt.yticks(y1)
        plt.gca().yaxis.set_minor_locator(tick.MultipleLocator(0.2))            # 小メモリの間隔（今回は0.2）
        plt.gca().xaxis.set_minor_locator(tick.MultipleLocator(0.1))   
        plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))     # y軸の小数点以下設定(今回は1桁)
        plt.xlim([int(Bond[0]),int(Bond[0])+2])
        plt.ylim([0.0, D+1.0])
        plt.tick_params(direction='in', length=5, bottom=True, top=True, left=True, right=True)
        plt.axhline(y=D,ls = '--')
        #plt.legend(prop={'size': 12})
        plt.tight_layout()
        plt.show()




def allcheck(A,B,C,D,E):
    r0,b = BVSparameter.bvsdict(E)
    plt.rcParams["font.size"] = 16
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.direction'] = 'in'
    name1 = "{}"
    name = name1.format(A)  #組成名
    parser = CifParser(name,occupancy_tolerance=1.0,site_tolerance=0.0001)
    structure = parser.get_structures()[0]
    structure.to(filename="poscar")
    num = structure.atomic_numbers
    data = np.genfromtxt(A,delimiter='\n',dtype=str)
    label1 = "{}-{}"
    label1 = label1.format(str(structure.species[num.index(B)]),str(structure.species[num.index(C)]))
    n = 0
    while n < len(data):
        if re.match(r'_chemical_formula_sum.*',data[n]):
            a = data[n]
        n += 1
    a0 = str(a)
    a1 = re.split('\s', a0)
    a1 = list(filter(lambda str:str != '', a1))
    del a1[0]
    name2 = ' '.join(a1).replace('\'', '')
    
    plt.figure(figsize=(6,4.5))
    for k in range(num.index(B),num.index(B)+num.count(B)):                         
        dis = np.array(structure.get_neighbors(structure[k],4))
        Bond = []
        for i in range(len(dis)):
            if '[PeriodicSite: %s' %str(structure.species[num.index(C)]) in str(dis[i]):
                Bond.append(dis[:,1][i])


        Bond.sort()
        BVS = []
        for j in range(len(Bond)):
            BVS.append(np.exp((r0-Bond[j])/b))
    
        plt.bar(Bond,BVS,width=0.01,color='b')
        for n in range(len(Bond)):
            if n < len(Bond)-1:
                BVS[n+1] = BVS[n] + BVS[n+1]
        
        plt.plot(Bond,BVS,'o')
    x = np.arange(int(Bond[0]),int(Bond[0])+2,0.01)
    y = np.exp((r0-x)/b)
    plt.plot(x,y,label=label1,color='b')
    plt.xlabel('Cation-O distance ($\mathrm{\AA}$)')
    plt.ylabel('Bond valence (-)')
    plt.title(name2)
    x1 = np.arange(int(Bond[0]),int(Bond[0])+2.1, 0.5)
    y1 = np.arange(0.0, D+1.1, 1.0)
    plt.xticks(x1)
    plt.yticks(y1)
    plt.gca().yaxis.set_minor_locator(tick.MultipleLocator(0.2))            # 小メモリの間隔（今回は0.2）
    plt.gca().xaxis.set_minor_locator(tick.MultipleLocator(0.1))   
    plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))     # y軸の小数点以下設定(今回は1桁)
    plt.xlim([int(Bond[0]),int(Bond[0])+2])
    plt.ylim([0.0, D+1.0])
    plt.tick_params(direction='in', length=5, bottom=True, top=True, left=True, right=True)
    plt.axhline(y=D,ls = '--')
    plt.legend(prop={'size': 12})
    plt.tight_layout()
    plt.show()


