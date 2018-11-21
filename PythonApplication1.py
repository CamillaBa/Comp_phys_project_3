import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

def plot_format(orientation):
    if orientation =='random':
        return '--'
    if orientation =='oriented':
        return '-'

def color(T):
    if T=='1.000000':
        return 'blue'
    if T=='2.400000':
        return 'red'


def load_and_plot(filename, 
                  default =['E', 'Mabs','CV', 'X', ], 
                  analytic_sol = False,
                  L=20,
                  name = '',
                  max_cycle=101):
    data = {}
    data['E'], data['M'], data['CV'], data['X_'], data['X'], data['Mabs'] = np.loadtxt(filename, delimiter=',', unpack=True)
    n = len(data['E'],)
    analytic = {}
    analytic['E'] = np.ones(n)*(-7.98392)
    analytic['CV'] = np.ones(n)*(0.128329)
    analytic['X'] = np.ones(n)*(15.97321)
    analytic['Mabs'] = np.ones(n)*(3.994643)

    ylabels = {}
    ylabels['E']=r"$\langle E \rangle\quad   [J]$"
    ylabels['CV']=r"$C_V \quad   [ J T^{-1} ]$" 
    ylabels['X']=r"$\chi \quad   [ J ]$"
    ylabels['Mabs']=r"$\langle| \mathscr{M}|  \rangle$"
    ylabels['M']=r"$\langle \mathscr{M} \rangle$"

    labels = {}
    labels['E'] = r'$E$'
    labels['Mabs'] = r'$|\mathscr{M}|$'
    labels['X'] = r'$\chi$'
    labels['CV'] = r'$C_V$'
    labels['M']=r'$\mathscr{M}$'

    for quantity in default:
        plt.figure("expectation_value_{}_L_{}".format(quantity,L))
        plt.plot(data[quantity][0:max_cycle],'-',label=(r'$\langle$'+labels[quantity]+r'$\rangle$ '+name))
        if analytic_sol == True:
            plt.plot(analytic[quantity][0:max_cycle],'-',label=(r'$\langle$'+labels[quantity]+r'$\rangle$ analytic'))
        plt.xlabel('Monte Carlo cycles')
        plt.ylabel(ylabels[quantity])
        plt.title(r"$L=${}".format(L))
        plt.legend(loc="best")
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.show(block=False)


#=================================================================================
# Exercise b)
#=================================================================================

#load_and_plot('small_lattice_2_2__cycles_100000_T_1.000000.txt',
#              analytic_sol=True,
#              default =['E','Mabs'],
#              max_cycle=201,
#              L=2)

#load_and_plot('small_lattice_2_2__cycles_100000_T_1.000000.txt',
#              analytic_sol=True,
#              default =['CV'],
#              max_cycle=1251,
#              L=2)

#load_and_plot('small_lattice_2_2__cycles_100000_T_1.000000.txt',
#              analytic_sol=True,
#              default =['X'],
#              max_cycle=15001,
#              L=2)

#=================================================================================
# Exercise c)
#=================================================================================

def load_and_plot_T(default =['E', 'Mabs','CV', 'X', ], 
                    max_cycle=101,
                    L=20,
                    orientation = 'random',
                    T = '1.000000'):
    data = {}
    data['E'], data['M'], data['CV'],  data['X'], data['X_'], data['Mabs'] = np.loadtxt(orientation+
                                                                           '_lattice_'+str(L)+'_'+str(L)+'__'+
                                                                           'cycles_100000_'+
                                                                           'T_'+T+'.txt', delimiter=',', unpack=True)
    n = len(data['E'],)
    analytic = {}

    ylabels = {}
    ylabels['E']=r"$\langle E \rangle\quad   [J]$"
    ylabels['CV']=r"$C_V \quad   [ J T^{-1} ]$" 
    ylabels['X']=r"$\chi \quad   [ J ]$"
    ylabels['Mabs']=r"$\langle| \mathscr{M}|  \rangle$"

    labels = {}
    labels['E'] = r'$E$'
    labels['Mabs'] = r'$|\mathscr{M}|$'
    labels['X'] = r'$\chi$'
    labels['CV'] = r'$C_V$'

    for quantity in default:
        plt.figure("expectation_value_{}_L_{}".format(quantity,L))
        plt.plot(data[quantity][0:max_cycle],plot_format(orientation),color = color(T),label=(r'$\langle$'+labels[quantity]+r'$\rangle$ '+orientation+r' $\quad \frac{TK}{J}=$'+str(T)))
        plt.xlabel('Monte Carlo cycles')
        plt.ylabel(ylabels[quantity])
        plt.title(r"$L=${}".format(L))
        plt.legend(loc="best")
        plt.show(block=False)

# temperature 1
#for orientation in ['oriented','random']:
#    for temperature in ['1.000000','2.400000']:
#        load_and_plot_T(default = ['E',],
#                        max_cycle = 7001,
#                        L=20,
#                        T = temperature,
#                        orientation = orientation)

#for orientation in ['oriented','random']:
#    for temperature in ['1.000000','2.400000']:
#        load_and_plot_T(default = ['Mabs',],
#                        max_cycle =35001,
#                        L=20,
#                        T = temperature,
#                        orientation = orientation)

# accepted configurations vs cycles


#plt.figure("passes_vs_orientation_temperature")
#data = {}
#for orientation in ['oriented','random']:
#    for temperature in ['1.000000','2.400000']:
#        data[orientation+'_'+temperature] = np.loadtxt('passes_'+orientation+'_lattice_20_20__cycles_100000_T_'+temperature+'.txt', delimiter=',', unpack=True)
#        plt.plot(data[orientation+'_'+temperature],plot_format(orientation),color = color(temperature),label=orientation+r'$\quad \frac{TK}{J}=$'+temperature)
#        plt.xlabel('Monte Carlo cycles')
#        plt.ylabel('Passing spin configurations')
#        plt.title(r"$L=$20")
#        plt.legend(loc="best")
#        plt.show(block=False)
#        plt.yscale('log')



#=================================================================================
# Exercise d)
#=================================================================================

#for T in ['1.000000','2.400000']:
#    E = np.loadtxt('random_energy_lattice_20_20__cycles_1000000_T_'+T+'.txt', delimiter=',', unpack=True)
#    E=E[200000:]
#    max_E = np.max(E)
#    min_E = np.min(E)
#    energies = np.arange(min_E,max_E,4)
#    dist = np.zeros(len(energies))
#    for i, target_energy in enumerate(energies):
#        for energy in E:
#            if energy == target_energy:
#                dist[i]+=1
#    fig = plt.figure('energy_distribution'+ T)
#    plt.plot(energies,dist,label=r'$\frac{TK}{J}=$'+T,color=color(T))
#    plt.xlabel('Energy [J]')
#    plt.ylabel('Frequency')
#    plt.legend(loc="best")
    
#    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

#    m = np.mean(E)
#    std = np.std(E)
#    plt.title('Energy distribution (800000 cycles)\n'+r'$\langle E \rangle =$'+str(round(m,3))+r'$\quad\sigma_E^2=$'+str(round(np.var(E),3)))

#    plt.errorbar(m, np.max(dist)/2, xerr=std , fmt='o',color=color(T))


#    plt.vlines(m+std, 0, np.max(dist),linestyle='--',color=color(T))
#    plt.vlines(m-std, 0, np.max(dist),linestyle='--',color=color(T))

#    print(m,'  ',np.var(E))

#=================================================================================
# Exercise e)
#=================================================================================

# check that the values converge
#load_and_plot('random_lattice_100_100__cycles_1500000_T_2.250000.txt',
#              max_cycle=1500000,
#              default =['E', 'M','Mabs','CV', 'X', ],
#              L=100)

##load_and_plot('random_lattice_40_40__cycles_500000_T_2.250000.txt',
#              default =['E', 'M','Mabs','CV', 'X', ],
#              max_cycle=500000,
#              L=40)


# plots
ylabels = {}
ylabels['E']=r"$\langle E/L^2 \rangle\quad   [J]$"
ylabels['CV']=r"$C_V/L^2 \quad   [ J T^{-1} ]$" 
ylabels['X']=r"$\chi/L^2 \quad   [ J ]$"
ylabels['Mabs']=r"$\langle| \mathscr{M}|/L^2  \rangle$"

labels = {}
labels['E'] = r'$E$'
labels['Mabs'] = r'$|\mathscr{M}|$'
labels['X'] = r'$\chi$'
labels['CV'] = r'$C_V$'

for L in [40,60,80,100]:
    L2 = L*L;
    data_T = {}; times = []
    data_T['E'], data_T['M'], data_T['CV'], data_T['X_'], data_T['X'], data_T['Mabs'] = [], [], [], [], [], []
    for T in ['2.000000', '2.050000', '2.100000', '2.150000', '2.200000', '2.250000','2.300000','2.350000']:
        data  = {}
        data['E'], data['M'], data['CV'], data['X'], data['X_'], data['Mabs'] = np.loadtxt('random_lattice_'+str(L)+'_'+str(L)+'__'+
                                                                           'cycles_1500000_'+
                                                                           'T_'+T+'.txt', delimiter=',', unpack=True)
        
        T_float = float(T);

        if  T ==  '2.250000':
            data_T['X'].append(sum(data['X'][1499900:])/(100*T_float*L2))
        else:
            data_T['X'].append(sum(data['X'][1499900:])/(100*T_float*L2))

        for quantity in ['E','Mabs']:
            data_T[quantity].append(sum(data[quantity][1499900:])/(100*L2))

        data_T['CV'].append(sum(data['CV'][1499900:])/(100*T_float*T_float*L2))

        times.append(T_float)

    for quantity in ['E','CV','X','Mabs']:
        plt.figure('expectation_value_'+quantity+'_versus_time')
        plt.plot(times,data_T[quantity],'o--',label=r'$L=$'+str(L))
        plt.xlabel('T')
        plt.ylabel(ylabels[quantity])
        plt.legend(loc="best")
        plt.show(block=False)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.show()