import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
import multiprocessing
import time

time_start = time.time()

'''
    impact parameter에 대한 N_k, N_ch 계산
'''

'''constants'''
# Radius = 0.8        #fm
# Diameter = 1.59       #fm
RA = 0.8
RB = 0.8
KP = 500.       #ratio 13 TeV
# KP = 367.       #ratio 7 TeV
# KP = 21.      #RHIC  200 GeV
# C_CMS = 0.269   #ratio 7 TeV
C_CMS = 0.252   #ratio 13 TeV
SIGMA = .14
T0 = 0.39        #fm/c  13 TeV
# T0 = 0.43       #fm/c  7 TeV
# T0 = 0.6        #RHIC 200 GeV
JETA = 0.50     #in eq(36) 13 TeV
# JETA = 0.20     #in eq(36) 7 TeV

#쪼개는 개수
nb = 100
nphis = 100
nb0 = 300
nl_step = 1000
nbp = 1000

#쪼개는 간격
db = 0.8*2/nb
dphis = np.pi/(2*nphis)
db0x = 0.8*2/nb0
db0y = 0.8*2/nb0
dl = (0.8*4-T0)/nl_step

# 0부터가 아니라 -0.8 ~ 0.8로 해야될거 같은데??

#함수 선언
def Nk_b(impact_param):
    result = 0.
    phis = np.linspace(0., np.pi/2, nphis)      # radian
    pool = multiprocessing.Pool(processes=16)
    b_array = np.zeros(len(phis))-impact_param
    result_dist = pool.starmap(Nk_phib, zip(phis, b_array))
    result = np.sum(result_dist)*dphis*(2/np.pi)
    return result


def Nk_phib(phis_host, b_scalar_host):   # scalar, vector
    b0y, l_scalar = cp.meshgrid(cp.linspace(-0.8, 0.8, nb0), cp.linspace(T0, 0.8*4, nl_step))
    b0x_arr = cp.linspace(-0.8, 0.8, nb0)
    # beta_host = cp.asnumpy(beta)
    phis = cp.asarray(phis_host)
    b_scalar = cp.asarray(b_scalar_host)

    numerator = 0
    denominator = 0
    for b0x in b0x_arr:
        Nk_dist = cp.sum(Nk_b0phib(b0x, b0y, l_scalar, phis, b_scalar)*dl, axis=0)
        Pjet_dist = Pjet(b_scalar, b0x, b0y[0])
        same = cp.exp(-JETA*Nk_dist)*Pjet_dist*db0x*db0y
        numerator += cp.sum(Nk_dist*same)
        denominator += cp.sum(same)
    numerator_host = cp.asnumpy(numerator) 
    denominator_host = cp.asnumpy(denominator)
    if denominator_host == 0:
        return 0
    else:
        return numerator_host/denominator_host


def Nk_b0phib(b0x, b0y, scalar_l, phis, b):   #scalar, scalar, scalar
    vec_bpA_x = b0x + scalar_l*cp.cos(phis) + b/2;    vec_bpA_y = b0y + scalar_l*cp.sin(phis)
    vec_bpB_x = b0x + scalar_l*cp.cos(phis) - b/2;    vec_bpB_y = b0y + scalar_l*cp.sin(phis)
    TA = TT(RA, vec_bpA_x, vec_bpA_y);    TB = TT(RB, vec_bpB_x, vec_bpB_y)
    step_bpA = cp.where(TA == 0, 0, 1);    step_bpB = cp.where(TB == 0, 0, 1)
    return step_bpA*step_bpB*SIGMA*KP*(TA + TB)/(2*(scalar_l))

def N_ch(impact_param):
    b = -impact_param
    bp_range = 0.8
    bpx, bpy = np.meshgrid(np.linspace(-bp_range,bp_range, nbp), np.linspace(-bp_range,bp_range, nbp))
    # print(np.shape(bpx), np.shape(bpy))
    dbpx = 2*bp_range/nbp;    dbpy = 2*bp_range/nbp
    bpxA = bpx+b/2;    bpyA = bpy;    bpxB = bpx-b/2;    bpyB = bpy
    TA = TT_host(RA, bpxA, bpyA);    TB = TT_host(RB, bpxB, bpyB)
    step_bA = np.where(TA==0, 0, 1);    step_bB = np.where(TB==0, 0, 1)
    result_dist = C_CMS*(2/3)*KP*np.sum(step_bA*step_bB*(TA+TB))*dbpx*dbpy
    return result_dist


def TT_host(R, bx, by):
    square_b = bx*bx+by*by
    stepfunction = np.heaviside(R-np.sqrt(square_b),0)
    squareroot = np.sqrt((R*R-square_b)*stepfunction)
    return (3*squareroot)/(2*np.pi*R*R*R)

def TT(R, bx, by):  # gpu
    square_b = bx*bx+by*by
    return (3*cp.sqrt(cp.clip(R*R-square_b, 0, )))/(2*cp.pi*R*R*R)

    
def Pjet(b, b0x, b0y):
    return TT(RA, b0x+b/2, b0y)*TT(RB, b0x-b/2, b0y)


if __name__ == '__main__':
    # High multiplicity인 Nk를 저장해서 평균해가지고 논문에 쓸 변수
    HighMulti_Average = 0
    HighMulti_Average_number = 0
    mpl.rcParams["text.usetex"] = True
    impact_param = np.linspace(0.,1.6, nb)
    Nk_result = np.zeros(nb)
    Nch_result = np.zeros(nb)
    print("Strating Calculate")

    print("impact parameter \t Nk \t Nch")
    for i in range(nb):
        Nch_result[i] = N_ch(impact_param[i])
        Nk_result[i] = Nk_b(impact_param[i])
        print(impact_param[i], Nk_result[i], Nch_result[i])
        if(Nch_result[i]>105):
            HighMulti_Average = HighMulti_Average + Nk_result[i]
            HighMulti_Average_number = HighMulti_Average_number + 1
        # Nch_result[i] = N_ch(impact_param[i])
        # print(impact_param[i], Nch_result[i])
    print(f"High multiplicity Average: {HighMulti_Average/HighMulti_Average_number}")
    # Graphs
    print("End")

    print("Saving File...")
    np.savetxt("impact_param_result.csv", list(zip(impact_param, Nk_result, Nch_result)), delimiter=',')
    print("End")

    print("Drawing Graphs...")

    mpl.rcParams["text.usetex"] = True

    fig = plt.figure()
    ax = plt.axes()
    fig.set_size_inches(35, 16.534, forward=True)

    impact_param = np.loadtxt("impact_param_result.csv", usecols=[0], delimiter=',')
    Nk_result = np.loadtxt("impact_param_result.csv", usecols=[1], delimiter=',')
    N_ch_result = np.loadtxt("impact_param_result.csv", usecols=[2], delimiter=',')

    plt.plot(impact_param, Nk_result, linewidth=7, label = r'$N_k$')
    # plt.plot(impact_param, Nch_result, linewidth=7, label = r'$N_ch$')

    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))
    plt.tick_params(axis='both',which='major',direction='in', width=2,length=30,labelsize=45, top='true')
    plt.tick_params(axis='both',which='minor',direction='in', width=2,length=15,labelsize=45, top='true')
    plt.xlabel(r'$impact \,\, parameter$',size=70)
    plt.ylabel(r'$\langle N_k \rangle$',size=70)

    plt.grid(color='silver',linestyle=':',linewidth=5)
    # plt.legend(fontsize=45,framealpha=False,loc='upper left')

    plt.tight_layout()


    plt.legend()
    # plt.show()
    fig.savefig('Nk.png')
    fig.clear()

    plt.plot(impact_param, Nch_result, linewidth=7, label = r'$N_ch$')
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))
    plt.tick_params(axis='both',which='major',direction='in', width=2,length=30,labelsize=45, top='true')
    plt.tick_params(axis='both',which='minor',direction='in', width=2,length=15,labelsize=45, top='true')
    plt.xlabel(r'$impact \,\, parameter$',size=70)
    plt.ylabel(r'$N_{ch}$',size=70)

    plt.grid(color='silver',linestyle=':',linewidth=5)
    # plt.legend(fontsize=45,framealpha=False,loc='upper left')

    plt.tight_layout()


    plt.legend()
    # plt.show()
    fig.savefig('Nch.png')
    fig.clear()


    plt.plot(Nch_result, np.nan_to_num(Nk_result), linewidth=7)
    # plt.legend()
    # plt.show()

    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))
    plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top='true')
    plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top='true')
    plt.grid(color='silver',linestyle=':',linewidth=5)
    plt.xlabel(r'$ N_{ch} $',size=70)
    plt.ylabel(r'$\langle N_k \rangle$',size=70)
    plt.tight_layout()

    fig.savefig('Nch_Nk.png')

    print("End")



    time_end = time.time()
    print("time : ",time_end-time_start, "sec")