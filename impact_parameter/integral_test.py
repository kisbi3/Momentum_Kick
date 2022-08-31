import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import multiprocessing
import time

# 그냥 아무생각없이 모두 cuda로 때려버렸음
# 이거 절대 안돌아감. 나중에 numpy를 일부 써야 할 것으로 보임.

time_start = time.time()

#constants
RA = 0.8        #fm
RB = 0.8        #fm
RR = 1.59       #fm
KP = 367.       #ratio
C_CMS = 0.269   #ratio
# SIGMA = 11.     #fm^2
SIGMA = 14.     #RHIC
# T0 = 0.43       #fm/c
T0 = 0.6        #RHIC
JETA = 0.20     #in eq(36)

#쪼개는 개수
nb = 53
nphis = 1000
nb0 = 1000
nbeta = 1000
nl_step = 500
nbp = 100
nalpha = 100


#쪼개는 간격
db = 1.59/nb
dphis = np.pi/(2*nphis)
db0 = 2./nb0
dl = (2.-0.)/nl_step
dbeta = cp.pi/(2*nbeta)
# np.set_printoptions(precision = 50)

#함수 선언
def Nk_b(impact_param):
    result = 0.
    phis = np.linspace(0., np.pi/2, nphis)
    pool = multiprocessing.Pool(processes=1)
    b_array = np.zeros(len(phis))-impact_param
    result_dist = pool.starmap(Nk_phib, zip(phis, b_array))
    # print(result_dist)
    result = np.sum(result_dist)*dphis*(2/np.pi)
    return result


def Nk_phib(phis_host, b_scalar_host):   # scalar, vector
    beta, l_scalar = cp.meshgrid(cp.linspace(0., cp.pi/2, nbeta), cp.linspace(0., 2., nl_step))
    b0_scalar_gpu = cp.linspace(0., 2., nb0)
    # beta_host = cp.asnumpy(beta)
    phis = cp.asarray(phis_host)
    b_scalar = cp.asarray(b_scalar_host)

    numerator = 0
    denominator = 0
    # dist_result = 0
    for b0 in b0_scalar_gpu:
        top = Nk_b0phib(b0, beta, l_scalar+dl, phis, b_scalar)
        bot = Nk_b0phib(b0, beta, l_scalar, phis, b_scalar)
        Nk_dist = cp.sum((top+bot)*dl/2, axis=0)
        
        # Nk_dist = cp.sum(Nk_b0phib(b0, beta, l_scalar, phis, b_scalar), axis=0)
        # print("\n Nk : ", Nk_dist, "\n Pjet : ", Pjet(b_scalar, b0, beta))
        same = b0*cp.exp(-JETA*Nk_dist)*Pjet(b_scalar, b0, beta)*db0*dbeta
        # print(same)
        numerator += cp.sum(Nk_dist*same)
        denominator += cp.sum(same)
        # 확률 함수에 의해 0이 나와서 nan이 등장하는 것으로 보임.
        # print("\n b0 : ", b0, "\n Nk_dist : ", Nk_dist, "\n exponential : ", cp.exp(-JETA*Nk_dist),"\nProbability : ", Pjet(b_scalar, b0, beta), "\n denominator ; ",  same)
    denominator = cp.where(denominator==0, 1, denominator)
    numerator = cp.where(denominator==0, 0, numerator)
    numerator_host = cp.asnumpy(numerator)
    denominator_host = cp.asnumpy(denominator)
    return numerator_host/denominator_host





def Nk_b0phib(b0_scalar, beta, scalar_l, phis, b):   #scalar, scalar, scalar
    vec_bpA_x = b0_scalar*cp.cos(beta) + scalar_l*cp.cos(phis) + b/2
    vec_bpA_y = b0_scalar*cp.sin(beta) + scalar_l*cp.sin(phis)
    vec_bpB_x = b0_scalar*cp.cos(beta) + scalar_l*cp.cos(phis) - b/2
    vec_bpB_y = b0_scalar*cp.sin(beta) + scalar_l*cp.sin(phis)
    bpA_size = vec_bpA_x*vec_bpA_x + vec_bpA_y*vec_bpA_y
    bpB_size = vec_bpB_x*vec_bpB_x + vec_bpB_y*vec_bpB_y

    # step_bpA = cp.clip(RA-cp.sqrt(bpA_size), 10**(-50), )/cp.absolute(RA-cp.sqrt(bpA_size))
    # step_bpB = cp.clip(RB-cp.sqrt(bpB_size), 10**(-50), )/cp.absolute(RB-cp.sqrt(bpB_size))
    step_bpA = cp.where((RA-cp.sqrt(bpA_size))<0, 0, 1)
    step_bpB = cp.where((RB-cp.sqrt(bpB_size))<0, 0, 1)
    # print("\n beta : ", beta, "\n phi : ", phis, "\n b0_scalar : " , b0_scalar, "\n scalar_l : " , scalar_l, "\n b : " , b)
    # print("\n TA : ", TT(RA, vec_bpA_x, vec_bpA_y), "\n TB : ", TT(RB, vec_bpB_x, vec_bpB_y))
    # print("\n TA : ", TT(RA, vec_bpA_x, vec_bpA_y), "\n TB : ", TT(RB, vec_bpB_x, vec_bpB_y), "\n Nkkkkk... : " , step_bpA*step_bpB*SIGMA*KP*(TT(RA, vec_bpA_x, vec_bpA_y) + TT(RB, vec_bpB_x, vec_bpB_y))/(2*(T0+scalar_l)))
    return step_bpA*step_bpB*SIGMA*KP*(TT(RA, vec_bpA_x, vec_bpA_y) + TT(RB, vec_bpB_x, vec_bpB_y))/(2*(T0+scalar_l))

def N_ch(impact_param):
    b = -impact_param
    bp_scalar, alpha = np.meshgrid(np.linspace(0.,1.59, nbp), np.linspace(0., np.pi*2, nalpha))
    # bpx = bp_scalar*np.cos(alpha)
    # bpy = bp_scalar*np.sin(alpha)

    bpxA = bp_scalar*np.cos(alpha)+b/2
    bpyA = bp_scalar*np.sin(alpha)
    bpxB = bp_scalar*np.cos(alpha)-b/2
    bpyB = bp_scalar*np.sin(alpha)
    bpA_size = bpxA*bpxA + bpyA*bpyA
    bpB_size = bpxB*bpxB + bpyB*bpyB
    step_bA = np.heaviside(RA-np.sqrt(bpA_size), 0)
    step_bB = np.heaviside(RB-np.sqrt(bpB_size), 0)
    # result_dist = np.sum(C_CMS*(2/3)*KP*(TT_host(RA, bpx+b/2, bpy)+TT_host(RB, bpx-b/2, bpy)))
    result_dist = np.sum(step_bA*step_bB*C_CMS*(2/3)*KP*(TT_host(RA, bpxA, bpyA)+TT_host(RB, bpxB, bpyB)))
    return result_dist*(1.59/nbp)*(2*np.pi/nalpha)

def TT_host(R, bx, by):
    square_b = bx*bx+by*by
    stepfunction = np.heaviside(R-np.sqrt(square_b),0)
    squareroot = np.sqrt((R*R-square_b)*stepfunction)
    return (3*squareroot)/(2*np.pi*R*R*R)

def TT(R, bx, by):  # gpu
    square_b = bx*bx+by*by
    return (3*cp.sqrt(cp.clip(R*R-square_b, 0, )))/(2*cp.pi*R*R*R)

    
def Pjet(b, scalar_b0, beta):
    b0x = scalar_b0*cp.cos(beta) 
    b0y = scalar_b0*cp.sin(beta)
    return TT(RA, b0x+b/2, b0y)*TT(RB, b0x-b/2, b0y)

def Pjet_norm_func(b_gpu, b0_norm, beta_norm):
    dist = np.sum(Pjet(-b_gpu, b0_norm, beta_norm) * db0 * dbeta, axis=0)
    return np.asnumpy(np.sum(dist, axis=1))



if __name__ == '__main__':
    impact_param = np.linspace(0.,1.59, nb)
    # Nk_result = Nk_b(impact_param)

    Nk_result = np.zeros(nb+1)
    Nch_result = np.zeros(nb+1)
    # Nk_result = np.zeros(nb)
    # print(int((1.59-0.)/db)+1)

    # pool = multiprocessing.Pool(processes=16)
    # Nk_result = pool.starmap(Nk_b, zip(impact_param))

    Nk_result = Nk_b(0.03)
    exit(1)

    # print("impact parameter \t Nk \t Nch")
    # for i in range(nb+1):
    #     Nch_result[i] = N_ch(impact_param[i])
    #     Nk_result[i] = Nk_b(impact_param[i])
    #     print(impact_param[i], Nk_result[i], Nch_result[i])

    # Graphs


    fig = plt.figure()
    ax = plt.axes()
    fig.set_size_inches(35, 16.534, forward=True)


    plt.plot(impact_param, Nk_result, linewidth=7, label = r'$N_k$')
    # plt.plot(impact_param, Nch_result, linewidth=7, label = r'$N_ch$')

    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))
    plt.tick_params(axis='both',which='major',direction='in', width=2,length=30,labelsize=45, top='true')
    plt.tick_params(axis='both',which='minor',direction='in', width=2,length=15,labelsize=45, top='true')
    plt.xlabel(r'$impact \,\, parameter$',size=70)
    # plt.ylabel(r'$N_ch, N_k $',size=70)

    plt.grid(color='silver',linestyle=':',linewidth=5)
    plt.legend(fontsize=45,framealpha=False,loc='upper left')

    plt.tight_layout()


    plt.legend()
    # plt.show()
    fig.savefig('impact_multi_cuda.png')
    fig.clear()

    plt.plot(Nch_result, np.nan_to_num(Nk_result), linewidth=7)
    # plt.legend()
    # plt.show()

    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '${:g}$'.format(y)))
    plt.tick_params(axis='both',which='major',direction='in',width=2,length=30,labelsize=45, top='true')
    plt.tick_params(axis='both',which='minor',direction='in',width=2,length=15,labelsize=45, top='true')
    plt.grid(color='silver',linestyle=':',linewidth=5)
    plt.xlabel(r'$ N_{ch} $',size=70)
    plt.ylabel(r'$ N_k $',size=70)
    plt.tight_layout()

    fig.savefig('NNN_multi_cuda.png')


    time_end = time.time()
    print("time : ",time_end-time_start, "sec")