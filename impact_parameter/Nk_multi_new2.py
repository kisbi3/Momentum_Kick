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
RB = 0.79       #fm
RR = 1.59       #fm
KP = 367.       #ratio
C_CMS = 0.269   #ratio
SIGMA = 11.     #fm^2
T0 = 0.43       #fm/c
JETA = 0.20     #in eq(36)

#쪼개는 간격
db = 0.03
dphis = 0.005
db0 = 0.005
dbeta = 0.005

#쪼개는 개수
nb = 100
nphis = 100
nb0 = 100
nbeta = 100
nl_step = 500
nbp = 1000
nalpha = 1000

np.set_printoptions(precision = 50)

#함수 선언
def Nk_b(impact_param):
    result = 0.
    dphis = np.pi/(2*nphis)
    phis = np.arange(0., np.pi/2, np.pi/(2*nphis))

    pool = multiprocessing.Pool(processes=16)
    b_array = np.zeros(len(phis))-impact_param
    result_dist = pool.starmap(Nk_phib, zip(phis, b_array))
    # print(result_dist)
    result = np.sum(result_dist)*dphis*(2/np.pi)
    # result=Nk_phib(0., -impact_param)

    # for i in range(nphis):
    #     phis = np.pi*i/(2*nphis)
    #     result += (2*np.sum(Nk_phib(phis, -impact_param))*np.pi/(2*nphis))/np.pi
    print(result)
    return result


def Nk_phib(phis_host, b_scalar_host):   # scalar, vector
    beta, l_scalar = cp.meshgrid(cp.arange(0., cp.pi/2, cp.pi/(2*nbeta)), cp.arange(0.43, 1., (1.-0.43)/nl_step))
    b0_scalar_gpu = cp.arange(0., 0.79, 0.79/nb0)
    # beta_host = cp.asnumpy(beta)
    phis = cp.asarray(phis_host)
    b_scalar = cp.asarray(b_scalar_host)

    numerator = 0
    denominator = 0
    # dist_result = 0
    for b0 in b0_scalar_gpu:
        top = Nk_b0phib(b0, beta, l_scalar+(1.-0.43)/nl_step, phis, b_scalar)
        bot = Nk_b0phib(b0, beta, l_scalar, phis, b_scalar)
        Nk_dist = cp.sum((top+bot)*((1.-0.43)/nl_step)/2, axis=0)
        print(Nk_dist)
        # print(Nk_dist, cp.exp(-JETA*Nk_dist))
        same = b0*cp.exp(-JETA*Nk_dist)*Pjet(b_scalar, b0, beta)*(0.79/nb0)*(cp.pi/(2*nbeta))
        # print(b0_scalar_gpu[i], Nk_dist, same, same*Nk_dist)
        numerator += cp.sum(Nk_dist*same)
        denominator += cp.sum(same)
        # print(b0_scalar_gpu[i], numerator, denominator)

    numerator_host = cp.asnumpy(numerator)
    denominator_host = cp.asnumpy(denominator)
    # print(b_scalar_host, numerator_host, denominator_host, numerator_host/denominator_host)
    return numerator_host/denominator_host
    # return numerator_host       #그냥 분자만 내보내게 해보자.


def Nk_b0phib(b0_scalar, beta, scalar_l, phis, b):   #scalar, scalar, scalar
    vec_bpA_x = b0_scalar*cp.cos(beta) + scalar_l*cp.cos(phis) + b/2
    vec_bpA_y = b0_scalar*cp.sin(beta) + scalar_l*cp.sin(phis)
    vec_bpB_x = b0_scalar*cp.cos(beta) + scalar_l*cp.cos(phis) - b/2
    vec_bpB_y = b0_scalar*cp.sin(beta) + scalar_l*cp.sin(phis)

    bpA_size = vec_bpA_x*vec_bpA_x + vec_bpA_y*vec_bpA_y
    bpB_size = vec_bpB_x*vec_bpB_x + vec_bpB_y*vec_bpB_y

    step_bpA = cp.clip(RA-cp.sqrt(bpA_size), 10**(-500), )/(RA-cp.sqrt(bpA_size))
    step_bpB = cp.clip(RB-cp.sqrt(bpB_size), 10**(-500), )/(RB-cp.sqrt(bpB_size))
    print(TT(RA, vec_bpA_x, vec_bpA_y), TT(RB, vec_bpB_x, vec_bpB_y))

    return step_bpA*step_bpB*(SIGMA*KP*(TT(RA, vec_bpA_x, vec_bpA_y) + TT(RB, vec_bpB_x, vec_bpB_y))/(2*(T0+scalar_l)))

def N_ch(impact_param):
    b = -impact_param
    bp_scalar, alpha = np.meshgrid(np.arange(0.,1.59, 1.59/nbp), np.arange(0., np.pi*2, 2*np.pi/nalpha))
    bpx = bp_scalar*np.cos(alpha)
    bpy = bp_scalar*np.sin(alpha)
    result_dist = np.sum(C_CMS*(2/3)*KP*(TT_host(RA, bpx+b/2, bpy)+TT_host(RA, bpx-b/2, bpy)))
    return result_dist*(1.59/nbp)*(2*np.pi/nalpha)

def TT_host(R, bx, by):
    square_b = bx*bx+by*by
    stepfunction = np.heaviside(R-np.sqrt(square_b),0)
    squareroot = np.sqrt((R*R-square_b)*stepfunction)
    return (3*squareroot)/(2*np.pi*R*R*R)

def TT(R, bx, by):  # gpu
    square_b = bx*bx+by*by
    # step = cp.clip(RA-cp.sqrt(square_b), 0, )/(RA-cp.sqrt(square_b))
    return (3*cp.sqrt(cp.clip(R*R-square_b, 10**(-500), )))/(2*cp.pi*R*R*R)

    
def Pjet(b, scalar_b0, beta):   #gpu
    b0x = scalar_b0*cp.cos(beta)
    b0y = scalar_b0*cp.sin(beta)
    return TT(RA, b0x+b/2, b0y)*TT(RB, b0x-b/2, b0y)

def Pjet_host(b, scalar_b0, beta):
    b0x = scalar_b0*np.cos(beta)
    b0y = scalar_b0*np.sin(beta)
    return TT_host(RA, b0x+b/2, b0y)*TT_host(RB, b0x-b/2, b0y)


if __name__ == '__main__':
    impact_param = np.arange(0.,1.59, db)
    # Nk_result = Nk_b(impact_param)

    Nk_result = np.zeros(int(1.59/db)+1)
    Nch_result = np.zeros(int(1.59/db)+1)
    # Nk_result = np.zeros(nb)
    # print(int(1.59/db)+1)

    # pool = multiprocessing.Pool(processes=16)
    # Nk_result = pool.starmap(Nk_b, zip(impact_param))
    print("impact parameter \t Nk \t Nch")
    for i in range(int(1.59/db)+1):
        Nch_result[i+0] = N_ch(impact_param[i+0])
        Nk_result[i+0] = Nk_b(impact_param[i+0])
        print(impact_param[i+0], Nk_result[i+0], Nch_result[i+0])




    # Graphs
    np.savetxt('Nk_result.csv', Nk_result, delimiter=',')
    np.savetxt('Nch_result.csv', Nch_result, delimiter=',')

    fig = plt.figure()
    ax = plt.axes()
    fig.set_size_inches(35, 16.534, forward=True)


    plt.plot(impact_param, Nk_result, linewidth=7, label = r'$N_k$')
    plt.plot(impact_param, Nch_result, linewidth=7, label = r'$N_ch$')

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