from scipy.optimize import curve_fit
import numpy as np
import cupy as cp
import scipy
import pprint as pp
import multiprocessing
import time

Nch_final = np.array([])
# multiple = 1


class function_gpu:
    '''fitting 하기 위해 함수가 실행된 횟수'''
    __count = 0
    
    '''N_k 계산에 사용할 코어 개수'''
    __Nk_core = 16

    '''상수 모음'''
    __RA = 0.8
    __RB = 0.8
    __KP = 367.       #ratio
    # KP = 21.      #RHIC
    __C_CMS = 0.269   #ratio
    __SIGMA = .14
    __T0 = 0.43       #fm/c
    # __T0 = 0.6        #RHIC
    __JETA = 0.20     #in eq(36)

    '''적분 bin 개수(Nk)'''
    __nphis = 50
    __nb0 = 300
    __nl_step = 1000
    '''적분 bin 개수(Nch)'''
    __nbp = 7500

    '''적분 범위(Nk)'''
    __phis_start = 0;   __phis_end = np.pi/2
    __b0_start = -0.8;    __b0_end = 0.8
    __l_start = __T0;   __l_end = __b0_end*4
    '''적분 범위(Nch)'''
    __bp_start = -0.8;  __bp_end = 0.8

    '''impact parameter를 numpy array 또는 list 형태로 대입하면 이에 대한 N_k 반환'''
    def Nk_b(self, impact_param):
        print("impact parameter \t Nk \t Nch")
        total_result = []
        for b in impact_param:
            dphis = (self.__phis_end - self.__phis_start)/self.__nphis
            phis = np.linspace(self.__phis_start, self.__phis_end, self.__nphis)      # radian
            pool = multiprocessing.Pool(processes = self.__Nk_core)
            b_array = np.zeros(len(phis)) - b
            result_dist = pool.starmap(self.Nk_phib, zip(phis, b_array))
            total_result.append(np.sum(result_dist)*dphis*(2/np.pi))
            print(b, total_result[-1])
        return total_result


    def Nk_phib(self, phis_host, b_scalar_host):   # scalar, vector
        b0y, l_scalar = cp.meshgrid(cp.linspace(self.__b0_start, self.__b0_end, self.__nb0), cp.linspace(self.__l_start, self.__l_end, self.__nl_step))
        b0x_arr = cp.linspace(self.__b0_start, self.__b0_end, self.__nb0)
        phis = cp.asarray(phis_host)
        b_scalar = cp.asarray(b_scalar_host)
        db0x = db0y = (self.__b0_end-self.__b0_start)/self.__nb0
        dl = (self.__l_end-self.__l_start)/self.__nl_step

        numerator = 0
        denominator = 0
        for b0x in b0x_arr:
            Nk_dist = cp.sum(self.__Nk_b0phib(b0x, b0y, l_scalar, phis, b_scalar)*dl, axis=0)
            Pjet_dist = self.__Pjet(b_scalar, b0x, b0y[0])
            same = cp.exp(-self.__JETA*Nk_dist)*Pjet_dist*db0x*db0y
            numerator += cp.sum(Nk_dist*same)
            denominator += cp.sum(same)
        numerator_host = cp.asnumpy(numerator) 
        denominator_host = cp.asnumpy(denominator)
        if denominator_host == 0:
            return 0
        else:
            return numerator_host/denominator_host


    def __Nk_b0phib(self, b0x, b0y, scalar_l, phis, b):   #scalar, scalar, scalar
        vec_bpA_x = b0x + scalar_l*cp.cos(phis) + b/2;    vec_bpA_y = b0y + scalar_l*cp.sin(phis)
        vec_bpB_x = b0x + scalar_l*cp.cos(phis) - b/2;    vec_bpB_y = b0y + scalar_l*cp.sin(phis)
        TA = self.__TTgpu(self.__RA, vec_bpA_x, vec_bpA_y);    TB = self.__TTgpu(self.__RB, vec_bpB_x, vec_bpB_y)
        step_bpA = cp.where(TA == 0, 0, 1);    step_bpB = cp.where(TB == 0, 0, 1)
        return step_bpA*step_bpB*self.__SIGMA*self.__KP*(TA + TB)/(2*(scalar_l))

    def N_ch(self, impact_param, multiple):
        # random = np.random.randint(0, len(xdata_dist))
        # xdata_dist[random] = xdata_dist[random] + impact_param
        # impact_param = xdata_dist
        # print(impact_param)

        # total_result = []
        total_result = np.array([])
        for b in impact_param:
            b = - b
            bpx, bpy = cp.meshgrid(cp.linspace(self.__bp_start,self.__bp_end, self.__nbp), cp.linspace(self.__bp_start,self.__bp_end, self.__nbp))
            # print(np.shape(bpx), np.shape(bpy))
            dbpx = (self.__bp_end-self.__bp_start)/self.__nbp;    dbpy = (self.__bp_end-self.__bp_start)/self.__nbp
            bpxA = bpx+b/2;    bpyA = bpy;    bpxB = bpx-b/2;    bpyB = bpy
            TA = self.__TTgpu(self.__RA, bpxA, bpyA);    TB = self.__TTgpu(self.__RB, bpxB, bpyB)
            step_bA = cp.where(TA==0, 0, 1);    step_bB = cp.where(TB==0, 0, 1)
            # total_result.append(cp.asnumpy(self.__C_CMS*(2/3)*self.__KP*cp.sum(step_bA*step_bB*(TA+TB))*dbpx*dbpy))
            total_result = np.append(total_result, cp.asnumpy(self.__C_CMS*(2/3)*self.__KP*cp.sum(step_bA*step_bB*(TA+TB))*dbpx*dbpy))
        self.__count = self.__count + 1

        # if self.__count % 50 == 0:
        #     print(f"{self.__count}회 \nimpact parameter")
        #     print(np.array(impact_param), "\nN_ch")
        #     print(total_result, "\n\n")

        global Nch_final
        Nch_final = total_result
        if multiple is None:
            ''' 그냥 결과 '''
            return total_result
        else:
            ''' fitting할 때 사용(Nch x 100)'''
            return total_result*multiple
    ''' fitting 하기 위해 wrapping '''
    def Nch_wrapper(self, *args):
        # time_start = time.time()
        impact_param = list(args[1::])
        # print(type(impact_param))
        result = self.N_ch(impact_param, multiple)
        # time_end = time.time()
        # print(f"Total end : {time_end-time_start:.3f} sec")
        return result
    
    '''multiprocessing'''
    '''def N_ch(self, impact_param):   
        b = - impact_param
        bpx, bpy = cp.meshgrid(cp.linspace(self.__bp_start,self.__bp_end, self.__nbp), cp.linspace(self.__bp_start,self.__bp_end, self.__nbp))
        # print(np.shape(bpx), np.shape(bpy))
        dbpx = (self.__bp_end-self.__bp_start)/self.__nbp;    dbpy = (self.__bp_end-self.__bp_start)/self.__nbp
        bpxA = bpx+b/2;    bpyA = bpy;    bpxB = bpx-b/2;    bpyB = bpy
        TA = self.__TTgpu(self.__RA, bpxA, bpyA);    TB = self.__TTgpu(self.__RB, bpxB, bpyB)
        step_bA = cp.where(TA==0, 0, 1);    step_bB = cp.where(TB==0, 0, 1)
        result_dist = self.__C_CMS*(2/3)*self.__KP*cp.sum(step_bA*step_bB*(TA+TB))*dbpx*dbpy

        return result_dist

    def Nch_wrapper(self, *args):
        time_start = time.time()
        impact_param = list(args[1::])
        pool = multiprocessing.Pool(processes = args[0])
        result = pool.starmap(self.N_ch, zip(impact_param))
        self.__count = self.__count + 1
        print(f"{self.__count}회 \n")
        print(impact_param, "\n")
        print(result, "\n\n")
        time_end = time.time()
        print(f"Total end : {time_end-time_start:.3f} sec")
        return result'''

    def __TTcpu(self, R, bx, by):
        square_b = bx*bx+by*by
        stepfunction = np.heaviside(R-np.sqrt(square_b),0)
        squareroot = np.sqrt((R*R-square_b)*stepfunction)
        return (3*squareroot)/(2*np.pi*R*R*R)

    def __TTgpu(self, R, bx, by):  # gpu
        square_b = bx*bx+by*by
        return (3*cp.sqrt(cp.clip(R*R-square_b, 0, )))/(2*cp.pi*R*R*R)

        
    def __Pjet(self, b, b0x, b0y):
        return self.__TTgpu(self.__RA, b0x+b/2, b0y)*self.__TTgpu(self.__RB, b0x-b/2, b0y)


class FITTING:
    '''N_k fitting에 사용할 코어 개수'''
    __Nch_fitting_core = 16

    ''' multiple : error를 크게 하기 위해 Nch x100 하려고 하는 것'''
    def __init__(self, Nch, initial, multiply):
        global multiple
        multiple = multiply
        self.Nch = Nch*multiply
        # print(Nch)
        # self.impact = np.linspace(0, 1.6, len(self.Nch))[::-1]
        if initial is None:
            self.impact = np.array([1.6    , 1.5315 , 1.50275, 1.48049, 1.46168, 1.44493, 1.42978, 1.41572, 1.4026 , 1.39019, 1.37849, 1.36727, 1.35652, 1.34612, 1.33613, 1.32645, 1.31705, 1.30792, 1.29905, 1.29036, 1.28186, 1.27361, 1.26546, 1.25751, 1.24966,
                                    1.24201, 1.23449, 1.22704, 1.21975, 1.21255, 1.20544, 1.19845, 1.19155, 1.18474, 1.17804, 1.17141, 1.16485, 1.15837, 1.15195, 1.14558, 1.1393 , 1.13306, 1.12688, 1.12078, 1.1147 , 1.10871, 1.10276, 1.0969 , 1.09107, 1.08524,
                                    1.07946, 1.07374, 1.06811, 1.06248, 1.05688, 1.05129, 1.0458 , 1.04034, 1.03488, 1.02949, 1.02411, 1.01873, 1.01344, 1.00816, 1.00291, 0.99768, 0.99247, 0.98729, 0.98219, 0.97701, 0.97195, 0.96687, 0.96181, 0.95682, 0.9518 ,
                                    0.94683, 0.94184, 0.93691, 0.93202, 0.92712, 0.92225, 0.91738, 0.91254, 0.90771, 0.90291, 0.89812, 0.89335, 0.88858, 0.88386, 0.87913, 0.87441, 0.86971, 0.86505, 0.86037, 0.85571, 0.85109, 0.84643, 0.84184, 0.83723, 0.83265,
                                    0.82805, 0.8235 , 0.81896, 0.81438, 0.80986, 0.80533, 0.8008 , 0.79633, 0.7918 , 0.78734, 0.78286, 0.77837, 0.77396, 0.76947, 0.76504, 0.76061, 0.75619, 0.75177, 0.74733, 0.74296, 0.73853, 0.73413, 0.72977, 0.72538, 0.72102,
                                    0.71664, 0.71227, 0.7079 , 0.70358, 0.6992 , 0.69489, 0.69052, 0.68618, 0.68186, 0.67752, 0.6732 , 0.66889, 0.66456, 0.66026, 0.65595, 0.65164, 0.64735, 0.64302, 0.63871, 0.63445, 0.63013, 0.62583, 0.62152, 0.61725, 0.61295,
                                    0.60863, 0.60438, 0.60007, 0.59578, 0.5915 , 0.58718, 0.58292, 0.57863, 0.5743 , 0.57003, 0.56575, 0.56142, 0.55715, 0.55284, 0.54852, 0.54425, 0.53991, 0.53559, 0.5313 , 0.52697, 0.52266, 0.51834, 0.514  , 0.5097 , 0.50535,
                                    0.501  , 0.49669, 0.49233, 0.48797, 0.48362, 0.47926, 0.47486, 0.47051, 0.4661 , 0.46173, 0.45732, 0.45292, 0.44849, 0.4441 , 0.43963, 0.43523, 0.43076, 0.42634, 0.42186, 0.41736, 0.41291, 0.40838, 0.4039 , 0.39937, 0.39481,
                                    0.39029, 0.38571, 0.38116, 0.37655, 0.37198, 0.36738, 0.36272, 0.35809, 0.35343, 0.34876, 0.34404, 0.33936, 0.3346 , 0.32986, 0.32506, 0.32029, 0.31545, 0.31064, 0.30576, 0.3009 , 0.29596, 0.29103, 0.28604, 0.28105, 0.27603,
                                    0.27098, 0.26586, 0.26074, 0.2556 , 0.2504 , 0.24515, 0.23991, 0.23455, 0.22919, 0.2238 , 0.21833, 0.21281, 0.20725, 0.20163, 0.19592, 0.19016, 0.18432, 0.17844, 0.17244, 0.16639, 0.16018, 0.15396, 0.14758, 0.14106, 0.13446,
                                    0.12768, 0.12075, 0.11362, 0.1063 , 0.09876, 0.09094, 0.0828 , 0.07427, 0.06526, 0.05562, 0.04512, 0.03322, 0.01874])
        else:
            self.impact = initial
        
        '''모든 boundary를 0~1.6으로 통일'''
        # bound_low = np.zeros(len(self.Nch))
        # bound_high = bound_low + 1.6
        # self.boundary = (bound_low, bound_high)

        '''각 boundary를 initial 값의 +-1로 통일'''
        bound_mid = initial
        bound_low = bound_mid - 0.000000000000001
        bound_high = bound_mid + 0.000000000000001
        self.boundary = (bound_low, bound_high)
        # print("boundary : \n", self.boundary)
    
    def fitting(self,):
        fitting_func = function_gpu()
        # popt, pcov = scipy.optimize.curve_fit(fitting_func.N_ch, xdata = self.impact, ydata = self.Nch, p0 = (self.impact))
        '''method = 'lm'(default), 'trf', 'dogbox' '''
        popt, pcov = scipy.optimize.curve_fit(fitting_func.Nch_wrapper, xdata = self.__Nch_fitting_core, ydata = self.Nch, bounds=self.boundary, p0 = (self.impact), method='dogbox')
        # popt, pcov = scipy.optimize.curve_fit(fitting_func.Nch_wrapper, xdata = self.__Nch_fitting_core, ydata = self.Nch, p0 = (self.impact))
        # print(popt)
        # return popt, np.sqrt(np.diag(pcov))
        return popt, Nch_final