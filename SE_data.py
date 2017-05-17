#!/usr/bin/python3

import math

# Solution of the equation by the golden section
tau = (math.sqrt(5.0) - 1.0) / 2.0 # ~0.618
precision = 1e-12
limit = 1e+12

def gold(f):
    a = 0.1
    b = limit
    x1, x2 = b - (b - a) * tau, a + (b - a) * tau
    y1, y2 = f(x1), f(x2)
    while math.fabs(b - a) > precision:
        if math.fabs(y1) >= math.fabs(y2):
            a = x1
            x1 = x2
            y1 = y2
            x2 = a + (b - a) * tau
            y2 = f(x2)
        else:
            b = x2
            x2 = x1
            y2 = y1
            x1 = b - (b - a) * tau
            y1 = f(x1)
    return (b + a) / 2.0

# Shock-expansion (SE) theory

t = 0.00508                                 # chort length of doublewerge airfoil [m]
c = 0.0508                                  # maximum thickness of doublewerge airfoil [m]
M = 1.86                                    # free stream Mach number
degrees = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0]  # anti-clockwise angle of attack [grad]
P = 101325.0                                # static pressure [Pa]
r = 1.4                                     # ratio of specific heats
R = 287.058                                 # specific gas constant in [J/(kg*K)]
T = 288.15                                  # free stream temperature [K]

# Function
_lambda = lambda _theta:(( M**2-1)**2-3.0*(1+((r-1)/2)*M**2)*(1+((r+1)/2)*M**2)*((math.tan(_theta))**2))
_chi = lambda _theta, _lambd: ((M**2-1)**3-9.0*(1+((r-1)/2)*(M**2))*(1+((r-1)/2)*(M**2)+((r+1)/4)*(M**4))*(math.tan(_theta)**2))/((_lambd)**3)
_beta = lambda _theta, _lambd, _chib: math.atan((M**2-1+2*(_lambd)*math.cos((4*math.pi*SWDv+math.acos(_chib))/3))/(3*(1+((r-1)/2)*(M**2))*math.tan(_theta)))
_P = lambda _betap: P*(1+((2.0*r)/(r+1))*((M*math.sin(_betap))**2-1))
_M = lambda _betam, _theta: math.sqrt(((M*math.sin(_betam)**2)+(2/(r-1)))/(((2.0*r)/(r-1))*M*math.sin(_betam)-1))/(math.sin(_betam-_theta))
_T = lambda _betat: T*(1+((2.0*r)/(r+1))*((M*math.sin(_betat))**2-1))*(((r-1)*(M*math.sin(_betat)**2)+2)/((r+1)*(M*math.sin(_betat)**2)))
_PT = lambda _M1, _M2: (1+((r-1)/2)*(_M1**2))/(1+((r-1)/2)*(_M2**2))
_fM = lambda _M: math.sqrt((r+1)/(r-1))*math.atan(math.sqrt(((r-1)/(r+1))*(_M**2-1)))-math.atan(math.sqrt(_M**2-1))

cs = math.sqrt(r*R*T)
l = 0.5*(math.sqrt(t**2+c**2))
it = 0
Cd_i = []
Cl_i = []
L_D_i = []

M2_ = [] # surface №1
M4_ = [] # surface №2
M5_ = [] # surface №3
M3_ = [] # surface №4

P2_ = [] # surface №1
P4_ = [] # surface №2
P5_ = [] # surface №3
P3_ = [] # surface №4

T2_ = [] # surface №1
T4_ = [] # surface №2
T5_ = [] # surface №3
T3_ = [] # surface №4
for step in range (len(degrees)):
    a = math.fabs(math.radians(degrees[it]))
    e = math.atan(t/c)
    Dv = _lambda((math.fabs(e)+math.fabs(a)))
    SWDv = 1 # weak shock case
    if Dv < 0:
        print('\n(λ < 0) detached shock exist or invalid data entry: Computation not possible !')
    else:
        Q2=math.fabs(e)-math.fabs(a)
        if Q2 > 0:
                Lambda = math.sqrt(_lambda(Q2))
                chi = _chi(Q2, Lambda)
                B2 = _beta(Q2, Lambda, chi)
                P2 = _P(B2)
                M2 = _M(B2, Q2)
                T2 = _T(B2)
        else:
                M2 = math.fabs(gold(lambda IM2: _fM(M) - _fM(IM2) + math.fabs(Q2)))
                P2 = P*_PT(M, M2)**(r/(r-1))
                T2 = T*_PT(M, M2)
    
    Q3 = math.fabs(e)+math.fabs(a)
    lambda1 = math.sqrt(_lambda(Q3))
    chi1 = _chi(Q3, lambda1)
    B3 = _beta(Q3, lambda1, chi1)
    P3 = _P(B3)
    M3 = _M(B3, Q3)
    T3 = _T(B3)
    Q4 = 2.0*math.fabs(e)
    M4 = math.fabs(gold(lambda IM4: _fM(M2) - _fM(IM4) + Q4))
    P4 = P2*_PT(M2, M4)**(r/(r-1))
    T4 = T2*_PT(M2, M4)
    Q5 = 2.0*math.fabs(e)
    M5 = math.fabs(gold(lambda IM5: _fM(M3) - _fM(IM5) + Q5))
    P5 = P3*_PT(M3, M5)**(r/(r-1))
    T5 = T3*_PT(M3, M5)

    L = l*((P5-P2)*math.cos(math.fabs(e)-math.fabs(a))+(P3-P4)*math.cos(math.fabs(e)+math.fabs(a)))
    D = l*((P5-P2)*math.sin(math.fabs(e)-math.fabs(a))+(P3-P4)*math.sin(math.fabs(e)+math.fabs(a)))
    if  (D == 0.0): RLD = 0.0
    else: RLD = L/D
    V = M*cs
    dp = 0.5*P/(R*T)*V**2
    Cl = L/(dp*c)
    Cd = D/(dp*c)
    
    M2_.append(M2)
    M4_.append(M4)
    M5_.append(M5)
    M3_.append(M3)

    P2_.append(P2)
    P4_.append(P4)
    P5_.append(P5)
    P3_.append(P3)

    T2_.append(T2)
    T4_.append(T4)
    T5_.append(T5)
    T3_.append(T3)
    
    Cd_i.append(Cd)
    Cl_i.append(Cl)
    L_D_i.append(RLD)
    it+=1

# Write data files
def ul(f, inf1, inf2):
	i = 0
	j = 0
	f.write(str(-c/2)+'\t')
	for _ in range (4*len(degrees)+1):
		if i < len(degrees):
			f.write(str(inf1[j])+'\t')
		elif i == len(degrees):
			f.write('\n' + '0\t')
			j = 0
		elif i < 2*len(degrees):
			f.write(str(inf1[j-1])+'\t')
		elif i == 2*len(degrees):
			f.write(str(inf1[j-1])+'\n' + '0\t')
			j = 0
		elif i < 3*len(degrees):
			f.write(str(inf2[j-1])+'\t')
		elif i == 3*len(degrees):
			f.write(str(inf2[j-1])+'\n' + str(c/2) +'\t')
			j = 0
		elif i < 4*len(degrees):
			f.write(str(inf2[j-1])+'\t')
		else:
			f.write(str(inf2[j-1]))
		j+=1
		i+=1
	return 0

Cd_ = open('Cd.txt', 'w')
Cl_ = open('Cl.txt', 'w')
L_D_ = open('L_D.txt', 'w')
Ma_upper = open('Ma_upper.txt', 'w')
Ma_lower = open('Ma_lower.txt', 'w')
P_upper = open('P_upper.txt', 'w')
P_lower = open('P_lower.txt', 'w')
T_upper = open('T_upper.txt', 'w')
T_lower = open('T_lower.txt', 'w')
P_m = open('P_m.txt', 'w')
Ma_m = open('Ma_m.txt', 'w')
T_m = open('T_m.txt', 'w')

Cd_.write('# grad	Cd\n')
Cl_.write('# grad	Cl\n')
L_D_.write('# grad	l/D\n')
Ma_upper.write('# x	2grad	4grad	 6grad	8grad	10grad 12grad\n')
Ma_lower.write('# x	2grad	4grad	 6grad	8grad	10grad 12grad\n')
P_upper.write('# x	2grad	4grad	 6grad	8grad	10grad 12grad\n')
P_lower.write('# x	2grad	4grad	 6grad	8grad	10grad 12grad\n')
T_upper.write('# x	2grad	4grad	 6grad	8grad	10grad 12grad\n')
T_lower.write('# x	2grad	4grad	 6grad	8grad	10grad 12grad\n')
P_m.write('# grad	surf1	surf2	 surtf3	surf4\n')
Ma_m.write('# grad	surf1	surf2	 surtf3	surf4\n')
T_m.write('# grad	surf1	surf2	 surtf3	surf4\n')

ul(Ma_upper, M2_, M4_)
ul(Ma_lower, M3_, M5_)
ul(P_upper, P2_, P4_)
ul(P_lower, P3_, P5_)
ul(T_upper, T2_, T4_)
ul(T_lower, T3_, T5_)

k = 0
for _ in range (len(degrees)):
	Cl_.write(str(degrees[k]) + '\t' + str(Cl_i[k]) + '\n')
	L_D_.write(str(degrees[k]) + '\t'+ str(L_D_i[k]) + '\n')
	Cd_.write(str(degrees[k]) + '\t' + str(Cd_i[k]) + '\n')
	
	P_m.write(str(degrees[k])+'\t'+str(P2_[k])+'\t'+str(P4_[k])+'\t'+str(P5_[k])+'\t'+str(P3_[k])+'\n')
	Ma_m.write(str(degrees[k])+'\t'+str(M2_[k])+'\t'+str(M4_[k])+'\t'+str(M5_[k])+'\t'+str(M3_[k])+'\n')
	T_m.write(str(degrees[k])+'\t'+str(T2_[k])+'\t'+str(T4_[k])+'\t'+str(T5_[k])+'\t'+str(T3_[k])+'\n')
	
	k+=1

# close files
Cd_.close()
Cl_.close()
L_D_.close()
Ma_upper.close()
Ma_lower.close()
P_upper.close()
P_lower.close()
T_upper.close()
T_lower.close()
P_m.close()
Ma_m.close()
T_m.close()
