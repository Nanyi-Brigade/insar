import numpy as np
import matplotlib.pyplot as plt

# R根据像素坐标求得
#参考椭球的长半轴，短半轴：a，b；地面目标大地高：h
#卫星的位置及速度（Xs,Ys,Zs），（Vx，Vy，Vz），电磁波波长：bc

def test(h , R, bc , Xs , Ys , Zs ,Vx ,Vy ,Vz , result):

    a = 6378137 
    b = 6356752.314 
    
    #初始化Xp,Yp,Zp
    Xp = result[0][0]
    Yp = result[1][0]
    Zp = result[2][0]
 
    #求逆矩阵，用的matlab，距离-多普勒方程
    ni = np.array([[-(Vy*Zp*a**2 + 2*Vy*Zp*a*h - Vz*Yp*b**2 + Vy*Zp*h**2)/(2*(Vy*Xs*Zp*a**2 - Vy*Xp*Zp*a**2 - Vz*Xs*Yp*b**2 + Vz*Xp*Ys*b**2 - Vx*Ys*Zp*a**2 + Vx*Yp*Zp*a**2 - Vy*Xp*Zs*b**2 + Vy*Xp*Zp*b**2 + Vx*Yp*Zs*b**2 - Vx*Yp*Zp*b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*a*h - 2*Vy*Xp*Zp*a*h - 2*Vx*Ys*Zp*a*h + 2*Vx*Yp*Zp*a*h)),\
  -(Ys*Zp*a**2 - Yp*Zp*a**2 - Yp*Zs*b**2 + Yp*Zp*b**2 + Ys*Zp*h**2 - Yp*Zp*h**2 + 2*Ys*Zp*a*h - 2*Yp*Zp*a*h)/(Vy*Xs*Zp*a**2 - Vy*Xp*Zp*a**2 - Vz*Xs*Yp*b**2 + Vz*Xp*Ys*b**2 - Vx*Ys*Zp*a**2 + Vx*Yp*Zp*a**2 - Vy*Xp*Zs*b**2 + Vy*Xp*Zp*b**2 + Vx*Yp*Zs*b**2 - Vx*Yp*Zp*b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*a*h - 2*Vy*Xp*Zp*a*h - 2*Vx*Ys*Zp*a*h + 2*Vx*Yp*Zp*a*h),\
  (b**2*(a + h)**2*(Vz*Ys - Vz*Yp - Vy*Zs + Vy*Zp))/(2*(Vy*Xs*Zp*a**2 - Vy*Xp*Zp*a**2 - Vz*Xs*Yp*b**2 + Vz*Xp*Ys*b**2 - Vx*Ys*Zp*a**2 + Vx*Yp*Zp*a**2 - Vy*Xp*Zs*b**2 + Vy*Xp*Zp*b**2 + Vx*Yp*Zs*b**2 - Vx*Yp*Zp*b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*a*h - 2*Vy*Xp*Zp*a*h - 2*Vx*Ys*Zp*a*h + 2*Vx*Yp*Zp*a*h))],\
 [(Vx*Zp*a**2 + 2*Vx*Zp*a*h - Vz*Xp*b**2 + Vx*Zp*h**2)/(2*(Vy*Xs*Zp*a**2 - Vy*Xp*Zp*a**2 - Vz*Xs*Yp*b**2 + Vz*Xp*Ys*b**2 - Vx*Ys*Zp*a**2 + Vx*Yp*Zp*a**2 - Vy*Xp*Zs*b**2 + Vy*Xp*Zp*b**2 + Vx*Yp*Zs*b**2 - Vx*Yp*Zp*b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*a*h - 2*Vy*Xp*Zp*a*h - 2*Vx*Ys*Zp*a*h + 2*Vx*Yp*Zp*a*h)),\
 (Xs*Zp*a**2 - Xp*Zp*a**2 - Xp*Zs*b**2 + Xp*Zp*b**2 + Xs*Zp*h**2 - Xp*Zp*h**2 + 2*Xs*Zp*a*h - 2*Xp*Zp*a*h)/(Vy*Xs*Zp*a**2 - Vy*Xp*Zp*a**2 - Vz*Xs*Yp*b**2 + Vz*Xp*Ys*b**2 - Vx*Ys*Zp*a**2 + Vx*Yp*Zp*a**2 - Vy*Xp*Zs*b**2 + Vy*Xp*Zp*b**2 + Vx*Yp*Zs*b**2 - Vx*Yp*Zp*b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*a*h - 2*Vy*Xp*Zp*a*h - 2*Vx*Ys*Zp*a*h + 2*Vx*Yp*Zp*a*h),\
 -(b**2*(a + h)**2*(Vz*Xs - Vz*Xp - Vx*Zs + Vx*Zp))/(2*(Vy*Xs*Zp*a**2 - Vy*Xp*Zp*a**2 - Vz*Xs*Yp*b**2 + Vz*Xp*Ys*b**2 - Vx*Ys*Zp*a**2 + Vx*Yp*Zp*a**2 - Vy*Xp*Zs*b**2 + Vy*Xp*Zp*b**2 + Vx*Yp*Zs*b**2 - Vx*Yp*Zp*b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*a*h - 2*Vy*Xp*Zp*a*h - 2*Vx*Ys*Zp*a*h + 2*Vx*Yp*Zp*a*h))],\
 [(b**2*(Vy*Xp - Vx*Yp))/(2*(Vy*Xs*Zp*a**2 - Vy*Xp*Zp*a**2 - Vz*Xs*Yp*b**2 + Vz*Xp*Ys*b**2 - Vx*Ys*Zp*a**2 + Vx*Yp*Zp*a**2 - Vy*Xp*Zs*b**2 + Vy*Xp*Zp*b**2 + Vx*Yp*Zs*b**2 - Vx*Yp*Zp*b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*a*h - 2*Vy*Xp*Zp*a*h - 2*Vx*Ys*Zp*a*h + 2*Vx*Yp*Zp*a*h)),\
 -(b**2*(Xs*Yp - Xp*Ys))/(Vy*Xs*Zp*a**2 - Vy*Xp*Zp*a**2 - Vz*Xs*Yp*b**2 + Vz*Xp*Ys*b**2 - Vx*Ys*Zp*a**2 + Vx*Yp*Zp*a**2 - Vy*Xp*Zs*b**2 + Vy*Xp*Zp*b**2 + Vx*Yp*Zs*b**2 - Vx*Yp*Zp*b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*a*h - 2*Vy*Xp*Zp*a*h - 2*Vx*Ys*Zp*a*h + 2*Vx*Yp*Zp*a*h),\
 (b**2*(a + h)**2*(Vy*Xs - Vy*Xp - Vx*Ys + Vx*Yp))/(2*(Vy*Xs*Zp*a**2 - Vy*Xp*Zp*a**2 - Vz*Xs*Yp*b**2 + Vz*Xp*Ys*b**2 - Vx*Ys*Zp*a**2 + Vx*Yp*Zp*a**2 - Vy*Xp*Zs*b**2 + Vy*Xp*Zp*b**2 + Vx*Yp*Zs*b**2 - Vx*Yp*Zp*b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*a*h - 2*Vy*Xp*Zp*a*h - 2*Vx*Ys*Zp*a*h + 2*Vx*Yp*Zp*a*h))]])
   
    
    f1 = (Xp - Xs) * (Xp - Xs) + (Yp - Ys) * (Yp - Ys) + (Zp - Zs) * (Zp - Zs) - R * R
    f2 = (Xp - Xs) * Vx + (Yp - Ys) * Vy + (Zp - Zs) * Vz 
    f3 = (Xp * Xp + Yp * Yp) / ((a + h) * (a + h)) + (Zp * Zp) / (b * b) - 1
    yuan_fc = np.array([f1 , f2 , f3]).reshape(3 , 1)
    
    #牛顿迭代公式
    result = result - np.dot(ni , yuan_fc)
    
    return result

result = np.random.uniform(low = 0 , high = 10000000, size = (3 , 1))

cancha_x = []
cancha_y = []
cancha_z = []

for i in range(10000):
    result_1 = result
    result = test(100 , 591470.2436 , 0.031 , -789763.527 , 5987037.583 , 3306282.176 , 1076.539 , 3798.485 ,-6597.234 ,result)
    
    cancha_x.append(result_1[0][0] - result[0][0])
    cancha_y.append(result_1[1][0] - result[1][0])
    cancha_z.append(result_1[2][0] - result[2][0])
   
    if abs(result_1[2][0] - result[2][0]) <= 0.0001:
        print(result)
        break
   
x = list(range(len(cancha_x)))
#print(x)
plt.subplot(131)
plt.plot(x, cancha_x, color='black')
plt.subplot(132)
plt.plot(x, cancha_y, color='black')
plt.subplot(133)
plt.plot(x, cancha_z, color='black')
plt.show()