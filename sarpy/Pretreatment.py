#构建sar类
import cv2 as cv
from osgeo import gdal
import numpy as np
from numpy.linalg import solve 
import math
import xml.dom.minidom
import time

class Sar_cls():
    
    a = 6378137 
    b = 6356752.314 
        
    def __init__(self):
        
        self.result = np.random.uniform(low = 0 , high = 10000000, size = (3 , 1))
    
    #读取sar影响的振幅及相位    
    def read_amp_phase(self , file_road):
        dataset = gdal.Open(file_road)
        im_width = dataset.RasterXSize #栅格矩阵的列数
        im_height = dataset.RasterYSize #栅格矩阵的行数
        im_bands = dataset.RasterCount #波段数
        im_data = dataset.ReadAsArray(0,0,im_width,im_height)#获取数据
        im_geotrans = dataset.GetGeoTransform()#获取仿射矩阵信息
        a_real = im_data.real
        a_imag = im_data.imag
        amplitude = abs(im_data) #振幅
        phase = np.arctan2(a_imag,a_real)  #相位
        return amplitude , phase
    
    #读取精轨数据的坐标及速度
    def read_jgdata(self , file_road):    
        dom = xml.dom.minidom.parse(file_road)  # 打开XML文件
        collection = dom.documentElement  # 获取元素对象
        object_list = collection.getElementsByTagName("OSV")  # 获取标签名信息
        self.X_LOC = []
        self.Y_LOC = []
        self.Z_LOC = []
        self.VX_VEL = []
        self.VY_VEL = []
        self.VZ_VEL = []
        self.UTC = []      
        for obj in object_list:
            x = obj.getElementsByTagName('X')[0].childNodes[0].data
            y = obj.getElementsByTagName('Y')[0].childNodes[0].data
            z = obj.getElementsByTagName('Z')[0].childNodes[0].data
            vx = obj.getElementsByTagName('VX')[0].childNodes[0].data
            vy = obj.getElementsByTagName('VY')[0].childNodes[0].data
            vz = obj.getElementsByTagName('VZ')[0].childNodes[0].data
            utc = obj.getElementsByTagName('UTC')[0].childNodes[0].data
            self.X_LOC.append(float(x))
            self.Y_LOC.append(float(y))
            self.Z_LOC.append(float(z))
            self.VX_VEL.append(float(vx))
            self.VY_VEL.append(float(vy))
            self.VZ_VEL.append(float(vz))
            self.UTC.append(utc)
        self.local_velocity = np.array([self.X_LOC,self.Y_LOC,self.Z_LOC,self.VX_VEL,self.VY_VEL,self.VZ_VEL])
        self.UTC_TIME = []
        for i in self.UTC:
            self.UTC_data = time.mktime(time.strptime(i ,'UTC=%Y-%m-%dT%H:%M:%S.000000'))
            self.UTC_TIME.append(self.UTC_data)
        self.UTC_TIME = np.array(self.UTC_TIME)
        print("---读取成功---")
    
    #基于多项式拟合求解SAR轨道参数
    #n次拟合曲线：y=a0+a1x+a2x^2+...+anx^n
    def Solu_parameters(self):  
        all_parameters = []
        for y_array in self.local_velocity:
            x_array = self.UTC_TIME
            n = 3 #方程的次数（n<5）
            m = len(x_array)#方程个数
            A = np.ones(m).reshape((m , 1))
            for i in range(n):
                A = np.hstack([A , (x_array**(i+1)).reshape((m , 1))])
            parameters = solve(np.dot(A.T , A) , np.dot(A.T , y_array.T))
            all_parameters.append(list(parameters))
        return np.array(all_parameters)
    
    #参考椭球的长半轴，短半轴：a，b；地面目标大地高：h
    #卫星的位置及速度（Xs,Ys,Zs），（Vx，Vy，Vz），电磁波波长：bc
    def __test(self , h , R , bc , Xs , Ys , Zs ,Vx ,Vy ,Vz):
        #初始化Xp,Yp,Zp
        Xp = self.result[0][0]
        Yp = self.result[1][0]
        Zp = self.result[2][0]
        #求逆矩阵，用的matlab，距离-多普勒方程
        ni = np.array([[-(Vy*Zp*self.a**2 + 2*Vy*Zp*self.a*h - Vz*Yp*self.b**2 + Vy*Zp*h**2)/(2*(Vy*Xs*Zp*self.a**2 - Vy*Xp*Zp*self.a**2 - Vz*Xs*Yp*self.b**2 + Vz*Xp*Ys*self.b**2 - Vx*Ys*Zp*self.a**2 + Vx*Yp*Zp*self.a**2 - Vy*Xp*Zs*self.b**2 + Vy*Xp*Zp*self.b**2 + Vx*Yp*Zs*self.b**2 - Vx*Yp*Zp*self.b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*self.a*h - 2*Vy*Xp*Zp*self.a*h - 2*Vx*Ys*Zp*self.a*h + 2*Vx*Yp*Zp*self.a*h)),\
      -(Ys*Zp*self.a**2 - Yp*Zp*self.a**2 - Yp*Zs*self.b**2 + Yp*Zp*self.b**2 + Ys*Zp*h**2 - Yp*Zp*h**2 + 2*Ys*Zp*self.a*h - 2*Yp*Zp*self.a*h)/(Vy*Xs*Zp*self.a**2 - Vy*Xp*Zp*self.a**2 - Vz*Xs*Yp*self.b**2 + Vz*Xp*Ys*self.b**2 - Vx*Ys*Zp*self.a**2 + Vx*Yp*Zp*self.a**2 - Vy*Xp*Zs*self.b**2 + Vy*Xp*Zp*self.b**2 + Vx*Yp*Zs*self.b**2 - Vx*Yp*Zp*self.b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*self.a*h - 2*Vy*Xp*Zp*self.a*h - 2*Vx*Ys*Zp*self.a*h + 2*Vx*Yp*Zp*self.a*h),\
      (self.b**2*(self.a + h)**2*(Vz*Ys - Vz*Yp - Vy*Zs + Vy*Zp))/(2*(Vy*Xs*Zp*self.a**2 - Vy*Xp*Zp*self.a**2 - Vz*Xs*Yp*self.b**2 + Vz*Xp*Ys*self.b**2 - Vx*Ys*Zp*self.a**2 + Vx*Yp*Zp*self.a**2 - Vy*Xp*Zs*self.b**2 + Vy*Xp*Zp*self.b**2 + Vx*Yp*Zs*self.b**2 - Vx*Yp*Zp*self.b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*self.a*h - 2*Vy*Xp*Zp*self.a*h - 2*Vx*Ys*Zp*self.a*h + 2*Vx*Yp*Zp*self.a*h))],\
     [(Vx*Zp*self.a**2 + 2*Vx*Zp*self.a*h - Vz*Xp*self.b**2 + Vx*Zp*h**2)/(2*(Vy*Xs*Zp*self.a**2 - Vy*Xp*Zp*self.a**2 - Vz*Xs*Yp*self.b**2 + Vz*Xp*Ys*self.b**2 - Vx*Ys*Zp*self.a**2 + Vx*Yp*Zp*self.a**2 - Vy*Xp*Zs*self.b**2 + Vy*Xp*Zp*self.b**2 + Vx*Yp*Zs*self.b**2 - Vx*Yp*Zp*self.b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*self.a*h - 2*Vy*Xp*Zp*self.a*h - 2*Vx*Ys*Zp*self.a*h + 2*Vx*Yp*Zp*self.a*h)),\
     (Xs*Zp*self.a**2 - Xp*Zp*self.a**2 - Xp*Zs*self.b**2 + Xp*Zp*self.b**2 + Xs*Zp*h**2 - Xp*Zp*h**2 + 2*Xs*Zp*self.a*h - 2*Xp*Zp*self.a*h)/(Vy*Xs*Zp*self.a**2 - Vy*Xp*Zp*self.a**2 - Vz*Xs*Yp*self.b**2 + Vz*Xp*Ys*self.b**2 - Vx*Ys*Zp*self.a**2 + Vx*Yp*Zp*self.a**2 - Vy*Xp*Zs*self.b**2 + Vy*Xp*Zp*self.b**2 + Vx*Yp*Zs*self.b**2 - Vx*Yp*Zp*self.b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*self.a*h - 2*Vy*Xp*Zp*self.a*h - 2*Vx*Ys*Zp*self.a*h + 2*Vx*Yp*Zp*self.a*h),\
     -(self.b**2*(self.a + h)**2*(Vz*Xs - Vz*Xp - Vx*Zs + Vx*Zp))/(2*(Vy*Xs*Zp*self.a**2 - Vy*Xp*Zp*self.a**2 - Vz*Xs*Yp*self.b**2 + Vz*Xp*Ys*self.b**2 - Vx*Ys*Zp*self.a**2 + Vx*Yp*Zp*self.a**2 - Vy*Xp*Zs*self.b**2 + Vy*Xp*Zp*self.b**2 + Vx*Yp*Zs*self.b**2 - Vx*Yp*Zp*self.b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*self.a*h - 2*Vy*Xp*Zp*self.a*h - 2*Vx*Ys*Zp*self.a*h + 2*Vx*Yp*Zp*self.a*h))],\
     [(self.b**2*(Vy*Xp - Vx*Yp))/(2*(Vy*Xs*Zp*self.a**2 - Vy*Xp*Zp*self.a**2 - Vz*Xs*Yp*self.b**2 + Vz*Xp*Ys*self.b**2 - Vx*Ys*Zp*self.a**2 + Vx*Yp*Zp*self.a**2 - Vy*Xp*Zs*self.b**2 + Vy*Xp*Zp*self.b**2 + Vx*Yp*Zs*self.b**2 - Vx*Yp*Zp*self.b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*self.a*h - 2*Vy*Xp*Zp*self.a*h - 2*Vx*Ys*Zp*self.a*h + 2*Vx*Yp*Zp*self.a*h)),\
     -(self.b**2*(Xs*Yp - Xp*Ys))/(Vy*Xs*Zp*self.a**2 - Vy*Xp*Zp*self.a**2 - Vz*Xs*Yp*self.b**2 + Vz*Xp*Ys*self.b**2 - Vx*Ys*Zp*self.a**2 + Vx*Yp*Zp*self.a**2 - Vy*Xp*Zs*self.b**2 + Vy*Xp*Zp*self.b**2 + Vx*Yp*Zs*self.b**2 - Vx*Yp*Zp*self.b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*self.a*h - 2*Vy*Xp*Zp*self.a*h - 2*Vx*Ys*Zp*self.a*h + 2*Vx*Yp*Zp*self.a*h),\
     (self.b**2*(self.a + h)**2*(Vy*Xs - Vy*Xp - Vx*Ys + Vx*Yp))/(2*(Vy*Xs*Zp*self.a**2 - Vy*Xp*Zp*self.a**2 - Vz*Xs*Yp*self.b**2 + Vz*Xp*Ys*self.b**2 - Vx*Ys*Zp*self.a**2 + Vx*Yp*Zp*self.a**2 - Vy*Xp*Zs*self.b**2 + Vy*Xp*Zp*self.b**2 + Vx*Yp*Zs*self.b**2 - Vx*Yp*Zp*self.b**2 + Vy*Xs*Zp*h**2 - Vy*Xp*Zp*h**2 - Vx*Ys*Zp*h**2 + Vx*Yp*Zp*h**2 + 2*Vy*Xs*Zp*self.a*h - 2*Vy*Xp*Zp*self.a*h - 2*Vx*Ys*Zp*self.a*h + 2*Vx*Yp*Zp*self.a*h))]])

        f1 = (Xp - Xs) * (Xp - Xs) + (Yp - Ys) * (Yp - Ys) + (Zp - Zs) * (Zp - Zs) - R * R
        f2 = (Xp - Xs) * Vx + (Yp - Ys) * Vy + (Zp - Zs) * Vz 
        f3 = (Xp * Xp + Yp * Yp) / ((self.a + h) * (self.a + h)) + (Zp * Zp) / (self.b * self.b) - 1
        yuan_fc = np.array([f1 , f2 , f3]).reshape(3 , 1)
      
        self.result = self.result - np.dot(ni , yuan_fc)    #牛顿迭代公式
    
    def get_loc(self , h , R, bc , Xs , Ys , Zs ,Vx ,Vy ,Vz):
        
        for i in range(10000):
            result_1 = self.result
            self.__test(h , R, bc , Xs , Ys , Zs ,Vx ,Vy ,Vz)
            if abs(result_1[2][0] - self.result[2][0]) <= 0.0001:
                #print(self.result)
                return self.result
        