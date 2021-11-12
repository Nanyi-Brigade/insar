#通过DEM计算每一像素对应的大地高
import numpy as np
import matplotlib.pyplot as plt

#SRTM1对应3601*3601，SRTM3对应1201*1201
SAMPLES = 3601  # Change this to 1201 for SRTM3

def read_hgt(f_name, lat, lon):#文件位置，起始经纬度
    with open(f_name, 'rb') as hgt_data:
        elevations = np.fromfile(hgt_data, np.dtype('>i2'), SAMPLES * SAMPLES).reshape((SAMPLES, SAMPLES))
    lat_range = np.arange(lat, lat + 1 / SAMPLES + 1, 1 / (SAMPLES - 1))
    lon_range = np.arange(lon, lon + 1 + 1 / SAMPLES, 1 / (SAMPLES - 1))
    return lat_range, lon_range, elevations
    
def plot_d(ax):
    lat, lon, ele = read_hgt(r"C:\Users\张少\Desktop\N28E106.hgt" , 27 , 28)
    lim = np.arange(1, 2401, 100)
    ax.contourf(lon, lat, ele, lim, cmap='binary')
    return ax

fig ,ax = plt.subplots()
ax = plot_d(ax)
plt.show()