from osgeo import gdal
import numpy as np
import math


def readSAR(tif_path):
    """ 读取sar图像的振幅和相位

    参数:
        tif_path (str): sar图像的路径.

    返回:
        tuple: 相位和振幅
    """
    dataset = gdal.Open(tif_path)
    im_width = dataset.RasterXSize  # 栅格矩阵的列数
    im_height = dataset.RasterYSize  # 栅格矩阵的行数
    im_bands = dataset.RasterCount  # 波段数
    im_data = dataset.ReadAsArray(0, 0, im_width, im_height)  # 获取数据
    im_geotrans = dataset.GetGeoTransform()  # 获取仿射矩阵信息
    a_real = im_data.real
    a_imag = im_data.imag
    amplitude = abs(im_data)  # 振幅
    phase = np.arctan2(a_imag, a_real)  # 相位
    return amplitude, phase