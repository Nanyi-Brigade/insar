import xml.dom.minidom
import time
import numpy as np


def readEOFParams(xml_path):
    """读取精密轨道文件计算参数

    Args:
        xml_path (str): 精密轨道文件路径.

    Returns:
        ndarray: 参数.
    """
    dom = xml.dom.minidom.parse(xml_path)  # 打开XML文件
    collection = dom.documentElement  # 获取元素对象
    object_list = collection.getElementsByTagName("OSV")  # 获取标签名信息
    # print(object_list)
    X_LOC = []
    Y_LOC = []
    Z_LOC = []
    VX_VEL = []
    VY_VEL = []
    VZ_VEL = []
    UTC = []
    for obj in object_list:
        x = obj.getElementsByTagName('X')[0].childNodes[0].data
        y = obj.getElementsByTagName('Y')[0].childNodes[0].data
        z = obj.getElementsByTagName('Z')[0].childNodes[0].data
        vx = obj.getElementsByTagName('VX')[0].childNodes[0].data
        vy = obj.getElementsByTagName('VY')[0].childNodes[0].data
        vz = obj.getElementsByTagName('VZ')[0].childNodes[0].data
        utc = obj.getElementsByTagName('UTC')[0].childNodes[0].data
        X_LOC.append(float(x))
        Y_LOC.append(float(y))
        Z_LOC.append(float(z))
        VX_VEL.append(float(vx))
        VY_VEL.append(float(vy))
        VZ_VEL.append(float(vz))
        UTC.append(utc)
    UTC_TIME = []
    for i in UTC:
        UTC_data = time.mktime(time.strptime(i, 'UTC=%Y-%m-%dT%H:%M:%S.000000'))
        UTC_TIME.append(UTC_data)
    UTC_TIME = np.array(UTC_TIME)
    x_parameters = soluParameters(X_LOC, UTC_TIME)
    y_parameters = soluParameters(Y_LOC, UTC_TIME)
    z_parameters = soluParameters(Z_LOC, UTC_TIME)
    vx_parameters = soluParameters(VX_VEL, UTC_TIME)
    vy_parameters = soluParameters(VY_VEL, UTC_TIME)
    vz_parameters = soluParameters(VZ_VEL, UTC_TIME)
    all_parameters = np.array([x_parameters, y_parameters, z_parameters,
                               vx_parameters, vy_parameters, vz_parameters])
    return all_parameters


def soluParameters(local_velocity , time, n=3):
    """基于多项式拟合求解SAR轨道参数
       n次拟合曲线：y=a0+a1x+a2x^2+...+anx^n

    Args:
        local_velocity (list): 速度.
        time (list): 时间.
        n (int, optional): 拟合次数. 默认为 3.

    Returns:
        ndarray: 拟合参数.
    """
    x_array = time
    y_array = np.array(local_velocity)
    m = len(x_array)  # 方程个数
    A = np.ones(m).reshape((m, 1))
    for i in range(n):
        A = np.hstack([A, (x_array**(i + 1)).reshape((m, 1))])
    from numpy.linalg import solve 
    parameters = solve(np.dot(A.T, A), np.dot(A.T, y_array.T))
    return parameters