import xml.dom.minidom
dom = xml.dom.minidom.parse('C:/Users/张少/Desktop/S1A_OPER_AUX_POEORB_OPOD_20210309T055829_V20151022T225943_20151024T005943.eof')  # 打开XML文件
collection = dom.documentElement  # 获取元素对象
object_list = collection.getElementsByTagName("OSV")  # 获取标签名信息
#print(object_list)

X = []
Y = []
Z = []
VX = []
VY = []
VZ = []

for obj in object_list:
    x = obj.getElementsByTagName('X')[0].childNodes[0].data
    y = obj.getElementsByTagName('Y')[0].childNodes[0].data
    z = obj.getElementsByTagName('Z')[0].childNodes[0].data
    vx = obj.getElementsByTagName('VX')[0].childNodes[0].data
    vy = obj.getElementsByTagName('VY')[0].childNodes[0].data
    vz = obj.getElementsByTagName('VZ')[0].childNodes[0].data
    
    X.append(float(x))
    Y.append(float(y))
    Z.append(float(z))
    VX.append(float(vx))
    VY.append(float(vy))
    VZ.append(float(vz))

print(X,Y,Z,VX,VY,VZ)