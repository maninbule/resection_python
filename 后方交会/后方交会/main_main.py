from numpy import * #引入Python的计算模块
#读入数据，并格式化保存到data_list
with open("坐标数据.txt",'r') as f:
    sou_data = f.read() #原始数据
    print(sou_data)
data = sou_data.split('\n')
data_list = []
for num in range(1,len(data)):
    list = data[num].split( )
    list2 = [float(it) for it in list]
    data_list.extend([list2])
#定义变量
pitch,roll,yaw = 0,0,0
ATAInv = None
delta = [0,0,0,0,0,0] #外方位元素所有改正数
scale = 40000  #摄影比例尺分母
f = 0.15324
num_ite = 0
x = [0,0,0,0]
y = [0,0,0,0]
A = mat(zeros((8,6)))
LL = mat(zeros((8,1)))
rotate = mat(zeros((3, 3)))#旋转矩阵
#从文件中录入像点坐标和地面点坐标
image_list = []
ground_list = []
for time3 in range(len(data_list)):
    image_list.append([data_list[time3][1]/1000, data_list[time3][2]/1000])
    ground_list.append([data_list[time3][3], data_list[time3][4],data_list[time3][5]])
Xs = (ground_list[0][0]+ground_list[1][0]+ground_list[2][0]+ground_list[3][0])/4
Ys = (ground_list[0][1]+ground_list[1][1]+ground_list[2][1]+ground_list[3][1])/4
Zs = scale*f
#打印外方位元素的初始值：
print('\n外方位元素的初始值为：\n')
print('%-15.f%-15.f%-15.f\n%-15.f%-15.f%-15.f\n'%(Xs,Ys,Zs,pitch,roll,yaw))
while(True):
    rotate[0,0] = cos(pitch)*cos(yaw)-sin(pitch)*sin(roll)*sin(yaw)
    rotate[0,1] = (-1.0) * cos(pitch) * sin(yaw) - sin(pitch) * sin(roll) * cos(yaw)
    rotate[0,2] = (-1.0) * sin(pitch) * cos(roll)
    rotate[1, 0] = cos(roll) * sin(yaw)
    rotate[1, 1] = cos(roll) * cos(yaw)
    rotate[1, 2] = (-1.0) * sin(roll)
    rotate[2, 0] = sin(pitch) * cos(yaw) + cos(pitch) * sin(roll) * sin(yaw)
    rotate[2, 1] = (-1.0) * sin(pitch) * sin(yaw) + cos(pitch) * sin(roll) * cos(yaw)
    rotate[2, 2] = cos(pitch) * cos(roll)
    for i in range(4):
        x[i] = (-1.0) * f * (rotate[0,0] * (ground_list[i][0] - Xs) + rotate[1,0] * (ground_list[i][1] - Ys)
                             + rotate[2,0] * (ground_list[i][2] - Zs)) / (rotate[0,2] * (ground_list[i][0] - Xs)+ rotate[1,2] * (ground_list[i][1] - Ys)
                              + rotate[2,2] * (ground_list[i][2] - Zs))
        y[i] = (-1.0) * f * (rotate[0,1] * (ground_list[i][0] - Xs) + rotate[1,1] * (ground_list[i][1] - Ys)
                             + rotate[2,1] * (ground_list[i][2] - Zs)) / (rotate[0,2] * (ground_list[i][0] - Xs)
                             + rotate[1,2] * (ground_list[i][1] - Ys) + rotate[2,2] * (ground_list[i][2] - Zs))
        LL[i * 2,0] = image_list[i][0] - x[i]
        LL[i * 2 + 1,0] = image_list[i][1] - y[i]
        A[i * 2,0] = (-1.0) * f / (Zs - ground_list[i][2])
        A[i * 2,1] = 0.0
        A[i * 2,2] = (-1.0) * x[i] / (Zs - ground_list[i][2])
        A[i * 2,3] = (-1.0) * f * (1 + (x[i] * x[i]) / (f * f))
        A[i * 2,4] = (-1.0) * x[i] * y[i] / f
        A[i * 2,5] = y[i]
        A[i * 2 + 1,0] = 0.0
        A[i * 2 + 1,1] = A[i * 2,0]
        A[i * 2 + 1,2] = (-1.0) * y[i] / (Zs - ground_list[i][2])
        A[i * 2 + 1,3] = A[i * 2,4]
        A[i * 2 + 1,4] = (-1.0) * f * (1 + (y[i] * y[i]) / (f * f))
        A[i * 2 + 1,5] = (-1.0) * x[i]

    AT = A.T #转置矩阵
    ATA = AT*A#矩阵相乘
    ATAInv = ATA.I#求逆矩阵
    ATAAT = ATAInv*AT
    delta = ATAAT*LL

    Xs += delta[0]
    Ys += delta[1]
    Zs += delta[2]
    pitch += delta[3]
    roll += delta[4]
    yaw += delta[5]

    num_ite+=1
    if (fabs(delta[3]) < 1e-6)and(fabs(delta[4]) < 1e-6)and(fabs(delta[5]) < 1e-6):
       break

print("迭代次数:%f\n"%num_ite)
rotate[0,0] = cos(pitch)*cos(yaw) - sin(pitch)*sin(roll)*sin(yaw)
rotate[0,1] = (-1.0)*cos(pitch)*sin(yaw) - sin(pitch)*sin(roll)*cos(yaw)
rotate[0,2] = (-1.0)*sin(pitch)*cos(roll)
rotate[1,0] = cos(roll)*sin(yaw)
rotate[1,1] = cos(roll)*cos(yaw)
rotate[1,2] = (-1.0)*sin(roll)
rotate[2,0] = sin(pitch)*cos(yaw) + cos(pitch)*sin(roll)*sin(yaw)
rotate[2,1] = (-1.0)*sin(pitch)*sin(yaw) + cos(pitch)*sin(roll)*cos(yaw)
rotate[2,2] = cos(pitch)*cos(roll)

AX = A*delta
v = AX - LL

vv = 0
for it in v:
    vv += it[0,0]*it[0,0]
m0 = sqrt(vv/2)
m = [0,0,0,0,0,0]
print(ATAInv)
for n in range(6):
    m[n] = m0 * sqrt(abs(ATAInv[n,n]))

print("\n像片的外方位元素为：\n")
print("Xs=%f  m=%f\n"%(Xs, m[0]))
print("Ys=%f  m=%f\n"%(Ys, m[1]))
print("Zs=%f   m=%f\n"% (Zs, m[2]))
print("pitch=%f  m=%f\n"% (pitch, m[3]))
print("roll=%f    m=%f\n"% (roll, m[4]))
print("yaw=%f    m=%f\n"%(yaw, m[5]))
print("\n旋转矩阵R为：\n")
print(rotate)
#写入数据到文件
with open('结果.txt','w') as f:
    f.writelines("\n像片的外方位元素为：\n")
    f.writelines("Xs=%f  m=%f\n" % (Xs, m[0]))
    f.writelines("Ys=%f  m=%f\n" % (Ys, m[1]))
    f.writelines("Zs=%f   m=%f\n" % (Zs, m[2]))
    f.writelines("pitch=%f  m=%f\n" % (pitch, m[3]))
    f.writelines("roll=%f    m=%f\n" % (roll, m[4]))
    f.writelines("yaw=%f    m=%f\n" % (yaw, m[5]))
    f.writelines("\n旋转矩阵R为：\n")
    for i in range(3):
        for j in range(3):
            f.write(str(rotate[i,j]).ljust(25,' '))
        f.write('\n')
