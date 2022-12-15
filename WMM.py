import math
from scipy.special import lpmv


def lp(phi, n, m):
  if m == 0:
    return ((-1)**m) * lpmv(m, n, math.sin(phi))
  return math.sqrt(2 * math.factorial(n - m) / math.factorial(n + m)) * (
    (-1)**m) * lpmv(m, n, math.sin(phi))  # полином Лежандра


def dlp(phi, n, m):
  return (n + 1) * math.tan(phi) * lp(phi, n, m) - math.sqrt(
    (n + 1)**2 - m**2) * lp(phi, n + 1, m) / math.cos(phi) # производная полинома Лежандра

lambd = 29.8344626  # долгота, "+" соответствует E
phi = 59.8798694  # широта, "+" соответствует N
t = 2022 # время в годах (не обязательно целое)
h_msl = 0  # высота в метрах над уровнем моря

phi = math.pi * phi / 180
lambd = math.pi * lambd / 180


A = 6378137
f = 1/298.257223563
e = math.sqrt(f*(2-f))
R_c = A / math.sqrt(1-e**2*(math.sin(phi)**2))
t_0 = 2020

p = (R_c+h_msl)*math.cos(phi)
z = (R_c*(1-e**2)+h_msl)*math.sin(phi)
r = math.sqrt(p**2+z**2)
phi_ = math.asin(z/r) # параметры Земли как эллипсоида


h = [[0] * 13 for i in range(13)]  #g_n^m = g[n][m]
g = [[0] * 13 for i in range(13)]
hdot = [[0] * 13 for i in range(13)]
gdot = [[0] * 13 for i in range(13)]

with open('WMM.txt', 'r') as file:
  for line in file:
    l = line.strip().split()
    if len(l) > 1:
      n = int(l[0])
      m = int(l[1])
      # print(n, m)
      gdot[n][m] = float(l[4])
      hdot[n][m] = float(l[5])
      g[n][m] = float(l[2]) + (t - t_0) * gdot[n][m]
      h[n][m] = float(l[3]) + (t - t_0) * hdot[n][m]

B_x_ = 0
B_y_ = 0
B_z_ = 0

a = 6371200

for n in range(1, 13):
  for m in range(n+1):
    B_x_ -= (
        (a / r)**(n + 2)) * (g[n][m] * math.cos(m * lambd) +
                           h[n][m] * math.sin(m * lambd)) * dlp(phi_, n, m) 
    B_y_ += (((a / r)**(n + 2)) * m *
             (g[n][m] * math.sin(m * lambd) - h[n][m] * math.cos(m * lambd)) *
             lp(phi_, n, m)) / math.cos(phi_)
    B_z_ -= (n + 1) * (
      (a / r)**(n + 2)) * (g[n][m] * math.cos(m * lambd) +
                           h[n][m] * math.sin(m * lambd)) * lp(phi_, n, m) # Ряд Гаусса


X = B_x_*math.cos(phi_ - phi) - B_z_*math.sin(phi_-phi)
Y = B_y_
Z = B_x_*math.sin(phi_ - phi) + B_z_*math.cos(phi_-phi) # переход в систему координат, связанную с поверхностью эллипсоида

H = math.sqrt(X**2+Y**2)
F = math.sqrt(H**2+Z**2)
I = math.atan2(Z, H)
D = math.atan2(Y, X) # вычисление горизонтальной компоненты, модуля В, наклонения и склонения поля 

print(H, F, I*180/math.pi, sep='\n')
