import math

def sclr(a,b):
  return sum(a[i]*b[i]  for i in range(3))

R_E = 6371.2 * 1e3  # радиус земли
mu = 7.70 * 10**22  # магнитный момент земли
mu_0 = 1.2566 * 10**-6  # магнитная проницаемость вакуума
pole_lat = 80.7  # координаты северного магнитного полюса
pole_long = 72.7

pole_lat *= math.pi / 180
pole_long *= math.pi / 180

coords = [59.8798694,-29.8344626]  # Координаты исследуемой точки на поверхности земли
coords[0] *= math.pi / 180
coords[1] *= math.pi / 180
height = 0 # высота над поверхностью земли (в метрах)
M = [
  mu * math.cos(pole_lat) * math.cos(pole_long),
  mu * math.cos(pole_lat) * math.sin(pole_long), mu * math.sin(pole_lat)
]  # вектор магнитного момента

R = [
  (R_E + height) * math.cos(coords[0]) * math.cos(coords[1]),
  (R_E + height) * math.cos(coords[0]) * math.sin(coords[1]), (R_E + height) * math.sin(coords[0])
]  # радиус - вектор

M_R = sclr(M, R)

B = [mu_0 / (4 * math.pi) * ((3 * M_R * R[i]) / ((R_E + height)**5) - M[i] / (R_E + height)**3) for i in range(3)]  # поле в исследуемой точке

B_abs = math.sqrt(sum([B[i]**2 for i in range(3)]))  # модуль вектора В

Sc = sclr(B, R)  # (B, R)

temp = math.acos(Sc / (B_abs * (R_E + height))) 
temp *= 180 / math.pi
D = 90 - temp  # наклонение

B_H = 1e6 * math.sqrt(B_abs**2 - Sc**2 / (R_E**2))
B_abs *= 1e6 # модуль и горизонтальная компонента

print('|B| = %.3f μT' % B_abs,
      'B_h = %.3f μT' % B_H,
      'Incl = %.3f °' % D,
      sep='\n')
