# oblicz kąt kierunkowy
alpha_1_2 = np.degrees(np.arctan2(p2_2000[0] - p1_2000[0], p2_2000[1] - p1_2000[1]))
alpha_2_3 = np.degrees(np.arctan2(p3_2000[0] - p2_2000[0], p3_2000[1] - p2_2000[1]))
alpha_3_4 = np.degrees(np.arctan2(p4_2000[0] - p3_2000[0], p4_2000[1] - p3_2000[1]))
alpha_4_1 = np.degrees(np.arctan2(p1_2000[0] - p4_2000[0], p1_2000[1] - p4_2000[1]))

# oblicz kąt kierunkowy odwrotny
alpha_2_1 = np.degrees(np.arctan2(p1_2000[0] - p2_2000[0], p1_2000[1] - p2_2000[1]))
alpha_3_2 = np.degrees(np.arctan2(p2_2000[0] - p3_2000[0], p2_2000[1] - p3_2000[1]))
alpha_4_3 = np.degrees(np.arctan2(p3_2000[0] - p4_2000[0], p3_2000[1] - p4_2000[1]))
alpha_1_4 = np.degrees(np.arctan2(p4_2000[0] - p1_2000[0], p4_2000[1] - p1_2000[1]))

# oblicz zbieżność południków γ
gamma = lambda_1 * np.sin(phi_1) + lambda_1**3 / 3 * np.sin(phi_1) * np.cos(phi_1)**2 * (1 + 3 * e2 + 2 * e2**2) + lambda_1**5 / 15 * np.sin(phi_1) * np.cos(phi_1)**4 * (2 - np.tan(phi_1)**2)

# oblicz redukcje kierunków δAB
delta_1_2 = (p2_2000[1] - p1_2000[1]) * (2 * p1_2000[0] + p2_2000[0]) / (6 * Rm1_2**2)
delta_2_1 = (p1_2000[1] - p2_2000[1]) * (2 * p2_2000[0] + p1_2000[0]) / (6 * Rm1_2**2)
delta_2_3 = (p3_2000[1] - p2_2000[1]) * (2 * p2_2000[0] + p3_2000[0]) / (6 * Rm2_3**2)
delta_3_2 = (p2_2000[1] - p3_2000[1]) * (2 * p3_2000[0] + p2_2000[0]) / (6 * Rm2_3**2)
delta_3_4 = (p4_2000[1] - p3_2000[1]) * (2 * p3_2000[0] + p4_2000[0]) / (6 * Rm3_4**2)
delta_4_3 = (p3_2000[1] - p4_2000[1]) * (2 * p4_2000[0] + p3_2000[0]) / (6 * Rm3_4**2)
delta_4_1 = (p1_2000[1] - p4_2000[1]) * (2 * p4_2000[0] + p1_2000[0]) / (6 * Rm4_1**2)
delta_1_4 = (p4_2000[1] - p1_2000[1]) * (2 * p1_2000[0] + p4_2000[0]) / (6 * Rm4_1**2)

# Oblicz azymut odcinka na powierzchni elipsoidy, jako:
# AAB = αAB + γA + δAB (12)
# i odpowiednio dla azymutu odwrotnego:
# ABA = αBA + γB + δBA (13)

A1_2 = alpha_1_2 + gamma + delta_1_2
A2_1 = alpha_2_1 + gamma + delta_2_1
A2_3 = alpha_2_3 + gamma + delta_2_3
A3_2 = alpha_3_2 + gamma + delta_3_2
A3_4 = alpha_3_4 + gamma + delta_3_4
A4_3 = alpha_4_3 + gamma + delta_4_3
A4_1 = alpha_4_1 + gamma + delta_4_1
A1_4 = alpha_1_4 + gamma + delta_1_4

# azymuty z zadania 3
A1_2_2 = geod.inv(lambda_1, phi_1, p2[0], p2[1])[0]
A2_3_2 = geod.inv(p2[0], p2[1], p3[0], p3[1])[0]
A3_4_2 = geod.inv(p3[0], p3[1], p4[0], p4[1])[0]
A4_1_2 = geod.inv(p4[0], p4[1], lambda_1, phi_1)[0]

print(f'Azymuty:\n')
print(f'Punkt 1-2:\nObliczony: {A1_2}\nZadanie 3: {A1_2_2}\n')
print(f'Punkt 2-3:\nObliczony: {A2_3}\nZadanie 3: {A2_3_2}\n')
print(f'Punkt 3-4:\nObliczony: {A3_4}\nZadanie 3: {A3_4_2}\n')
print(f'Punkt 4-1:\nObliczony: {A4_1}\nZadanie 3: {A4_1_2}\n')
