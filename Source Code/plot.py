import matplotlib.pyplot as plt

# Data Single Tank Drainage
waktu1 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
tinggi_rk4 = [5, 4.802, 4.6078, 4.4174, 4.2308, 4.0481, 3.8693, 3.6943, 3.5232, 3.3559, 3.1926]
tinggi_analitik = [5, 4.802, 4.6078, 4.4174, 4.2308, 4.0481, 3.8693, 3.6943, 3.5232, 3.3559, 3.1926]

plt.figure()
plt.plot(waktu1, tinggi_rk4, marker='o', label='RK4')
plt.plot(waktu1, tinggi_analitik, marker='x', label='Analitik')
plt.title('Single Tank Drainage')
plt.xlabel('Waktu (s)')
plt.ylabel('Tinggi Air (m)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Data Sistem Dua Tangki
waktu2 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
tinggi_t1 = [1, 1.4901, 1.8569, 2.1394, 2.3639, 2.5462, 2.6968, 2.8221, 2.927, 3.0151, 3.0894]
tinggi_t2 = [0.5, 0.6842, 0.8416, 0.9857, 1.1221, 1.2534, 1.3812, 1.5064, 1.6295, 1.7508, 1.8707]
laju_aliran = [0.002978, 0.003633, 0.004054, 0.004351, 0.004573, 0.004747, 0.004885, 0.004996, 0.005087, 0.005163, 0.005227]

fig, ax1 = plt.subplots()
ax1.plot(waktu2, tinggi_t1, marker='o', label='Tangki 1')
ax1.plot(waktu2, tinggi_t2, marker='s', label='Tangki 2')
ax1.set_xlabel('Waktu (s)')
ax1.set_ylabel('Tinggi Air (m)')
ax1.grid(True)
ax1.legend(loc='upper left')

ax2 = ax1.twinx()
ax2.plot(waktu2, laju_aliran, marker='^', linestyle='--', color='red', label='Laju Aliran')
ax2.set_ylabel('Laju Aliran (m³/s)', color='red')
ax2.tick_params(axis='y', labelcolor='red')
ax2.legend(loc='upper right')

plt.title('Sistem Dua Tangki')
plt.tight_layout()
plt.show()