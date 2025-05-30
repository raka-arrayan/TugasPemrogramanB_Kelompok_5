# Kelompok 5

# PENERAPAN METODE RUNGE-KUTTA UNTUK ANALISIS SISTEM DRAINASE TANK

## Anggota Kelompok

| Nama | NPM |
|------|-----|
| Abednego Zebua| 2306161883 |
| Syahmi Hamdani| 2306220532 |
| Raka Arrayan Muttaqien| 2306161800 |
| Wilman Saragih Sitio| 2306161776 |
| Matthew Immanuel | 2306221024 |

----

## Deskripsi Proyek

### Metode Runge-Kutta Orde 4 (RK4)

Metode Runge-Kutta Orde 4 (RK4) adalah metode numerik yang umum digunakan untuk menyelesaikan persamaan diferensial biasa (Ordinary Differential Equations / ODE) dalam bentuk awal (initial value problem). Metode ini dikenal karena keseimbangan optimal antara akurasi tinggi dan efisiensi komputasi.

#### Rumus Iterasi RK4

Untuk sistem ODE satu variabel:

```
yᵢ₊₁ = yᵢ + (1/6)(k₁ + 2k₂ + 2k₃ + k₄)
```

Dengan:

```
k₁ = h ⋅ f(tᵢ, yᵢ)
k₂ = h ⋅ f(tᵢ + h/2, yᵢ + k₁/2)
k₃ = h ⋅ f(tᵢ + h/2, yᵢ + k₂/2)
k₄ = h ⋅ f(tᵢ + h, yᵢ + k₃)
```


- `h`: step size
- `f(t, y)`: fungsi turunan dari y terhadap t

#### Keunggulan RK4

- Error per langkah: O(h⁵)
- Error total: O(h⁴)
- Stabilitas: Baik untuk banyak kasus
- Evaluasi fungsi: 4 kali per langkah

---

## Tujuan Proyek

Proyek ini bertujuan untuk menerapkan metode numerik Runge-Kutta orde 4 (RK4) dalam menyelesaikan sistem ODE yang muncul pada permasalahan drainase tank. Dua skenario utama dianalisis:

1. Proses pengosongan pada sebuah tank silinder tunggal
2. Dinamika fluida dalam sistem dua tank yang saling terhubung dengan aliran masuk konstan

Implementasi dilakukan menggunakan bahasa C++

---

## Studi Kasus

### 1. Drainase Tank Tunggal

#### Persamaan ODE:

```
dh/dt = - (Cd ⋅ Ahole ⋅ √(2gh)) / Atank
```

#### Parameter:

- Cd = 0.6
- Ahole = 0.01 m²
- Atank = 3.0 m²
- h₀ = 5.0 m
- g = 9.81 m/s²

#### Disederhanakan menjadi:

```
dh/dt = -0.00894√h
```

#### Solusi Analitik:

```
h(t) = (5.0 - 0.00447t)²
```

---

### 2. Sistem Dua Tank yang Terhubung

#### Parameter:

- Tank 1:
  - A₁ = 2.0 m²
  - A₁₀ = 0.008 m²
- Tank 2:
  - A₂ = 1.5 m²
  - A₂₀ = 0.006 m²
- Debit masuk ke Tank 1:
  - Qin = 0.02 m³/s
- Kondisi awal:
  - h₁(0) = 1.0 m
  - h₂(0) = 0.5 m

#### Sistem ODE:

```
dh₁/dt = (Qin / A₁) - (Cd ⋅ A₁₀ ⋅ √(2gh₁)) / A₁

dh₂/dt = (Cd ⋅ A₁₀ ⋅ √(2gh₁)) / A₂ - (Cd ⋅ A₂₀ ⋅ √(2gh₂)) / A₂
```

---

## Tujuan Kasus Dua Tank

- Menentukan waktu mencapai steady-state
- Menganalisis dampak dinamika aliran antar tank
- Menunjukkan kapabilitas RK4 dalam sistem ODE nonlinier

---

## Hasil yang Diharapkan

- Validasi metode Runge-Kutta orde 4 (RK4) dengan solusi analitik pada model tangki tunggal

- Visualisasi dinamika aliran fluida dalam sistem dua tangki hingga mencapai kondisi setimbang

- Analisis sensitivitas ukuran langkah (step size) terhadap akurasi hasil simulasi dan efisiensi komputasi

- Penerapan metode numerik untuk simulasi sistem teknik dalam pemodelan dinamika fluida
****

<div align="center">
  <p>© 2025 Tugas Pemrograman B Kelompok 5</p>
</div>

