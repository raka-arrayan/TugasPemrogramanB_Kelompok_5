#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

class RungeKuttaSolver {
private:
    double h;  // Ukuran langkah (step size) untuk integrasi numerik
    
public:
    RungeKuttaSolver(double step_size) : h(step_size) {}
    
    /**
     * Melakukan satu langkah integrasi Runge-Kutta orde 4 untuk ODE tunggal
     * 
     * t = Waktu saat ini
     * y = Nilai fungsi saat ini
     * f = Pointer ke fungsi diferensial dy/dt = f(t, y)
     * Nilai y pada langkah berikutnya (t + h)
     * 
     * Proses:
     * 1. Hitung k1: gradien di titik awal
     * 2. Hitung k2: gradien di titik tengah menggunakan k1
     * 3. Hitung k3: gradien di titik tengah menggunakan k2
     * 4. Hitung k4: gradien di titik akhir menggunakan k3
     * 5. Kombinasi weighted average untuk mendapat nilai berikutnya
     */
    double rk4_step(double t, double y, double (*f)(double, double)) {
        // k1: Gradien di awal interval
        double k1 = h * f(t, y);
        
        // k2: Gradien di tengah interval dengan prediksi k1
        double k2 = h * f(t + h/2.0, y + k1/2.0);
        
        // k3: Gradien di tengah interval dengan prediksi k2 (lebih akurat)
        double k3 = h * f(t + h/2.0, y + k2/2.0);
        
        // k4: Gradien di akhir interval dengan prediksi k3
        double k4 = h * f(t + h, y + k3);
        
        // Weighted average dengan bobot 1:2:2:1 (Simpson's rule)
        return y + (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
    }
    
    /**
     * Menyelesaikan ODE tunggal dari waktu awal hingga akhir
     * 
     * t0 Waktu awal
     * tf Waktu akhir
     * y0 Kondisi awal y(t0)
     * f Fungsi diferensial dy/dt = f(t, y)
     * Vector berisi pasangan (waktu, nilai) sepanjang solusi
     */
    vector<pair<double, double>> solve_ode(double t0, double tf, double y0, 
                                          double (*f)(double, double)) {
        vector<pair<double, double>> solution;
        double t = t0;
        double y = y0;
        
        // Simpan kondisi awal
        solution.push_back({t, y});
        
        // Iterasi dari t0 sampai tf dengan langkah h
        while (t < tf) {
            y = rk4_step(t, y, f);  // Hitung nilai berikutnya
            t += h;                 // Maju satu langkah waktu
            solution.push_back({t, y});
        }
        
        return solution;
    }
    
    /**
     * RK4
     * Untuk sistem: dy1/dt = f1(t, y1, y2, ...), dy2/dt = f2(t, y1, y2, ...), ...
     * 
     * t0 Waktu awal
     * tf Waktu akhir
     * y0 Vector kondisi awal [y1(0), y2(0), ...]
     * functions Vector fungsi diferensial untuk setiap variabel
     * @return Matrix solusi: setiap baris = [waktu, y1, y2, ...]
     */
    vector<vector<double>> solve_system(double t0, double tf, 
                                       vector<double> y0,
                                       vector<double (*)(double, const vector<double>&)> functions) {
        vector<vector<double>> solution;
        double t = t0;
        vector<double> y = y0;
        int n = y0.size();  // Jumlah variabel dalam sistem
        
        // Simpan kondisi awal: [t0, y1(0), y2(0), ...]
        vector<double> point = {t};
        point.insert(point.end(), y.begin(), y.end());
        solution.push_back(point);
        
        // Iterasi waktu
        while (t < tf) {
            // Inisialisasi koefisien k untuk setiap variabel
            vector<double> k1(n), k2(n), k3(n), k4(n);
            vector<double> y_temp(n);
            
            // LANGKAH 1: Hitung k1 untuk semua variabel
            // k1[i] = h * f[i](t, y_current)
            for (int i = 0; i < n; i++) {
                k1[i] = h * functions[i](t, y);
            }
            
            // LANGKAH 2: Hitung k2 untuk semua variabel
            // k2[i] = h * f[i](t + h/2, y + k1/2)
            for (int i = 0; i < n; i++) {
                y_temp[i] = y[i] + k1[i]/2.0;  // Prediksi nilai tengah
            }
            for (int i = 0; i < n; i++) {
                k2[i] = h * functions[i](t + h/2.0, y_temp);
            }
            
            // LANGKAH 3: Hitung k3 untuk semua variabel
            // k3[i] = h * f[i](t + h/2, y + k2/2)
            for (int i = 0; i < n; i++) {
                y_temp[i] = y[i] + k2[i]/2.0;  // Prediksi nilai tengah (lebih akurat)
            }
            for (int i = 0; i < n; i++) {
                k3[i] = h * functions[i](t + h/2.0, y_temp);
            }
            
            // LANGKAH 4: Hitung k4 untuk semua variabel
            // k4[i] = h * f[i](t + h, y + k3)
            for (int i = 0; i < n; i++) {
                y_temp[i] = y[i] + k3[i];  // Prediksi nilai akhir
            }
            for (int i = 0; i < n; i++) {
                k4[i] = h * functions[i](t + h, y_temp);
            }
            
            // LANGKAH 5: Update semua variabel dengan weighted average
            // y_next[i] = y[i] + (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6
            for (int i = 0; i < n; i++) {
                y[i] = y[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
            }
            
            t += h;  // Maju satu langkah waktu
            
            // Simpan solusi: [t, y1, y2, ...]
            point = {t};
            point.insert(point.end(), y.begin(), y.end());
            solution.push_back(point);
        }
        
        return solution;
    }
    
    // Getter dan Setter untuk ukuran langkah
    void set_step_size(double new_h) { h = new_h; }
    double get_step_size() const { return h; }
};

/**
 * =====================================================================================
 * PARAMETER FISIK SISTEM TANGKI
 * =====================================================================================
 */

// Parameter tangki tunggal
const double A_tank = 3.0;      // Luas penampang tangki (m²)
const double A_hole = 0.01;     // Luas lubang drainase (m²)  
const double g = 9.81;          // Percepatan gravitasi (m/s²)
const double Cd = 0.6;          // Koefisien discharge (tanpa dimensi)

/**
 * Persamaan diferensial untuk drainase tangki tunggal
 * t Waktu (s) - tidak digunakan karena persamaan autonomous
 * h Tinggi air dalam tangki (m)
 * Laju perubahan tinggi air dh/dt (m/s)
 */
double tank_drainage_ode(double t, double h) {
    if (h <= 0) return 0.0;  // Tangki kosong, tidak ada aliran
    return -(Cd * A_hole * sqrt(2.0 * g * h)) / A_tank;
}

/**
 * Solusi analitis untuk drainase tangki (untuk validasi)
 * 
 * Derivasi:
 * dh/dt = -k*sqrt(h), dimana k = (Cd*A_hole*sqrt(2*g))/A_tank
 * Pemisahan variabel: dh/sqrt(h) = -k*dt
 * Integrasi: 2*sqrt(h) = -k*t + C
 * Kondisi awal h(0) = h0: C = 2*sqrt(h0)
 * Hasil: h(t) = (sqrt(h0) - k*t/2)²
 * 
 * t Waktu (s)
 * h0 Tinggi awal (m)
 * Tinggi air pada waktu t (m)
 */
double analytical_solution(double t, double h0) {
    double k = (Cd * A_hole * sqrt(2.0 * g)) / A_tank;
    double temp = sqrt(h0) - 0.5 * k * t;
    return (temp > 0) ? temp * temp : 0.0;  // Pastikan non-negatif
}

/**
 * =====================================================================================
 * SISTEM DUA TANGKI TERHUBUNG
 * =====================================================================================
 * 
 * [Inflow] -> [Tank 1] -> [Tank 2] -> [Outflow]
 * 
 * Persamaan:
 * Tank 1: dh1/dt = (Qin - Qout1) / A1
 * Tank 2: dh2/dt = (Qin2 - Qout2) / A2, dimana Qin2 = Qout1
 * 
 * Flow rates:
 * Qout1 = Cd * A1_hole * sqrt(2*g*h1)
 * Qout2 = Cd * A2_hole * sqrt(2*g*h2)
 */

// Parameter sistem dua tangki
const double A1 = 2.0, A2 = 1.5;           // Luas tangki 1 dan 2 (m²)
const double A1_hole = 0.008, A2_hole = 0.006; // Luas lubang tangki 1 dan 2 (m²)
const double Qin = 0.02;                    // Laju aliran masuk (m³/s)

/**
 * Persamaan diferensial untuk tinggi tangki 1
 * 
 * t Waktu (s) - tidak digunakan
 * y Vector state [h1, h2] (m)
 * @return dh1/dt (m/s)
 */
double two_tank_h1(double t, const vector<double>& y) {
    double h1 = max(0.0, y[0]);  // Pastikan tinggi non-negatif
    
    // Hitung laju aliran keluar tangki 1 menggunakan Torricelli's law
    double Q_out1 = (h1 > 0) ? Cd * A1_hole * sqrt(2.0 * g * h1) : 0.0;
    
    // Massa tangki 1: dh1/dt = (Qin - Qout1) / A1
    return (Qin - Q_out1) / A1;
}

/**
 * Persamaan diferensial untuk tinggi tangki 2
 * 
 * t Waktu (s) - tidak digunakan
 * y Vector state [h1, h2] (m)
 * @return dh2/dt (m/s)
 */
double two_tank_h2(double t, const vector<double>& y) {
    double h1 = max(0.0, y[0]);  // Pastikan tinggi non-negatif
    double h2 = max(0.0, y[1]);  // Pastikan tinggi non-negatif
    
    // Aliran masuk tangki 2 = aliran keluar tangki 1
    double Q_in2 = (h1 > 0) ? Cd * A1_hole * sqrt(2.0 * g * h1) : 0.0;
    
    // Aliran keluar tangki 2
    double Q_out2 = (h2 > 0) ? Cd * A2_hole * sqrt(2.0 * g * h2) : 0.0;
    
    // Massa tangki 2: dh2/dt = (Qin2 - Qout2) / A2
    return (Q_in2 - Q_out2) / A2;
}

/**
 * Menyimpan data solusi ODE tunggal ke file CSV
 * 
 * data Vector pasangan (waktu, tinggi)
 * filename Nama file output
 */
void saveToCSV(const vector<pair<double, double>>& data, const string& filename) {
    ofstream file(filename);
    
    // Header CSV
    file << "Time(s),Height(m)" << endl;
    
    // Data points
    for (const auto& point : data) {
        file << point.first << "," << point.second << endl;
    }
    
    file.close();
    cout << "Data berhasil disimpan ke: " << filename << endl;
}

/**
 * Menyimpan data solusi sistem ODE ke file CSV
 * 
 * data = Matrix solusi [waktu, h1, h2]
 * filename = Nama file output
 */
void saveSystemToCSV(const vector<vector<double>>& data, const string& filename) {
    ofstream file(filename);
    
    // Header CSV
    file << "Time(s),Tank1_Height(m),Tank2_Height(m)" << endl;
    
    // Data points
    for (const auto& point : data) {
        file << point[0] << "," << point[1] << "," << point[2] << endl;
    }
    
    file.close();
    cout << "Data sistem berhasil disimpan ke: " << filename << endl;
}

/**
 * =====================================================================================
 * MAIN FUNCTION - PROGRAM UTAMA
 * =====================================================================================
 */
int main() {
    cout << "======================================" << endl;
    cout << "RUNGE-KUTTA METHOD (4th ORDER)" << endl;
    cout << "Tank Drainage System Analysis" << endl;
    cout << "======================================" << endl;

    /*
     * ---------------------------------------------------------------------------------
     * PROBLEM 1: ANALISIS DRAINASE TANGKI TUNGGAL
     * ---------------------------------------------------------------------------------
     */
    cout << "\n--- PROBLEM 1: SINGLE TANK DRAINAGE ---" << endl;
    cout << "\nTank Parameters:" << endl;
    cout << "- Tank area: " << A_tank << " m²" << endl;
    cout << "- Drain hole area: " << A_hole << " m²" << endl;
    cout << "- Discharge coefficient: " << Cd << endl;
    cout << "- Initial water height: 5.0 m" << endl;

    // Parameter simulasi
    double h0 = 5.0;        // Tinggi awal air (m)
    double t0 = 0.0;        // Waktu awal (s)
    double tf = 100.0;      // Waktu akhir (s)
    double h_step = 0.5;    // Ukuran langkah (s)

    // Inisialisasi solver
    RungeKuttaSolver solver(h_step);
    
    cout << "\nSolving ODE: dh/dt = -(Cd * A_hole * sqrt(2*g*h)) / A_tank" << endl;
    
    // Selesaikan ODE menggunakan RK4
    auto solution = solver.solve_ode(t0, tf, h0, tank_drainage_ode);
    
    // Tampilkan hasil dan bandingkan dengan solusi analitis
    cout << "\nResults (every 10 seconds):" << endl;
    cout << setw(8) << "Time(s)" << setw(15) << "RK4 Height(m)" 
         << setw(18) << "Analytical(m)" << setw(15) << "Error(%)" << endl;
    cout << string(56, '-') << endl;
    
    // Tampilkan setiap 20 data points (≈ 10 detik dengan h=0.5)
    for (size_t i = 0; i < solution.size(); i += 20) {
        double t = solution[i].first;
        double h_rk4 = solution[i].second;
        double h_analytical = analytical_solution(t, h0);
        
        // Hitung error relatif (dalam persen)
        double error = (h_analytical > 1e-10) ? 
                      abs(h_rk4 - h_analytical) / h_analytical * 100.0 : 0.0;
        
        cout << setw(8) << fixed << setprecision(1) << t 
             << setw(15) << setprecision(4) << h_rk4
             << setw(18) << setprecision(4) << h_analytical
             << setw(15) << setprecision(2) << error << endl;
    }
    
    // Cari waktu ketika tangki kosong (tinggi < 0.01 m)
    double empty_time = tf;  // Default jika tidak kosong dalam waktu simulasi
    for (const auto& point : solution) {
        if (point.second < 0.01) {
            empty_time = point.first;
            break;
        }
    }
    
    cout << "\nTank emptying time: " << fixed << setprecision(1) 
         << empty_time << " seconds" << endl;

    /*
     * ---------------------------------------------------------------------------------
     * PROBLEM 2: ANALISIS SISTEM DUA TANGKI
     * ---------------------------------------------------------------------------------
     */
    cout << "\n--- PROBLEM 2: TWO-TANK SYSTEM ---" << endl;
    cout << "\nSystem Parameters:" << endl;
    cout << "- Tank 1 area: " << A1 << " m²" << endl;
    cout << "- Tank 2 area: " << A2 << " m²" << endl;
    cout << "- Inflow rate: " << Qin << " m³/s" << endl;
    cout << "- Initial heights: h1 = 1.0 m, h2 = 0.5 m" << endl;

    // Kondisi awal sistem: h1(0), h2(0)
    vector<double> y0_system = {1.0, 0.5};
    
    // Vector fungsi diferensial untuk sistem
    vector<double (*)(double, const vector<double>&)> funcs = {two_tank_h1, two_tank_h2};
    
    // Gunakan langkah lebih kecil untuk sistem yang lebih kompleks
    solver.set_step_size(0.2);
    tf = 50.0;  // Waktu simulasi lebih pendek untuk analisis sistem
    
    // Selesaikan sistem ODE
    auto system_solution = solver.solve_system(t0, tf, y0_system, funcs);
    
    // Tampilkan dinamika sistem
    cout << "\nSystem dynamics (every 5 seconds):" << endl;
    cout << setw(8) << "Time(s)" << setw(15) << "Tank1 Height(m)" 
         << setw(18) << "Tank2 Height(m)" << setw(15) << "Flow Rate(m³/s)" << endl;
    cout << string(66, '-') << endl;
    
    // Tampilkan setiap 25 data points (5 detik dengan h = 0.2)
    for (size_t i = 0; i < system_solution.size(); i += 25) {
        double t = system_solution[i][0];
        double h1 = system_solution[i][1];
        double h2 = system_solution[i][2];
        
        // Hitung laju aliran antar tangki
        double flow_rate = (h1 > 0) ? Cd * A1_hole * sqrt(2.0 * g * h1) : 0.0;
        
        cout << setw(8) << fixed << setprecision(1) << t 
             << setw(15) << setprecision(4) << h1
             << setw(18) << setprecision(4) << h2
             << setw(15) << setprecision(6) << flow_rate << endl;
    }

    /*
     * ---------------------------------------------------------------------------------
     * ANALISIS STEADY-STATE SISTEM DUA TANGKI
     * ---------------------------------------------------------------------------------
     */
    cout << "\nSteady-State Analysis:" << endl;
    double h1_final = system_solution.back()[1];
    double h2_final = system_solution.back()[2];
    
    // Hitung steady state
    // Kondisi steady: Qin = Qout1 = Qout2
    // Untuk tangki 2: Qin = Cd * A2_hole * sqrt(2*g*h2_steady)
    // Untuk tangki 1: Qin = Cd * A1_hole * sqrt(2*g*h1_steady) 
    double h2_theoretical = pow(Qin / (Cd * A2_hole), 2.0) / (2.0 * g);
    double h1_theoretical = pow(Qin / (Cd * A1_hole), 2.0) / (2.0 * g);
    
    cout << "Final heights: h1 = " << setprecision(4) << h1_final 
         << " m, h2 = " << h2_final << " m" << endl;
    cout << "Theoretical steady state: h1 = " << h1_theoretical 
         << " m, h2 = " << h2_theoretical << " m" << endl;

    /*
     * ---------------------------------------------------------------------------------
     * ANALISIS ERROR DENGAN BERBAGAI UKURAN STEP
     * ---------------------------------------------------------------------------------
     */
    cout << "\n--- ERROR ANALYSIS ---" << endl;
    cout << "Comparison of step sizes for single tank problem:" << endl;
    cout << setw(12) << "Step Size" << setw(15) << "Height at t=20s" 
         << setw(15) << "Error(%)" << setw(15) << "CPU Steps" << endl;
    cout << string(57, '-') << endl;
    
    // Test berbagai ukuran step
    vector<double> step_sizes = {2.0, 1.0, 0.5, 0.25, 0.1};
    double h_analytical_20 = analytical_solution(20.0, h0);
    
    for (double step : step_sizes) {
        RungeKuttaSolver temp_solver(step);
        auto temp_solution = temp_solver.solve_ode(0.0, 20.0, h0, tank_drainage_ode);
        double h_numerical = temp_solution.back().second;
        double error = abs(h_numerical - h_analytical_20) / h_analytical_20 * 100.0;
        int steps = temp_solution.size();
        
        cout << setw(12) << step << setw(15) << setprecision(4) << h_numerical
             << setw(15) << setprecision(2) << error 
             << setw(15) << steps << endl;
    }

    /*
     * ---------------------------------------------------------------------------------
     * EXPORT DATA UNTUK VISUALISASI
     * ---------------------------------------------------------------------------------
     */
    cout << "\n--- DATA EXPORT ---" << endl;
    saveToCSV(solution, "single_tank_drainage.csv");
    saveSystemToCSV(system_solution, "two_tank_system.csv");
    
    cout << "\nData files generated:" << endl;
    cout << "- single_tank_drainage.csv (untuk plotting drainase tangki tunggal)" << endl;
    cout << "- two_tank_system.csv (untuk plotting dinamika sistem dua tangki)" << endl;

    cout << "\n======================================" << endl;
    cout << "Runge-Kutta Analysis Complete!" << endl;
    cout << "Metode RK4 berhasil diterapkan untuk:" << endl;
    cout << "✓ Drainase tangki tunggal dengan validasi analitis" << endl;
    cout << "✓ Sistem dua tangki terhubung" << endl;
    cout << "✓ Analisis konvergensi dan akurasi" << endl;
    cout << "✓ Export data untuk visualisasi" << endl;
    cout << "======================================" << endl;

    return 0;
}