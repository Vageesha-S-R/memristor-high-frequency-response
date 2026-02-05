# High Frequency Response of Non-Volatile Memristors

This project reproduces and studies the dynamic behavior of TaOx-based non-volatile memristors following:

**Ioannis Messaris – _High Frequency Response of Non-Volatile Memristors_**

The device is modeled using a physics-based state equation and simulated under different excitation signals to analyze:

• Frequency response  
• Fading memory behavior  
• DC offset switching dynamics  

---

## Experiments

### 1. Sine Wave – Frequency Sweep (Python)
State evolution for multiple frequencies:
30 kHz, 150 kHz, 500 kHz, 100 MHz

→ Demonstrates transition from full switching to high-frequency averaging.

### 2. Fading Memory (Python)
Different initial conditions converge to the same steady state.

→ Shows non-volatile fading memory property.

### 3. Square Wave with DC Offset (MATLAB)
State convergence for different DC biases.

→ Demonstrates programmable state levels.

---

## How to Run

### Python simulation
```bash
pip install -r requirements.txt
python python/sine_frequency_sweep.py
python python/fading_memory.py
```
### MATLAB simulation
```bash
Open MATLAB and run:
matlab/square_wave_dc_offset.m
```

## Requirements

### Python
```bash
Python 3.9+
numpy
matplotlib
scipy
```
### MATLAB
```bash
MATLAB R2020+ (or compatible)
```
