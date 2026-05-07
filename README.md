# Cairo University – Faculty of Engineering

**Electronics and Electrical Communications**  
Mainstream – Second Year  
Fields and Wave Propagation - EECG252

---

# Assignment 1

Submitted for EECG252 Assignment  
Under the supervision of **Dr. Mohamed A. Nasr**

| Code | Sec | BN | Name |
|------|-----|----|------|
| 91240804 | 4 | 13 | مينا وائل تناغو فهمى سدراك |

---

## Contents

1. [Analytical Solution](#1-analytical-solution)
   - 1.1 [Hand Analysis](#11-hand-analysis)
     - 1.1.1 [Finding Electrical Potential Everywhere](#111-finding-electrical-potential-everywhere)
     - 1.1.2 [Finding Surface Charge Density on the Infinite Plane](#112-finding-surface-charge-density-on-the-infinite-plane)
   - 1.2 [MATLAB Results](#12-matlab-results)
2. [CST Solution](#2-cst-solution)
   - 2.1 [Setting Up](#21-setting-up)
   - 2.2 [Finding Electrical Potential Everywhere](#22-finding-electrical-potential-everywhere)
   - 2.3 [Finding Surface Charge Density on the Infinite Plane](#23-finding-surface-charge-density-on-the-infinite-plane)
3. [Comparison Between Results](#3-comparison-between-results)

---

## 1 Analytical Solution

> **Assuming Q = 2 nC, D = 5 m.**

### 1.1 Hand Analysis

#### 1.1.1 Finding Electrical Potential Everywhere

We will start the hand analysis by performing the method of images.

We know that the Electrical Field of a point charge (Q) at any point is:

```
        Q
E = ———————— r̂        ...(1.01)
     4πε₀r²
```

For the main point charge (Q), the Electrical Field is:

```
         Q
E₁ = ———————— r̂₁       ...(1.02)
      4πε₀r₁²
```

Where:

```
|r₁| = √( x² + y² + (z − D)² )      ...(1.03)
```

For the image point charge (−Q), the Electrical Field is:

```
         −Q
E₂ = ———————— r̂₂       ...(1.04)
      4πε₀r₂²
```

Where:

```
|r₂| = √( x² + y² + (z + D)² )      ...(1.05)
```

To get the total Electrical Field:

```
E_T = E₁ + E₂           ...(1.06)

        Q     ⎛  1      1  ⎞
E_T = ———— × ⎜ ——— − ——— ⎟ Ê    ...(1.07)
       4πε₀  ⎝  r₁²   r₂² ⎠
```

We know that:

```
φ = − ∫ E_T · dl         ...(1.08)

        Q     ⎛  1      1  ⎞
φ = − ∫ ———— ⎜ ——— − ——— ⎟ dr   ...(1.09)
       4πε₀  ⎝  r₁²   r₂² ⎠

      Q    ⎛  1     1  ⎞
φ = ———— ⎜ ——— − ——— ⎟          ...(1.10)
     4πε₀ ⎝  r₁    r₂ ⎠
```

By substituting, we get the potential everywhere:

```
      Q    ⎛          1                        1           ⎞
φ = ———— ⎜ ———————————————————— − ———————————————————— ⎟   ...(1.11)
     4πε₀ ⎝ √(x²+y²+(z−D)²)      √(x²+y²+(z+D)²)    ⎠
```

**Values of φ at different points:**

- At (0, 0, 0):  
  `φ = 0 V`  ...(1.12)

- At (−10, −10, 5):  
  `φ = 0.2332386743 V`  ...(1.13)

- At (7, 7, 7):  
  `φ = 0.624317043 V`  ...(1.14)

---

#### 1.1.2 Finding Surface Charge Density on the Infinite Plane

By using Gauss's law:

```
∯ E_T · dS = ΣQ_enclosed / ε₀      ...(1.15)
 S
```

By knowing:

```
Q_enclosed = σ × S      ...(1.16)
```

Where **σ** is the surface charge density.

We will only take the Z-component of the electrical field as it is the only component parallel to the area vector; the other components are perpendicular and zeroed by the dot product:

```
E_zT × S = (σ × S) / ε₀      ...(1.17)

σ = ε₀ · E_zT               ...(1.18)
```

Getting that:

```
E_z1 = |E₁| cos(θ) (−ẑ)     ...(1.19)
E_z2 = |E₂| cos(θ) (−ẑ)     ...(1.20)
```

Where:

```
cos(θ) = D / |r|             ...(1.21)
```

By substituting into equation (1.18) from equations (1.21), (1.19), and (1.07):

```
         1    ⎛   Q       |−Q|  ⎞   D
σ = ε₀ ———— ⎜ ———— + ———— ⎟ × ——— (−ẑ)    ...(1.22)
        4πε₀ ⎝  |r₁|²   |r₂|² ⎠   |r|
```

As y = 0, z = 0 on the x-axis:

```
|r₁| = |r₂| = |r| = √(x² + D²)     ...(1.23)
```

By updating equation (1.22):

```
       −QD
σ = —————————      ...(1.24)
      2π|r|³

          −QD
σ = ———————————————      ...(1.25)
     2π √(x² + D²)³
```

**Values of σ at different points:**

- At (0, 0, 0):  
  `σ = −1.273239545 × 10⁻¹¹ C/m²`  ...(1.26)

- At (−5, 0, 0):  
  `σ = −4.501581581 × 10⁻¹² C/m²`  ...(1.27)

- At (10, 0, 0):  
  `σ = −1.138820069 × 10⁻¹² C/m²`  ...(1.28)

---

### 1.2 MATLAB Results

**Listing 1.1:** MATLAB code for calculating the Electric Potential at certain points.

```matlab
% Auther  : Mina Wael Tanagho Fahmy Sedrak
% Section : 4
% BN      : 13
% ID      : 91240804

clc;
clear all;
close all;

%% Calculations

Q  = 2e-9;
D  = 5;
e0 = 8.854e-12;

X = [0, -10, 7];
Y = [0, -10,  7];
Z = [0,   5,  7];
V = zeros(1, 3);

for i = 1:3
    V(i) = ...
    (Q/(4*pi*e0))*((1./(sqrt(X(i).^2 + Y(i).^2 + (Z(i) - D).^2)))-...
         (1./(sqrt(X(i).^2 + Y(i).^2 + (Z(i) + D).^2))))
end

display(V);
```

MATLAB output:
```
V =
    0    0.2332    0.6243
```

We can compare between the hand analysis at equations (1.12), (1.13), and (1.14) and the MATLAB results — they agree with each other.

---

**Listing 1.2:** MATLAB code for plotting the surface charge density and marking certain values.

```matlab
% Auther  : Mina Wael Tanagho Fahmy Sedrak
% Section : 4
% BN      : 13
% ID      : 91240804

%% plot code
clc;
clear all;
close all;

Q = 2e-9;
D = 5.04;

x_axis = -100:1:100;

surface_charge_density = (-Q * D) ./ (2 * pi * (x_axis.^2 + D^2).^(1.5));

figure();
plot(x_axis, surface_charge_density, 'LineWidth', 2);
grid on;

x_point = [0, -5, 10];
y_point = (-Q * D) ./ (2 * pi * (x_point.^2 + D^2).^(1.5));
hold on;

scatter(x_point, y_point, 'filled');

for i = 1:3
    text(x_point(i), y_point(i), sprintf(' (%.f, %.4e)', x_point(i), ...
    y_point(i)), 'FontWeight', 'bold');
end

xlabel("x-axis");
ylabel("Surface Charge Density (C/m^2)");
```

We can compare between the hand analysis at equations (1.26), (1.27), and (1.28) and the MATLAB results — they agree with each other.

---

## 2 CST Solution

### 2.1 Setting Up

Starting with placing a sphere with radius = 0.05 m as a point charge having Q = 2×10⁻⁹ C and material type PEC, at a distance D = 5 m from the x-y plane.

Then placed the infinite surface from −100 to +100 in the x direction and −100 to +100 in the y direction, thickness = 1 m, with material type PEC.

Then added a straight line along the x-axis from −100 to +100 for the calculation of the surface charge density.

---

### 2.2 Finding Electrical Potential Everywhere

By running the solver and going to the 2D/3D Results, we get the potential everywhere.

Then went to the Post-Processing Result Templates and calculated the potential at:
(0, 0, 0), (−10, −10, 5), (7, 7, 7).

**Table 2.1:** CST Result Templates

| Result Name | Value |
|---|---|
| Potential (Eₛ) @ (0, 0, 0) | 0 |
| Potential (Eₛ) @ (−10, −10, 5) | 0.226764068 |
| Potential (Eₛ) @ (7, 7, 7) | 0.6174885631 |

The results are very close to the hand analysis and MATLAB results.

---

### 2.3 Finding Surface Charge Density on the Infinite Plane

After placing the infinite straight line, the surface charge density was plotted against the length in metres. Value markers were added to compare with hand analysis and MATLAB results.

The CST plot matches the MATLAB plot. Comparing values from hand analysis equations (1.26), (1.27), and (1.28) confirms they are very close.

---

## 3 Comparison Between Results

Calculating the error in the potential at point (−10, −10, 5) between CST and hand analysis:

```
         |0.226764068 − 0.2332386743|
Error = ————————————————————————————— × 100 = 2.855%    ...(3.1)
                0.226764068
```

Calculating the error in the surface charge density at point (−5, 0, 0) between CST and hand analysis:

```
         |−4.7094238×10⁻¹² − (−4.501581581×10⁻¹²)|
Error = ————————————————————————————————————————————— × 100 = 4.413%    ...(3.2)
                    −4.7094238×10⁻¹²
```

Both errors are very small and therefore acceptable.

**Reason for the errors:**

- **Math (Hand Analysis):** Infinite, continuous, and perfect.
- **CST:** Finite, meshed, and limited by the size of the simulation box.

> **Note:** The plate used in CST is not truly infinite — it has an area of 40,000 m² (200 × 200 metres).
