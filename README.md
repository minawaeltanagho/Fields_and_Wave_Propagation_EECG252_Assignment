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

> **Assuming $Q = 2\ \text{nC},\ D = 5\ \text{m}$**

### 1.1 Hand Analysis

#### 1.1.1 Finding Electrical Potential Everywhere

We will start the hand analysis by performing the method of images.

We know that the Electrical Field of a point charge $Q$ at any point is:

$$\vec{E} = \frac{Q}{4\pi\varepsilon_0 r^2}\,\hat{r} \tag{1.01}$$

For the main point charge $+Q$:

$$\vec{E}_1 = \frac{Q}{4\pi\varepsilon_0 r_1^2}\,\hat{r}_1 \tag{1.02}$$

$$|r_1| = \sqrt{x^2 + y^2 + (z-D)^2} \tag{1.03}$$

For the image point charge $-Q$:

$$\vec{E}_2 = \frac{-Q}{4\pi\varepsilon_0 r_2^2}\,\hat{r}_2 \tag{1.04}$$

$$|r_2| = \sqrt{x^2 + y^2 + (z+D)^2} \tag{1.05}$$

Total Electrical Field:

$$\vec{E}_T = \vec{E}_1 + \vec{E}_2 \tag{1.06}$$

$$E_T = \frac{Q}{4\pi\varepsilon_0}\left(\frac{1}{r_1^2} - \frac{1}{r_2^2}\right)\hat{E} \tag{1.07}$$

Since $\varphi = -\int \vec{E}_T \cdot d\vec{l}$:

$$\varphi = -\int \frac{Q}{4\pi\varepsilon_0}\left(\frac{1}{r_1^2} - \frac{1}{r_2^2}\right)dr \tag{1.09}$$

$$\varphi = \frac{Q}{4\pi\varepsilon_0}\left(\frac{1}{r_1} - \frac{1}{r_2}\right) \tag{1.10}$$

Substituting, the potential everywhere is:

$$\boxed{\varphi = \frac{Q}{4\pi\varepsilon_0}\left(\frac{1}{\sqrt{x^2+y^2+(z-D)^2}} - \frac{1}{\sqrt{x^2+y^2+(z+D)^2}}\right)} \tag{1.11}$$

**Values of $\varphi$ at different points:**

| Point | $\varphi$ (V) | Eq. |
|-------|--------------|-----|
| $(0,\ 0,\ 0)$ | $0$ | (1.12) |
| $(-10,\ -10,\ 5)$ | $0.2332386743$ | (1.13) |
| $(7,\ 7,\ 7)$ | $0.624317043$ | (1.14) |

---

#### 1.1.2 Finding Surface Charge Density on the Infinite Plane

By Gauss's law:

$$\oiint_S \vec{E}_T \cdot d\vec{S} = \frac{\sum Q_{\text{enclosed}}}{\varepsilon_0} \tag{1.15}$$

$$Q_{\text{enclosed}} = \sigma \times S \tag{1.16}$$

where $\sigma$ is the surface charge density. Taking only the Z-component (perpendicular to the plane):

$$E_{z_T} \times S = \frac{\sigma \times S}{\varepsilon_0} \tag{1.17}$$

$$\sigma = \varepsilon_0\, E_{z_T} \tag{1.18}$$

The Z-components from each charge are:

$$\vec{E}_{z_1} = |E_1|\cos\theta\,(-\hat{z}) \tag{1.19}$$

$$\vec{E}_{z_2} = |E_2|\cos\theta\,(-\hat{z}) \tag{1.20}$$

$$\cos\theta = \frac{D}{|r|} \tag{1.21}$$

Substituting equations (1.21), (1.19), and (1.07) into (1.18):

$$\sigma = \varepsilon_0 \cdot \frac{1}{4\pi\varepsilon_0}\left(\frac{Q}{|r_1|^2} + \frac{|-Q|}{|r_2|^2}\right)\frac{D}{|r|}(-\hat{z}) \tag{1.22}$$

On the x-axis where $y=0,\ z=0$:

$$|r_1| = |r_2| = |r| = \sqrt{x^2 + D^2} \tag{1.23}$$

Simplifying:

$$\sigma = \frac{-QD}{2\pi|r|^3} \tag{1.24}$$

$$\boxed{\sigma = \frac{-QD}{2\pi\sqrt{(x^2+D^2)^3}}} \tag{1.25}$$

**Values of $\sigma$ at different points:**

| Point | $\sigma\ (\text{C/m}^2)$ | Eq. |
|-------|--------------------------|-----|
| $(0,\ 0,\ 0)$ | $-1.273239545 \times 10^{-11}$ | (1.26) |
| $(-5,\ 0,\ 0)$ | $-4.501581581 \times 10^{-12}$ | (1.27) |
| $(10,\ 0,\ 0)$ | $-1.138820069 \times 10^{-12}$ | (1.28) |

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
         (1./(sqrt(X(i).^2 + Y(i).^2 + (Z(i) + D).^2))));
end

display(V);
```

**Output:**
```
V =
    0    0.2332    0.6243
```

The MATLAB results match the hand analysis values from equations (1.12), (1.13), and (1.14).

---

**Listing 1.2:** MATLAB code for plotting the surface charge density.

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

The MATLAB results match the hand analysis values from equations (1.26), (1.27), and (1.28).

---

## 2 CST Solution

### 2.1 Setting Up

- Placed a **sphere** (radius = 0.05 m) as a point charge with $Q = 2 \times 10^{-9}\ \text{C}$, material type **PEC**, at distance $D = 5\ \text{m}$ from the x-y plane.
- Placed an **infinite surface** spanning $x \in [-100, +100]$, $y \in [-100, +100]$, thickness = 1 m, material type **PEC**.
- Added a **straight line** along the x-axis from $-100$ to $+100$ for surface charge density calculation.

---

### 2.2 Finding Electrical Potential Everywhere

After running the solver, the potential was evaluated at three points using Post-Processing Result Templates:

**Table 2.1:** CST Potential Results

| Result Name | Value (V) |
|---|---|
| Potential $(E_s)$ @ $(0,\ 0,\ 0)$ | $0$ |
| Potential $(E_s)$ @ $(-10,\ -10,\ 5)$ | $0.226764068$ |
| Potential $(E_s)$ @ $(7,\ 7,\ 7)$ | $0.6174885631$ |

Results are in close agreement with hand analysis and MATLAB.

---

### 2.3 Finding Surface Charge Density on the Infinite Plane

The surface charge density was plotted along the x-axis using the Post-Processing Result Templates. Value markers were added to compare against hand analysis and MATLAB results. The CST plot closely matches the MATLAB plot.

---

## 3 Comparison Between Results

**Error in potential** at point $(-10,\ -10,\ 5)$:

$$\text{Error}\ (\%) = \frac{|0.226764068 - 0.2332386743|}{0.226764068} \times 100 = 2.855\% \tag{3.1}$$

**Error in surface charge density** at point $(-5,\ 0,\ 0)$:

$$\text{Error}\ (\%) = \frac{|-4.7094238\times10^{-12} - (-4.501581581\times10^{-12})|}{|-4.7094238\times10^{-12}|} \times 100 = 4.413\% \tag{3.2}$$

Both errors are very small and therefore **acceptable**.

**Reason for the errors:**

| Source | Nature |
|--------|--------|
| Hand Analysis | Infinite, continuous, and perfect |
| CST Simulation | Finite, meshed, limited by simulation box size |

> **Note:** The plate in CST is not truly infinite — it has an area of $40{,}000\ \text{m}^2$ ($200 \times 200$ metres).
