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

*Figure 1.1: The problem illustration.*

### 1.1 Hand Analysis

#### 1.1.1 Finding Electrical Potential Everywhere

We will start the hand analysis by performing the method of images.

*Figure 1.2: Method of Images.*

We know that the Electrical Field of point charge (Q) at any point is:

$$\vec{E} = \frac{Q}{4\pi\varepsilon_0 r^2} \hat{r} \tag{1.01}$$

For the main point charge (Q), the Electrical Field is:

$$\vec{E_1} = \frac{Q}{4\pi\varepsilon_0 r_1^2} \hat{r_1} \tag{1.02}$$

Where:

$$|r_1| = \sqrt{x^2 + y^2 + (z-D)^2} \tag{1.03}$$

For the image point charge (–Q), the Electrical Field is:

$$\vec{E_2} = \frac{-Q}{4\pi\varepsilon_0 r_2^2} \hat{r_2} \tag{1.04}$$

Where:

$$|r_2| = \sqrt{x^2 + y^2 + (z+D)^2} \tag{1.05}$$

To get the total Electrical Field:

$$\vec{E_T} = \vec{E_1} + \vec{E_2} \tag{1.06}$$

$$E_T = \frac{Q}{4\pi\varepsilon_0} \left(\frac{1}{r_1^2} - \frac{1}{r_2^2}\right) \hat{E} \tag{1.07}$$

We know that:

$$\varphi = -\int E_T \cdot dl \tag{1.08}$$

$$\varphi = -\int \frac{Q}{4\pi\varepsilon_0} \left(\frac{1}{r_1^2} - \frac{1}{r_2^2}\right) dr \tag{1.09}$$

$$\varphi = \frac{Q}{4\pi\varepsilon_0} \left(\frac{1}{r_1} - \frac{1}{r_2}\right) \tag{1.10}$$

By substituting, we get the potential everywhere:

$$\varphi = \frac{Q}{4\pi\varepsilon_0} \left(\frac{1}{\sqrt{x^2+y^2+(z-D)^2}} - \frac{1}{\sqrt{x^2+y^2+(z+D)^2}}\right) \tag{1.11}$$

Getting values of **φ** at different points:

- At (0, 0, 0):
$$\varphi = 0 \text{ V} \tag{1.12}$$

- At (−10, −10, 5):
$$\varphi = 0.2332386743 \text{ V} \tag{1.13}$$

- At (7, 7, 7):
$$\varphi = 0.624317043 \text{ V} \tag{1.14}$$

---

#### 1.1.2 Finding Surface Charge Density on the Infinite Plane

*Figure 1.3: Added an infinite straight line below the point charge to calculate surface charge density on it.*

By using Gauss's law:

$$\oiint_S \vec{E_T} \cdot d\vec{S} = \frac{\sum Q_{enclosed}}{\varepsilon_0} \tag{1.15}$$

By knowing:

$$Q_{enclosed} = \sigma \times S \tag{1.16}$$

Where **σ** is the surface charge density.

We will only take the Z-component of the electrical field as it is the only component parallel to the area vector; the other components are perpendicular to it and are zeroed by the dot product:

$$E_{z_T} \times S = \frac{\sigma \times S}{\varepsilon_0} \tag{1.17}$$

$$\sigma = \varepsilon_0 \, E_{z_T} \tag{1.18}$$

*Figure 1.4: Getting the Z component of E₁ and E₂.*

Getting that:

$$\vec{E}_{z_1} = |E_1| \cos(\theta)(-\hat{z}) \tag{1.19}$$

$$\vec{E}_{z_2} = |E_2| \cos(\theta)(-\hat{z}) \tag{1.20}$$

Where:

$$\cos(\theta) = \frac{D}{|r|} \tag{1.21}$$

By substituting into equation (1.18) from equations (1.21), (1.19), and (1.07):

$$\sigma = \varepsilon_0 \frac{1}{4\pi\varepsilon_0} \left(\frac{Q}{|r_1|^2} + \frac{|-Q|}{|r_2|^2}\right) \frac{D}{|r|} (-\hat{z}) \tag{1.22}$$

As y = 0, z = 0 on the x-axis:

$$|r_1| = |r_2| = |r| = \sqrt{x^2 + D^2} \tag{1.23}$$

By updating equation (1.22):

$$\sigma = \frac{-QD}{2\pi|r|^3} \tag{1.24}$$

$$\sigma = \frac{-QD}{2\pi\sqrt{(x^2+D^2)^3}} \tag{1.25}$$

Getting values of **σ** at different points:

- At (0, 0, 0):
$$\sigma = -1.273239545 \times 10^{-11} \text{ C/m}^2 \tag{1.26}$$

- At (−5, 0, 0):
$$\sigma = -4.501581581 \times 10^{-12} \text{ C/m}^2 \tag{1.27}$$

- At (10, 0, 0):
$$\sigma = -1.138820069 \times 10^{-12} \text{ C/m}^2 \tag{1.28}$$

---

### 1.2 MATLAB Results

**Listing 1.1:** *MATLAB code for calculating the Electric Potential at certain points.*

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
```

*Figure 1.5: MATLAB results part 1.*

We can compare between the hand analysis at equations (1.12), (1.13) and (1.14) and the MATLAB results; they agree with each other.

**Listing 1.2:** *MATLAB code for plotting the surface charge density and marking certain values.*

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

*Figure 1.6: MATLAB results part 2.*

We can compare between the hand analysis at equations (1.26), (1.27) and (1.28) and the MATLAB results; they agree with each other.

---

## 2 CST Solution

### 2.1 Setting Up

Starting with placing a sphere with radius equals 0.05 m as a point charge having a charge Q = 2×10⁻⁹ C and a material type of PEC at a distance D = 5 m from the x-y plane.

*Figure 2.1: Placing the point charge.*

Then placed the infinite surface from −100 to +100 in the x direction and from −100 to +100 in the y direction, with a thickness of 1 m and a material type of PEC.

*Figure 2.2: Placing the infinite plate.*

Then added the straight line for the calculation of the surface charge density at the x-axis from −100 to +100.

*Figure 2.3: Placing the infinite straight line.*

---

### 2.2 Finding Electrical Potential Everywhere

By running the Set up solver and going to the 2D/3D Results, we get the potential everywhere as shown in **Figure 2.4** and **Figure 2.5**.

*Figure 2.4: Potential Everywhere.*

*Figure 2.5: Potential Everywhere Section view.*

Then went to the Post-Processing Result Templates and calculated the potential at certain points: (0, 0, 0), (−10, −10, 5), (7, 7, 7).

**Table 2.1:** CST Result Templates

| Result Name | Value |
|---|---|
| Potential (Eₛ) @ (0, 0, 0) | 0 |
| Potential (Eₛ) @ (−10, −10, 5) | 0.226764068 |
| Potential (Eₛ) @ (7, 7, 7) | 0.6174885631 |

As we can see, the results are very close to the hand analysis and MATLAB results.

---

### 2.3 Finding Surface Charge Density on the Infinite Plane

After placing the infinite straight line, we went to the Post-Processing Result Templates and plotted the surface charge density graph against the length in metres, as in **Figure 2.7**. Then added some value markers to compare with hand analysis and MATLAB results, as in **Figure 2.8**.

*Figure 2.7: Plot of surface charge density vs length in m.*

*Figure 2.8: Marking the values we are interested in and enabling the legend.*

As we can see, the plot is the same as the MATLAB plot. Also, by comparing values from hand analysis equations at (1.26), (1.27) and (1.28), they are very close.

---

## 3 Comparison Between Results

Calculating the error in the potential result at point (−10, −10, 5) from CST and hand analysis:

$$Error\ (\%) = \frac{|0.226764068 - 0.2332386743|}{0.226764068} \times 100 = 2.855\% \tag{3.1}$$

Calculating the error in the surface charge density result at point (−5, 0, 0) from CST and hand analysis:

$$Error\ (\%) = \frac{|-4.7094238 \times 10^{-12} - (-4.501581581 \times 10^{-12})|}{-4.7094238 \times 10^{-12}} \times 100 = 4.413\% \tag{3.2}$$

Both errors are very small and therefore acceptable.

The reason these errors occur:

- **Math (Hand Analysis):** Infinite, continuous, and perfect.
- **CST:** Finite, meshed, and limited by the size of the box.

> **Note:** The plate used in CST is not infinite — it has an area of 40,000 m² (200 × 200 metres).
