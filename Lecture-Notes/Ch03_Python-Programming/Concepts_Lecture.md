# Chapter 3 — Python Programming

[Repository Guide](../../README.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **Last Updated:** 2026-11-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 3

> **Learning Objectives:**
> 1. Write structured programs for numerical computation with loops, conditionals, and functions
> 2. Implement the bungee jumper problem with Euler's method and compare the result with the analytical solution
> 3. Solve the same ODE with SciPy's `solve_ivp` and compare it in one plot with Euler's method and the analytical solution

---

## Table of Contents

- [1. Control Structures](#1-control-structures)
  - [1.1 For Loop](#11-for-loop)
  - [1.2 While Loop](#12-while-loop)
  - [1.3 Conditional Statements](#13-conditional-statements)
- [2. Functions](#2-functions)
- [3. Example: Bungee Jumper in Python](#3-example-bungee-jumper-in-python)
  - [3.1 Implementing Euler's Method](#31-implementing-eulers-method)
  - [3.2 Plotting Numerical vs. Analytical Solution](#32-plotting-numerical-vs-analytical-solution)
  - [3.3 Using scipy.integrate.solve_ivp](#33-using-scipyintegratesolve_ivp)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Control Structures

### 1.1 For Loop

The `for` loop iterates over a sequence of values. The `range()` function generates integer sequences:

```python
# Sum of 1 to 10
total = 0
for i in range(1, 11):
    total += i
print(total)  # 55
```

`range(start, stop, step)` generates integers from `start` in increments of `step` up to (but not including) `stop`.

### 1.2 While Loop

The `while` loop repeats as long as a condition is `True`:

```python
count = 10
while count > 0:
    print(count, end=' ')
    count -= 1
# Output: 10 9 8 7 6 5 4 3 2 1
```

While loops are especially useful when the number of iterations is not known in advance (e.g., iterating until convergence in a numerical algorithm).

### 1.3 Conditional Statements

The `if`/`elif`/`else` structure controls branching:

```python
x = 42

if x > 0:
    print('Positive')
elif x == 0:
    print('Zero')
else:
    print('Negative')
```

---

<br>

## 2. Functions

Functions encapsulate reusable blocks of code. They are defined with the `def` keyword:

```python
def myfunction(x):
    return x ** 2

print(myfunction(5))  # 25
```

Functions can return multiple values as a tuple:

```python
def circle_properties(r):
    area = 3.14159 * r ** 2
    circumference = 2 * 3.14159 * r
    return area, circumference

a, c = circle_properties(5)
print(f'Area = {a:.2f}, Circumference = {c:.2f}')
# Area = 78.54, Circumference = 31.42
```

Default arguments allow flexible function calls:

```python
def greet(name, greeting='Hello'):
    return f'{greeting}, {name}!'

print(greet('World'))            # Hello, World!
print(greet('World', 'Hi'))      # Hi, World!
```

---

<br>

## 3. Example: Bungee Jumper in Python

This section implements the bungee jumper problem from Chapter 1 in Python, demonstrating how to translate mathematical formulations into working code.

### 3.1 Implementing Euler's Method

The ODE for a bungee jumper with linear drag is:

```math
\frac{dv}{dt} = g - \frac{c}{m}v
```

Applying Euler's method gives:

```math
v(t_{i+1}) = v(t_i) + \left(g - \frac{c}{m}v(t_i)\right) \cdot h
```

```python
import numpy as np

# Parameters
g = 9.81      # m/s^2
m = 68.1      # kg
c = 12.5      # kg/s
h = 2.0       # step size (s)
t_end = 12.0  # final time (s)

# Initialize
t = np.arange(0, t_end + h, h)
v = np.zeros(len(t))
v[0] = 0  # initial condition

# Euler's method
for i in range(len(t) - 1):
    dvdt = g - (c / m) * v[i]
    v[i + 1] = v[i] + dvdt * h

print('t (s)  |  v (m/s)')
for ti, vi in zip(t, v):
    print(f'{ti:5.1f}  |  {vi:.4f}')
```

This code prints the same values as the hand calculation in Section 3.2 of Chapter 1; the last line gives 50.0102 m/s at $`t = 12`$ s.

### 3.2 Plotting Numerical vs. Analytical Solution

The analytical solution derived in Chapter 1 is:

```math
v(t) = \frac{gm}{c}\left(1 - e^{-(c/m)t}\right)
```

```python
import matplotlib.pyplot as plt

# Analytical solution on a fine grid
t_fine = np.linspace(0, t_end, 200)
v_analytical = (g * m / c) * (1 - np.exp(-(c / m) * t_fine))

# Plot comparison
plt.plot(t_fine, v_analytical, 'b-', label='Analytical')
plt.plot(t, v, 'ro--', label=f'Euler (h={h}s)')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Bungee Jumper: Euler vs. Analytical')
plt.legend()
plt.grid(True)
plt.show()
```

### 3.3 Using scipy.integrate.solve_ivp

For more accurate ODE solving, SciPy provides `solve_ivp`, which uses adaptive step-size methods:

```python
from scipy.integrate import solve_ivp

# Define the ODE as a function: dv/dt = f(t, v)
def bungee_ode(t, v):
    return g - (c / m) * v

# Solve
sol = solve_ivp(bungee_ode, [0, t_end], [0], t_eval=t_fine)

plt.plot(t_fine, v_analytical, 'b-', label='Analytical')
plt.plot(sol.t, sol.y[0], 'g--', label='solve_ivp (RK45)')
plt.plot(t, v, 'ro--', label=f'Euler (h={h}s)')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Bungee Jumper: Comparison of Methods')
plt.legend()
plt.grid(True)
plt.show()
```

> **Note:** `scipy.integrate.solve_ivp` uses an adaptive Runge-Kutta method (RK45) by default, which is much more accurate than Euler's method. Runge-Kutta methods are covered in detail in Chapter 22, and adaptive step-size control in Chapter 23.

---

<br>

## Summary

| Concept | Key Summary |
|:------|:----------|
| Control Structures | `for`, `while`, and `if`/`elif`/`else` control repetition and branching. |
| Functions | Functions are defined as `def func(args): return result` and can return multiple values as a tuple. |
| ODE Solving | Euler's method is a simple method implemented by hand, while `scipy.integrate.solve_ivp` gives a more accurate solution with adaptive steps. |

---

<br>

## Self-Check Questions

1. **End value of range:** Write all the integers that `range(2, 11, 3)` generates.

   <details>
   <summary>Answer</summary>

   `range(start, stop, step)` generates integers starting from `start`, which is 2, in steps of `step`, which is 3, up to but not including `stop`, which is 11. The result is therefore 2, 5, and 8, and the next value 11 is not included because it equals `stop`.

   </details>

2. **Tracing the Euler code:** Find the velocities that the code in Section 3.1 prints at $`t = 2`$ s and $`t = 4`$ s.

   <details>
   <summary>Answer</summary>

   $`v_1 = 0 + (9.81 - 0) \times 2 = 19.62`$. $`v_2 = 19.62 + (9.81 - (12.5/68.1) \times 19.62) \times 2 = 19.62 + 6.2087 \times 2 \approx 32.0374`$. The printed values are therefore 19.6200 and 32.0374.

   </details>
