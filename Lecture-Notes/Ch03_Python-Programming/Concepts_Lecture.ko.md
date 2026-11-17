# 3장 — Python 프로그래밍

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **최종 수정일:** 2026-11-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 3

> **학습 목표:**
> 1. 반복문, 조건문, 함수를 사용하여 수치 계산을 위한 구조화된 프로그램을 작성할 수 있다
> 2. 번지 점프 문제를 오일러 방법으로 구현하고 결과를 해석해와 비교할 수 있다
> 3. SciPy의 `solve_ivp`로 같은 ODE를 풀어 오일러 방법, 해석해와 함께 그래프로 비교할 수 있다

---

## 목차

- [1. 제어 구조](#1-제어-구조)
  - [1.1 For 반복문](#11-for-반복문)
  - [1.2 While 반복문](#12-while-반복문)
  - [1.3 조건문](#13-조건문)
- [2. 함수](#2-함수)
- [3. 예제: Python으로 구현하는 번지 점프](#3-예제-python으로-구현하는-번지-점프)
  - [3.1 오일러 방법 구현](#31-오일러-방법-구현)
  - [3.2 수치해 vs. 해석해 그래프](#32-수치해-vs-해석해-그래프)
  - [3.3 scipy.integrate.solve_ivp 사용](#33-scipyintegratesolve_ivp-사용)
- [요약](#요약)
- [점검 문제](#점검-문제)

---

<br>

## 1. 제어 구조

### 1.1 For 반복문

`for` 반복문은 값의 시퀀스를 순회한다. `range()` 함수는 정수 시퀀스를 생성한다:

```python
# Sum of 1 to 10
total = 0
for i in range(1, 11):
    total += i
print(total)  # 55
```

`range(start, stop, step)`은 `start`부터 `step` 간격으로 `stop` 직전까지 정수를 생성하며, `stop`은 포함하지 않는다.

### 1.2 While 반복문

`while` 반복문은 조건이 `True`인 동안 반복한다:

```python
count = 10
while count > 0:
    print(count, end=' ')
    count -= 1
# Output: 10 9 8 7 6 5 4 3 2 1
```

While 반복문은 반복 횟수를 미리 알 수 없을 때 특히 유용하다(예: 수치 알고리즘에서 수렴할 때까지 반복).

### 1.3 조건문

`if`/`elif`/`else` 구조는 분기를 제어한다:

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

## 2. 함수

함수는 재사용 가능한 코드 블록을 캡슐화한다. `def` 키워드로 정의한다:

```python
def myfunction(x):
    return x ** 2

print(myfunction(5))  # 25
```

함수는 튜플로 여러 값을 반환할 수 있다:

```python
def circle_properties(r):
    area = 3.14159 * r ** 2
    circumference = 2 * 3.14159 * r
    return area, circumference

a, c = circle_properties(5)
print(f'Area = {a:.2f}, Circumference = {c:.2f}')
# Area = 78.54, Circumference = 31.42
```

기본 인수를 사용하면 유연한 함수 호출이 가능하다:

```python
def greet(name, greeting='Hello'):
    return f'{greeting}, {name}!'

print(greet('World'))            # Hello, World!
print(greet('World', 'Hi'))      # Hi, World!
```

---

<br>

## 3. 예제: Python으로 구현하는 번지 점프

이 절에서는 1장의 번지 점프 문제를 Python으로 구현하여, 수학적 공식화를 작동하는 코드로 변환하는 방법을 보인다.

### 3.1 오일러 방법 구현

선형 항력을 가진 번지 점프의 ODE는 다음과 같다:

```math
\frac{dv}{dt} = g - \frac{c}{m}v
```

오일러 방법을 적용하면 다음과 같다:

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

이 코드는 1장 3.2절의 손 계산과 같은 값을 출력하며, 마지막 줄은 $`t = 12`$ s에서 50.0102 m/s이다.

### 3.2 수치해 vs. 해석해 그래프

1장에서 구한 해석해는 다음과 같다:

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

### 3.3 scipy.integrate.solve_ivp 사용

더 정확한 ODE 풀이를 위해 SciPy는 적응적 스텝 크기 방법을 사용하는 `solve_ivp`를 제공한다:

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

> **참고:** `scipy.integrate.solve_ivp`는 기본적으로 적응적 룽게-쿠타 방법(RK45)을 사용하며, 오일러 방법보다 훨씬 정확하다. 룽게-쿠타 방법은 22장에서, 적응적 스텝 크기 제어는 23장에서 자세히 다룬다.

---

<br>

## 요약

| 개념 | 핵심 요약 |
|:------|:----------|
| 제어 구조 | `for`, `while`, `if`/`elif`/`else`로 반복과 분기를 제어한다. |
| 함수 | `def func(args): return result` 형태로 정의하며, 튜플로 여러 값을 반환할 수 있다. |
| ODE 풀이 | 오일러 방법은 직접 구현하는 간단한 방법이고, `scipy.integrate.solve_ivp`는 적응적 스텝 크기로 더 정확한 해를 준다. |

---

<br>

## 점검 문제

1. **range의 끝값:** `range(2, 11, 3)`이 만드는 정수를 모두 쓰라.

   <details>
   <summary>정답</summary>

   `range(start, stop, step)`은 `start`인 2부터 `step`인 3씩 더하며 `stop`인 11 직전까지 정수를 만든다. 따라서 결과는 2, 5, 8이며, 그다음 값 11은 `stop`과 같으므로 포함되지 않는다.

   </details>

2. **오일러 코드 추적:** 3.1절의 코드가 $`t = 2`$ s와 $`t = 4`$ s에서 출력하는 속도를 구하라.

   <details>
   <summary>정답</summary>

   $`v_1 = 0 + (9.81 - 0) \times 2 = 19.62`$이다. $`v_2 = 19.62 + (9.81 - (12.5/68.1) \times 19.62) \times 2 = 19.62 + 6.2087 \times 2 \approx 32.0374`$이다. 따라서 출력은 19.6200과 32.0374이다.

   </details>
