# 2장 — Python 기초

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **최종 수정일:** 2026-10-20
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 2, 3

> **학습 목표:**
> 1. Python과 NumPy, SciPy, Matplotlib으로 이루어진 수치 계산 환경의 구성을 설명할 수 있다
> 2. NumPy 배열로 요소별 연산과 행렬 연산을 구분하여 수행할 수 있다
> 3. 반복문, 조건문, 함수를 사용하여 수치 계산을 위한 구조화된 프로그램을 작성할 수 있다
> 4. 번지 점프 문제를 오일러 방법으로 구현하고 결과를 해석해와 비교할 수 있다

---

## 목차

- [1. 프로그래밍 환경](#1-프로그래밍-환경)
- [2. Python 기초](#2-python-기초)
  - [2.1 변수와 대입](#21-변수와-대입)
  - [2.2 데이터 타입](#22-데이터-타입)
  - [2.3 NumPy 배열](#23-numpy-배열)
- [3. 기본 연산](#3-기본-연산)
  - [3.1 요소별 연산](#31-요소별-연산)
  - [3.2 행렬 연산](#32-행렬-연산)
- [4. 유용한 NumPy/SciPy 함수](#4-유용한-numpyscipy-함수)
- [5. Matplotlib을 이용한 그래프](#5-matplotlib을-이용한-그래프)
- [6. 제어 구조](#6-제어-구조)
  - [6.1 For 반복문](#61-for-반복문)
  - [6.2 While 반복문](#62-while-반복문)
  - [6.3 조건문](#63-조건문)
- [7. 함수](#7-함수)
- [8. 예제: Python으로 구현하는 번지 점프](#8-예제-python으로-구현하는-번지-점프)
  - [8.1 오일러 방법 구현](#81-오일러-방법-구현)
  - [8.2 수치해 vs. 해석해 그래프](#82-수치해-vs-해석해-그래프)
  - [8.3 scipy.integrate.solve_ivp 사용](#83-scipyintegratesolve_ivp-사용)
- [요약](#요약)
- [점검 문제](#점검-문제)

---

<br>

## 1. 프로그래밍 환경

수치 계산에는 **Python**을 주요 프로그래밍 언어로 사용하며, 과학 계산 생태계를 함께 활용한다:

- **NumPy**: 수치 배열, 선형대수, 수학 함수를 제공한다.
- **SciPy**: 고급 과학 계산(ODE 풀이, 최적화, 보간 등)을 제공한다.
- **Matplotlib**: 2D/3D 그래프와 시각화를 제공한다.

| 특징 | MATLAB | Python |
|:--------|:-------|:-------|
| 라이선스 | 상용 (고가) | 무료 및 오픈 소스 |
| 배열 처리 | 내장 행렬 타입 | NumPy `ndarray` |
| 생태계 | 툴박스 (유료 애드온) | 패키지 (pip/conda를 통해 무료) |
| 문법 스타일 | 1 기반 인덱싱, end 키워드 | 0 기반 인덱싱, 들여쓰기 기반 |
| 커뮤니티 | 공학/학계 | 폭넓음 (데이터 과학, 웹, ML 등) |

*표 1. MATLAB과 Python의 비교*

> **참고:** 수치 해법에 대한 대부분의 교재는 역사적으로 MATLAB용으로 작성되었다. 개념과 알고리즘은 프로그래밍 언어에 관계없이 동일하며, 문법만 다르다. 교재는 알고리즘을 Python으로 구현한다.

---

<br>

## 2. Python 기초

### 2.1 변수와 대입

Python에서 변수는 단순한 대입으로 생성된다. 타입을 명시적으로 선언할 필요가 없으며, Python이 대입된 값으로부터 타입을 추론한다:

```python
x = 5        # int
y = 3.14     # float
name = "NA"  # string
```

변수 이름은 대소문자를 구분하며(`x`와 `X`는 다른 변수), 문자 또는 밑줄로 시작해야 한다.

### 2.2 데이터 타입

Python은 과학 계산과 관련된 여러 내장 데이터 타입을 제공한다:

| 타입 | 예시 | 설명 |
|:-----|:--------|:------------|
| `int` | `5`, `-3` | 정수 (임의 정밀도) |
| `float` | `3.14`, `1e-6` | 부동 소수점 수 (64비트 배정밀도) |
| `str` | `'hello'` | 문자열 |
| `list` | `[1, 2, 3]` | 가변(Mutable) 순서열 |
| `tuple` | `(1, 2, 3)` | 불변(Immutable) 순서열 |
| `dict` | `{'a': 1}` | 키-값 매핑 |

*표 2. 과학 계산과 관련된 Python 내장 데이터 타입*

수치 작업에서 Python의 내장 `list`는 벡터화 산술을 지원하지 않으므로 충분하지 않다. 이것이 **NumPy 배열**을 사용하는 이유이다.

### 2.3 NumPy 배열

NumPy의 `ndarray`는 Python에서 수치 계산을 위한 기본 데이터 구조이다. 빠르고, 메모리 효율적이며, 다차원 배열에 대해 벡터화 연산을 제공한다:

```python
import numpy as np

# 1D array (vector)
a = np.array([1, 2, 3])
print(a)       # [1 2 3]
print(a.shape) # (3,)

# 2D array (matrix)
A = np.array([[1, 2],
              [3, 4]])
print(A)
# [[1 2]
#  [3 4]]
print(A.shape) # (2, 2)
```

Python 리스트와 NumPy 배열의 주요 차이점은 다음과 같다:

| 특징 | Python `list` | NumPy `ndarray` |
|:--------|:-------------|:----------------|
| 요소 타입 | 혼합 허용 | 동일 타입(Homogeneous) |
| 산술 연산 | 연결(Concatenation) (`+`) | 요소별 덧셈 |
| 성능 | 느림 (인터프리터 루프) | 빠름 (컴파일된 C 백엔드) |
| 메모리 | 높은 오버헤드 | 컴팩트, 연속적 |

*표 3. Python 리스트와 NumPy 배열의 차이*

---

<br>

## 3. 기본 연산

### 3.1 요소별 연산

NumPy 배열에 대한 표준 산술 연산자는 **요소별(Element-wise)** 로 동작한다:

```python
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])

print(a + b)   # [5 7 9]     — element-wise addition
print(a - b)   # [-3 -3 -3]  — element-wise subtraction
print(a * b)   # [4 10 18]   — element-wise multiplication
print(a / b)   # [0.25 0.4 0.5] — element-wise division
print(a ** 2)  # [1 4 9]     — element-wise power
```

### 3.2 행렬 연산

행렬(선형대수) 연산에는 `@` 연산자 또는 전용 NumPy 함수를 사용한다:

```python
A = np.array([[1, 2],
              [3, 4]])
B = np.array([[5, 6],
              [7, 8]])

# Matrix multiplication
C = A @ B            # preferred syntax (Python 3.5+)
C = np.dot(A, B)     # equivalent

# Transpose
print(A.T)
# [[1 3]
#  [2 4]]
```

> **주의:** `*`(요소별 곱셈)과 `@`(행렬 곱셈)을 구분해야 한다. 이것은 수치 Python 코드에서 가장 흔한 버그 원인 중 하나이다.

---

<br>

## 4. 유용한 NumPy/SciPy 함수

다음 함수들은 이후의 수치 계산 전반에서 자주 사용된다.

**배열 생성:**

| 함수 | 설명 | 예시 |
|:---------|:-----------|:--------|
| `np.linspace(start, stop, num)` | 균등 간격 점 | `np.linspace(0, 1, 5)` $`\to`$ `[0, 0.25, 0.5, 0.75, 1]` |
| `np.arange(start, stop, step)` | 스텝 크기로 균등 간격 | `np.arange(0, 1, 0.2)` $`\to`$ `[0, 0.2, 0.4, 0.6, 0.8]` |
| `np.zeros(n)` | 0 배열 | `np.zeros(3)` $`\to`$ `[0, 0, 0]` |
| `np.ones(n)` | 1 배열 | `np.ones(3)` $`\to`$ `[1, 1, 1]` |
| `np.eye(n)` | $`n \times n`$ 단위 행렬 | `np.eye(3)` $`\to`$ $`I_3`$ |

*표 4. NumPy의 배열 생성 함수*

**선형대수 (`np.linalg`):**

| 함수 | 설명 |
|:---------|:-----------|
| `np.linalg.solve(A, b)` | $`Ax = b`$에서 $`x`$ 구하기 |
| `np.linalg.inv(A)` | $`A^{-1}`$ (역행렬) 계산 |
| `np.linalg.det(A)` | $`\det(A)`$ (행렬식) 계산 |
| `np.linalg.norm(x)` | 벡터/행렬 노름 계산 |
| `np.linalg.eig(A)` | 고유값과 고유벡터 계산 |

*표 5. np.linalg의 선형대수 함수*

```python
A = np.array([[2, 1],
              [5, 3]])
b = np.array([4, 7])

x = np.linalg.solve(A, b)  # Solve Ax = b
print(x)                    # [5. -6.]

print(np.linalg.det(A))    # 1.0000000000000002
print(np.linalg.inv(A))    # [[ 3. -1.] [-5.  2.]]
```

이 행렬의 행렬식은 정확히 $`2 \times 3 - 1 \times 5 = 1`$이지만, `np.linalg.det`는 LU 분해를 거쳐 부동 소수점으로 계산하므로 1.0000000000000002처럼 마지막 자리에 반올림 오차가 나타난다. 반올림 오차는 3장에서 다룬다.

---

<br>

## 5. Matplotlib을 이용한 그래프

Matplotlib은 Python의 표준 그래프 라이브러리이다. `pyplot` 인터페이스는 MATLAB과 유사한 플로팅 명령을 제공한다:

```python
import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 2 * np.pi, 100)
y = np.sin(x)

plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('sin(x)')
plt.title('Sine Function')
plt.grid(True)
plt.show()
```

같은 그림에 여러 곡선을 그릴 수 있다:

```python
plt.plot(x, np.sin(x), label='sin(x)')
plt.plot(x, np.cos(x), label='cos(x)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Trigonometric Functions')
plt.legend()
plt.grid(True)
plt.show()
```

일반적인 커스터마이징에는 선 스타일(`'-'`, `'--'`, `':'`), 색상(`'r'`, `'b'`, `'g'`), 마커(`'o'`, `'s'`, `'^'`)가 포함된다.

---

<br>

## 6. 제어 구조

### 6.1 For 반복문

`for` 반복문은 값의 시퀀스를 순회한다. `range()` 함수는 정수 시퀀스를 생성한다:

```python
# Sum of 1 to 10
total = 0
for i in range(1, 11):
    total += i
print(total)  # 55
```

`range(start, stop, step)`은 `start`부터 `step` 간격으로 `stop` 직전까지 정수를 생성하며, `stop`은 포함하지 않는다.

### 6.2 While 반복문

`while` 반복문은 조건이 `True`인 동안 반복한다:

```python
count = 10
while count > 0:
    print(count, end=' ')
    count -= 1
# Output: 10 9 8 7 6 5 4 3 2 1
```

While 반복문은 반복 횟수를 미리 알 수 없을 때 특히 유용하다(예: 수치 알고리즘에서 수렴할 때까지 반복).

### 6.3 조건문

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

## 7. 함수

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

## 8. 예제: Python으로 구현하는 번지 점프

이 절에서는 1장의 번지 점프 문제를 Python으로 구현하여, 수학적 공식화를 작동하는 코드로 변환하는 방법을 보인다.

### 8.1 오일러 방법 구현

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

### 8.2 수치해 vs. 해석해 그래프

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

### 8.3 scipy.integrate.solve_ivp 사용

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
| 프로그래밍 환경 | Python과 NumPy, SciPy, Matplotlib을 함께 사용하며, 이 조합은 MATLAB의 무료 오픈 소스 대안이다. |
| 변수 및 타입 | Python은 동적 타이핑을 사용하며, `int`, `float`, `str`, `list`와 NumPy `ndarray`를 주로 쓴다. |
| NumPy 배열 | `np.array(...)`로 벡터와 행렬을 만들며, 동일 타입의 요소를 저장하여 빠르고 메모리 효율적이다. |
| 요소별 연산 | `+`, `-`, `*`, `/`는 배열에 대해 요소별로 동작한다. |
| 행렬 연산 | `A @ B`로 행렬 곱셈을, `A.T`로 전치를 계산한다. |
| 주요 함수 | `np.linspace`, `np.zeros`, `np.ones`, `np.eye`, `np.linalg.solve` 등을 자주 사용한다. |
| 그래프 | `plt.plot(x, y)`로 그리고 `plt.xlabel(...)`, `plt.title(...)`로 꾸민 뒤 `plt.show()`로 표시한다. |
| 제어 구조 | `for`, `while`, `if`/`elif`/`else`로 반복과 분기를 제어한다. |
| 함수 | `def func(args): return result` 형태로 정의하며, 튜플로 여러 값을 반환할 수 있다. |
| ODE 풀이 | 오일러 방법은 직접 구현하는 간단한 방법이고, `scipy.integrate.solve_ivp`는 적응적 스텝 크기로 더 정확한 해를 준다. |

---

<br>

## 점검 문제

1. **linspace와 arange:** `np.linspace(0, 2, 5)`와 `np.arange(0, 2, 0.5)`의 결과를 쓰라.

   <details>
   <summary>정답</summary>

   `np.linspace(0, 2, 5)`는 끝점 2를 포함한 5개의 점 `[0, 0.5, 1, 1.5, 2]`를 만든다. `np.arange(0, 2, 0.5)`는 `stop`인 2를 포함하지 않으므로 `[0, 0.5, 1, 1.5]`이다.

   </details>

2. **리스트와 배열의 `+`:** `[1, 2] + [3, 4]`와 `np.array([1, 2]) + np.array([3, 4])`의 결과를 쓰라.

   <details>
   <summary>정답</summary>

   리스트의 `+`는 연결이므로 `[1, 2, 3, 4]`이다. NumPy 배열의 `+`는 요소별 덧셈이므로 `[4 6]`이다.

   </details>

3. **요소별 곱셈과 행렬 곱셈:** 3.2절의 `A`와 `B`에 대해 `A * B`와 `A @ B`를 구하라.

   <details>
   <summary>정답</summary>

   `A * B`는 요소별 곱이므로 `[[5, 12], [21, 32]]`이다. `A @ B`는 행렬 곱이므로 1×5 + 2×7 = 19, 1×6 + 2×8 = 22, 3×5 + 4×7 = 43, 3×6 + 4×8 = 50에서 `[[19, 22], [43, 50]]`이다.

   </details>

4. **연립방정식 풀이:** $`3x_1 + 2x_2 = 7`$, $`x_1 + 2x_2 = 5`$를 `np.linalg.solve`로 풀 때의 입력과 결과를 쓰라.

   <details>
   <summary>정답</summary>

   `A = np.array([[3, 2], [1, 2]])`, `b = np.array([7, 5])`로 두고 `np.linalg.solve(A, b)`를 호출하면 `[1. 2.]`가 나온다. 즉 $`x_1 = 1`$, $`x_2 = 2`$이며, $`3 + 4 = 7`$, $`1 + 4 = 5`$로 확인할 수 있다.

   </details>

5. **오일러 코드 추적:** 8.1절의 코드가 $`t = 2`$ s와 $`t = 4`$ s에서 출력하는 속도를 구하라.

   <details>
   <summary>정답</summary>

   $`v_1 = 0 + (9.81 - 0) \times 2 = 19.62`$이다. $`v_2 = 19.62 + (9.81 - (12.5/68.1) \times 19.62) \times 2 = 19.62 + 6.2087 \times 2 \approx 32.0374`$이다. 따라서 출력은 19.6200과 32.0374이다.

   </details>
