# 제1장 실습 — 과학 계산을 위한 Python 기초

> **최종 수정일:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 1

> **선수 지식**: [프로그래밍언어] MATLAB/Python. 수치해석 사전 지식 불필요. [미적분학] 기본 미분과 적분.
>
> **학습 목표**:
> 1. 공학에서 수학적 모델링의 역할을 설명할 수 있다
> 2. 수치 계산에서 오차의 원인을 식별할 수 있다
> 3. 모델 복잡도와 정확도 간의 트레이드오프를 설명할 수 있다

---

<br>

## 목차

- [1. 변수와 대입](#1-변수와-대입)
  - [1.1 동적 타이핑](#11-동적-타이핑)
  - [1.2 수학 상수와 서식 출력](#12-수학-상수와-서식-출력)
- [2. NumPy 배열, 벡터, 행렬](#2-numpy-배열-벡터-행렬)
  - [2.1 왜 NumPy인가?](#21-왜-numpy인가)
  - [2.2 배열 생성](#22-배열-생성)
  - [2.3 인덱싱과 슬라이싱](#23-인덱싱과-슬라이싱)
  - [2.4 시퀀스 생성: arange, linspace, logspace](#24-시퀀스-생성-arange-linspace-logspace)
- [3. 수학 연산](#3-수학-연산)
  - [3.1 산술 연산자와 우선순위](#31-산술-연산자와-우선순위)
  - [3.2 요소별 곱셈 vs. 행렬 곱셈](#32-요소별-곱셈-vs-행렬-곱셈)
  - [3.3 내적, 외적, 전치](#33-내적-외적-전치)
- [4. 내장 함수](#4-내장-함수)
  - [4.1 주요 내장 함수](#41-주요-내장-함수)
  - [4.2 은행원 반올림](#42-은행원-반올림)
  - [4.3 응용: 번지 점프 속도](#43-응용-번지-점프-속도)
- [5. Matplotlib을 이용한 시각화](#5-matplotlib을-이용한-시각화)
  - [5.1 기본 선 그래프](#51-기본-선-그래프)
  - [5.2 범례가 있는 다중 곡선](#52-범례가-있는-다중-곡선)
- [6. NumPy에서의 뷰와 복사본](#6-numpy에서의-뷰와-복사본)
  - [6.1 뷰 (슬라이싱)](#61-뷰-슬라이싱)
  - [6.2 복사본 (팬시 인덱싱)](#62-복사본-팬시-인덱싱)
- [요약](#요약)

---

<br>

## 1. 변수와 대입

### 1.1 동적 타이핑

Python은 **동적 타이핑(Dynamically Typed)** 언어로, 변수를 사용하기 전에 타입을 선언할 필요가 없다. 타입은 대입하는 값에 따라 실행 시점에 자동으로 결정된다. 변수 이름은 **대소문자를 구분**한다 (`myVar`와 `myvar`는 서로 다른 변수이다).

내장 함수 `type()`을 사용하여 변수의 타입을 확인할 수 있다. Python은 한 줄에서 **다중 대입**도 지원한다:

```python
a, b, c = 1, 2, 3  # assigns 1 to a, 2 to b, 3 to c
```

```python
a = 4
print(type(a))  # <class 'int'>

c = 2 + 4j
print(type(c))  # <class 'complex'>
```

> **[Python]** Python의 동적 타이핑은 변수에 명시적 타입 선언이 필요 없음을 의미한다. 타입은 대입된 값에 따라 실행 시점에 결정된다. 이는 C++이나 Java와 같은 정적 타입 언어와 다르다.

### 1.2 수학 상수와 서식 출력

Python의 `math` 모듈은 기본적인 수학 상수와 함수를 제공한다. 과학 계산에 사용되는 물리 상수에는 `scipy.constants` 모듈이 매우 유용하다.

```python
import math
print(math.pi)   # 3.141592653589793
print(math.e)    # 2.718281828459045

# Formatted output using str.format()
# {0:7.4f} means: argument index 0, total width 7, 4 decimal places, fixed-point
print('{0:7.4f}'.format(math.pi))  # ' 3.1416'

# Physical constants from scipy
import scipy.constants as pc
print('{0:12.4e}'.format(pc.h))  # Planck's constant: 6.6261e-34
```

> **참고:** 서식 지정자 `{0:7.4f}`는 다음과 같이 분해된다: `0` = 첫 번째 인수, `7` = 최소 전체 너비(소수점 포함), `.4` = 소수점 이하 4자리, `f` = 고정 소수점 표기법. 과학적 표기법에는 `f` 대신 `e`를 사용한다.

---

<br>

## 2. NumPy 배열, 벡터, 행렬

### 2.1 왜 NumPy인가?

**NumPy**(Numerical Python)는 Python에서 과학 계산을 위한 기반 라이브러리이다. 주요 장점은 다음과 같다:

- **속도**: NumPy 배열은 C 코드로 지원되며 연속 메모리 블록에 저장되므로, Python 리스트보다 수십~수백 배 빠르다
- **벡터화(Vectorization)**: 명시적인 Python 루프 없이 전체 배열에 요소별로 연산이 적용된다
- **브로드캐스팅(Broadcasting)**: 잘 정의된 규칙에 따라 서로 다른 형상의 배열 간 산술 연산을 허용한다

> **[Python]** NumPy 배열은 (Python 리스트와 달리) 연속 메모리 블록에 저장되어 C 속도로 실행되는 벡터화 연산을 가능하게 한다. 브로드캐스팅은 명시적 루프 없이 서로 다른 형상의 배열 간 연산을 허용한다.

### 2.2 배열 생성

NumPy는 배열을 생성하는 여러 함수를 제공한다:

| 함수 | 설명 | 예시 |
|:---------|:-----------|:--------|
| `np.array(list)` | Python 리스트로부터 배열 생성 | `np.array([1, 2, 3])` |
| `np.zeros(shape)` | 0으로 채워진 배열 | `np.zeros((3, 4))` |
| `np.ones(shape)` | 1로 채워진 배열 | `np.ones((2, 3))` |
| `np.eye(n)` | 단위 행렬 ($n \times n$) | `np.eye(3)` |
| `np.arange(start, stop, step)` | 균등 간격 값 (stop 미포함) | `np.arange(0, 1, 0.1)` |
| `np.linspace(start, stop, n)` | $n$개의 균등 간격 값 (stop 포함) | `np.linspace(0, 1, 11)` |
| `np.logspace(start, stop, n)` | $10^{\text{start}}$에서 $10^{\text{stop}}$까지 로그 간격 값 | `np.logspace(-1, 2, 4)` |

> **참고:** `np.matrix`는 더 이상 사용되지 않으며(deprecated), 항상 `np.array`를 대신 사용해야 한다. 2D 배열은 `@` 연산자를 통해 모든 행렬 연산을 지원한다.

```python
import numpy as np

# 1D array (vector)
data = np.array([12.2, 10.9, 13.6, 8.4, 11.1])
print(data.shape)  # (5,)

# 2D array (matrix)
A = np.array([[2, 4], [1, 3]])
print(A.shape)  # (2, 2)

# Special arrays
Z = np.zeros((5, 3))   # 5x3 matrix of zeros
O = np.ones((2, 3))    # 2x3 matrix of ones
I = np.eye(3)          # 3x3 identity matrix
```

### 2.3 인덱싱과 슬라이싱

NumPy는 **0 기반 인덱싱(Zero-Based Indexing)**을 사용한다. 2D 배열의 경우, 인덱스는 `[행, 열]`로 지정한다. 콜론 `:`은 전체 행 또는 열을 선택한다:

```python
A = np.array([[2, 4], [1, 3]])

print(A[0, 1])   # 4 (row 0, col 1)
print(A[0, :])   # [2, 4] (entire row 0)
print(A[:, 1])   # [4, 3] (entire column 1)
```

### 2.4 시퀀스 생성: arange, linspace, logspace

이 세 함수는 숫자의 시퀀스를 생성하지만, 간격을 지정하는 방식이 다르다:

**`np.arange`** — 스텝 크기로 지정, 끝점은 포함**하지 않음**:

```python
x = np.arange(0, 1, 0.1)  # [0.0, 0.1, 0.2, ..., 0.9] — 10 values, 1.0 NOT included
```

**`np.linspace`** — 점의 개수로 지정, 끝점은 기본적으로 **포함**:

```python
x = np.linspace(0, 1, 11)  # [0.0, 0.1, 0.2, ..., 1.0] — 11 values, 1.0 included
x = np.linspace(0, 1, 10, endpoint=False)  # [0.0, 0.1, ..., 0.9] — excludes 1.0
```

**`np.logspace`** — 로그 간격 값:

```python
# Default base=10: from 10^(-1) to 10^2
x = np.logspace(-1, 2, 4)  # [0.1, 1, 10, 100]

# Custom base: from 2^1 to 2^6
x = np.logspace(1, 6, 6, base=2.0)  # [2, 4, 8, 16, 32, 64]
```

---

<br>

## 3. 수학 연산

### 3.1 산술 연산자와 우선순위

Python은 다음과 같은 산술 연산자를 제공한다:

| 연산자 | 설명 | 예시 | 결과 |
|:---------|:-----------|:--------|:-------|
| `+` | 덧셈 | `3 + 2` | `5` |
| `-` | 뺄셈 | `3 - 2` | `1` |
| `*` | 곱셈 (배열은 요소별) | `3 * 2` | `6` |
| `/` | 나눗셈 | `7 / 2` | `3.5` |
| `//` | 정수 나눗셈(Floor Division) | `7 // 2` | `3` |
| `%` | 나머지(Modulus) | `7 % 2` | `1` |
| `**` | 거듭제곱 | `2 ** 3` | `8` |
| `@` | 행렬 곱셈 | `A @ B` | 행렬 곱 |

**연산자 우선순위**는 표준 수학 규칙(PEMDAS)을 따른다. 중요한 미묘한 점: **거듭제곱은 오른쪽에서 왼쪽으로 결합**한다:

```python
x, y, z = 2, 3, 2
print(x ** y ** z)  # 2 ** (3 ** 2) = 2 ** 9 = 512, NOT (2 ** 3) ** 2 = 64
```

### 3.2 요소별 곱셈 vs. 행렬 곱셈

이것은 NumPy에서 가장 중요한 구분 중 하나이다. `*` 연산자는 **요소별(Element-wise)** 곱셈을 수행하고, `@`는 **행렬 곱셈(Matrix Multiplication)**을 수행한다:

```python
import numpy as np

# 1D arrays
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])

print(a * b)  # [4, 10, 18] — element-wise: [1*4, 2*5, 3*6]
print(a @ b)  # 32 — dot product: 1*4 + 2*5 + 3*6

# 2D arrays (matrices)
A = np.array([[2, 4], [1, 3]])

print(A @ A)              # matrix multiplication
# [[8, 20],
#  [5, 13]]

print(np.multiply(A, A))  # element-wise multiplication
# [[4, 16],
#  [1,  9]]
```

### 3.3 내적, 외적, 전치

```python
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])

# Dot product: sum of element-wise products
print(np.dot(a, b))    # 32 (same as a @ b for 1D arrays)

# Cross product: vector perpendicular to both inputs
print(np.cross(a, b))  # [-3, 6, -3]

# Transpose: swap rows and columns
A = np.array([[2, 4], [1, 3]])
print(A.T)
# [[2, 1],
#  [4, 3]]
```

---

<br>

## 4. 내장 함수

### 4.1 주요 내장 함수

Python은 일상적인 계산에 필수적인 여러 내장 함수를 제공한다:

| 함수 | 설명 | 예시 |
|:---------|:-----------|:--------|
| `abs(x)` | 절댓값 (복소수의 경우 크기) | `abs(-6 + 4j)` → `7.211...` |
| `round(x, n)` | 소수점 이하 $n$자리로 반올림 | `round(3.14159, 2)` → `3.14` |
| `min(...)` | 최솟값 | `min(3, 1, 4)` → `1` |
| `max(...)` | 최댓값 | `max(3, 1, 4)` → `4` |
| `len(x)` | 시퀀스의 길이 | `len([1, 2, 3])` → `3` |
| `type(x)` | 객체의 타입 | `type(3.14)` → `<class 'float'>` |
| `float(x)` | 실수로 변환 | `float('3.14')` → `3.14` |
| `int(x)` | 정수로 변환 (절삭) | `int(3.7)` → `3` |
| `help(x)` | 문서 표시 | `help(abs)` |

```python
print(abs(-6 + 4j))  # 7.211102550927978 (magnitude: sqrt(36 + 16))
```

### 4.2 은행원 반올림

Python의 `round()` 함수는 **은행원 반올림(Banker's Rounding)**(또는 "짝수 반올림(Round Half to Even)")을 사용한다. 값이 정확히 두 정수의 중간일 때, **가장 가까운 짝수**로 반올림한다:

```python
print(round(4.5))     # 4 (rounds to nearest even)
print(round(4.51))    # 5 (not exactly halfway, rounds normally)
print(round(0.5))     # 0 (rounds to nearest even)
print(round(1.5))     # 2 (rounds to nearest even)
print(round(2.5))     # 2 (rounds to nearest even)
print(round(3.5))     # 4 (rounds to nearest even)
```

> **참고:** Python은 "은행원 반올림"(짝수 반올림)을 사용한다: `round(0.5) = 0`, `round(1.5) = 2`. 이는 통계 계산에서 누적 반올림 편향을 최소화한다.

### 4.3 응용: 번지 점프 속도

이차 항력을 가진 번지 점프의 해석해를 사용하여, 모든 시간 점에서 벡터화된 속도를 계산할 수 있다:

$$v(t) = \sqrt{\frac{mg}{c_d}} \tanh\!\left(\sqrt{\frac{c_d\, g}{m}}\, t\right)$$

```python
import math
import numpy as np

m, g, cd = 68.1, 9.81, 0.25
tm = np.linspace(0, 20, 11)  # time from 0 to 20 s, 11 points

v = math.sqrt(m * g / cd) * np.tanh(math.sqrt(cd * g / m) * tm)
print(v)
# [ 0.          9.63...  17.85...  23.53...  27.29...  29.67...
#  31.09...  31.87...  32.26...  32.44...  32.52...]
```

> **[Python]** `np.tanh()` 함수는 전체 `tm` 배열에 요소별로 작동하여, 명시적 루프 없이 동일한 형상의 결과 배열을 생성한다. 이것이 벡터화의 실제 동작이다.

---

<br>

## 5. Matplotlib을 이용한 시각화

### 5.1 기본 선 그래프

**Matplotlib**은 Python의 표준 그래프 라이브러리이다. `pyplot` 모듈은 MATLAB과 유사한 인터페이스를 제공한다:

```python
import matplotlib.pyplot as plt
import numpy as np
import math

m, g, cd = 68.1, 9.81, 0.25
tm = np.linspace(0, 20, 11)
v = math.sqrt(m * g / cd) * np.tanh(math.sqrt(cd * g / m) * tm)

plt.plot(tm, v)
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Bungee Jumper Velocity')
plt.grid(True)
plt.show()
```

### 5.2 범례가 있는 다중 곡선

같은 그래프에서 서로 다른 시나리오를 비교하려면, `plt.plot()`을 여러 번 호출하고 `label`과 `plt.legend()`를 사용한다:

```python
import matplotlib.pyplot as plt
import numpy as np
import math

m, g = 68.1, 9.81
tm = np.linspace(0, 20, 11)

cd1, cd2 = 0.25, 0.30
v1 = math.sqrt(m * g / cd1) * np.tanh(math.sqrt(cd1 * g / m) * tm)
v2 = math.sqrt(m * g / cd2) * np.tanh(math.sqrt(cd2 * g / m) * tm)

plt.plot(tm, v1, ls='-', marker='o', label=f'cd = {cd1}')
plt.plot(tm, v2, ls='--', marker='D', label=f'cd = {cd2}')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Bungee Jumper: Effect of Drag Coefficient')
plt.legend()
plt.grid(True)
plt.savefig('sample.png', bbox_inches='tight')
plt.show()
```

> **[Python]** `ls` 매개변수는 선 스타일을 제어하고(`'-'` 실선, `'--'` 파선, `':'` 점선), `marker`는 데이터 점 마커를 제어하며(`'o'` 원, `'D'` 다이아몬드, `'s'` 사각형), `label`은 범례 항목의 텍스트를 제공한다. `bbox_inches='tight'`는 저장 시 여분의 공백을 제거한다.

---

<br>

## 6. NumPy에서의 뷰와 복사본

**뷰(View)**와 **복사본(Copy)**의 차이를 이해하는 것은 과학 계산에서 매우 중요하다. 실수로 뷰를 수정하면 원본 데이터가 변경되어, 미묘한 버그가 발생할 수 있다.

### 6.1 뷰 (슬라이싱)

**기본 슬라이싱**(예: `x[1:3]`)을 사용하여 배열의 부분집합을 만들면, NumPy는 **뷰(View)**를 반환한다 — 동일한 기본 메모리를 가리키는 창이다. 뷰의 변경은 원본 배열에 영향을 준다:

```python
import numpy as np

x = np.array([1, 2, 3, 4, 5])
y = x[1:3]       # y is a VIEW of x — shares the same memory
print(y)          # [2, 3]

x[2] = 99         # modify the original
print(y)          # [2, 99] — y changed too, because it's a view!

# Verify that y is a view
print(y.base)     # [1, 2, 99, 4, 5] — shows the original (base) array
```

### 6.2 복사본 (팬시 인덱싱)

**팬시 인덱싱(Fancy Indexing)**(인덱스의 리스트나 배열로 인덱싱)을 사용하면, NumPy는 **복사본(Copy)**을 반환한다 — 완전히 독립적인 배열이다:

```python
import numpy as np

x = np.array([1, 2, 3, 4, 5])
y = x[[1, 2]]    # y is a COPY — independent memory
print(y)          # [2, 3]

x[2] = 99         # modify the original
print(y)          # [2, 3] — y is unchanged, because it's a copy

# Verify that y is a copy
print(y.base is None)  # True — no base array, it's independent
```

> **참고:** 이 구분은 과학 계산에서 매우 중요하다. 실수로 뷰를 수정하면 원본 데이터가 변경된다. 독립적인 복사본이 필요할 때는 `.copy()`를 명시적으로 사용하라: `y = x[1:3].copy()`.

---

<br>

## 요약

| 주제 | 핵심 함수 / 개념 |
|:------|:------------------------|
| 변수 | `type()`, 동적 타이핑, `math.pi`, `scipy.constants` |
| NumPy 배열 | `np.array`, `np.zeros`, `np.ones`, `np.eye`, `np.linspace` |
| 수학 연산 | `*` (요소별), `@` (행렬), `np.dot`, `np.cross` |
| 내장 함수 | `abs`, `round` (은행원 반올림), `min`, `max`, `help()` |
| 그래프 | `plt.plot`, `plt.xlabel`, `plt.legend`, `plt.savefig` |
| 뷰 vs. 복사본 | 슬라이싱 → 뷰 (공유 메모리), 팬시 인덱싱 → 복사본 (독립적) |

---
