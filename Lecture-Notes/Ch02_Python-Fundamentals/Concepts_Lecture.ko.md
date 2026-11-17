# 2장 — Python 기초

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **최종 수정일:** 2026-11-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 2

> **학습 목표:**
> 1. Python과 NumPy, SciPy, Matplotlib으로 이루어진 수치 계산 환경의 구성을 설명할 수 있다
> 2. NumPy 배열로 요소별 연산과 행렬 연산을 구분하여 수행할 수 있다
> 3. NumPy와 SciPy 함수로 배열을 만들고 연립방정식을 풀며, Matplotlib으로 그래프를 그릴 수 있다

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

이 행렬의 행렬식은 정확히 $`2 \times 3 - 1 \times 5 = 1`$이지만, `np.linalg.det`는 LU 분해를 거쳐 부동 소수점으로 계산하므로 1.0000000000000002처럼 마지막 자리에 반올림 오차가 나타난다. 반올림 오차는 4장에서 다룬다.

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
