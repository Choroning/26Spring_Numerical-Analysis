# 2장 실습 — Python 프로그래밍

[레포지토리 안내](../../README.ko.md) | [English](Concepts_Lab.md) | [한국어](Concepts_Lab.ko.md)

> **최종 수정일:** 2026-10-20
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 3

> **학습 목표:**
> 1. 수치 계산을 위한 구조화된 프로그램을 작성할 수 있다
> 2. 모듈화 프로그래밍 원칙을 적용할 수 있다
> 3. 예외 처리로 수치 계산 중에 생기는 오류를 감지하고 처리할 수 있다
> 4. LEGB 규칙에 따라 변수 이름이 어느 스코프에서 해석되는지 판단할 수 있다
> 5. 명시적 반복문을 NumPy 벡터화 연산으로 바꿀 수 있다

---

## 목차

- [1. 스크립트 파일과 함수](#1-스크립트-파일과-함수)
  - [1.1 번지 점프: 인라인 스크립트](#11-번지-점프-인라인-스크립트)
  - [1.2 함수로 리팩토링](#12-함수로-리팩토링)
  - [1.3 다중 반환값](#13-다중-반환값)
- [2. 변수 스코프: LEGB 규칙](#2-변수-스코프-legb-규칙)
  - [2.1 지역 스코프와 전역 스코프](#21-지역-스코프와-전역-스코프)
  - [2.2 둘러싸는 스코프 (중첩 함수)](#22-둘러싸는-스코프-중첩-함수)
  - [2.3 global 키워드](#23-global-키워드)
- [3. 인수](#3-인수)
  - [3.1 키워드 인수와 기본값](#31-키워드-인수와-기본값)
  - [3.2 가변 위치 인수 (*args)](#32-가변-위치-인수-args)
  - [3.3 가변 키워드 인수 (**kwargs)](#33-가변-키워드-인수-kwargs)
- [4. NumPy를 이용한 파일 입출력](#4-numpy를-이용한-파일-입출력)
  - [4.1 텍스트 파일 (CSV/TXT)](#41-텍스트-파일-csvtxt)
  - [4.2 바이너리 파일 (NPY/NPZ)](#42-바이너리-파일-npynpz)
- [5. 구조적 프로그래밍](#5-구조적-프로그래밍)
  - [5.1 조건문](#51-조건문)
  - [5.2 For 및 While 반복문](#52-for-및-while-반복문)
  - [5.3 벡터화 vs. 반복문](#53-벡터화-vs-반복문)
- [6. 일급 객체로서의 함수](#6-일급-객체로서의-함수)
  - [6.1 람다 함수](#61-람다-함수)
  - [6.2 클로저와 고차 함수](#62-클로저와-고차-함수)
  - [6.3 functools.partial](#63-functoolspartial)
  - [6.4 zip, filter, 컴프리헨션](#64-zip-filter-컴프리헨션)
- [7. 예외 처리](#7-예외-처리)
- [8. 데코레이터](#8-데코레이터)
- [9. OOP 기초](#9-oop-기초)
  - [9.1 `__init__`, `__repr__`, `__str__`을 가진 클래스](#91-__init__-__repr__-__str__을-가진-클래스)
  - [9.2 데이터클래스](#92-데이터클래스)
- [10. 타입 힌트](#10-타입-힌트)
- [11. 제너레이터](#11-제너레이터)
  - [11.1 메모리 효율성](#111-메모리-효율성)
  - [11.2 yield를 사용한 제너레이터 함수](#112-yield를-사용한-제너레이터-함수)
  - [11.3 대용량 데이터를 위한 배치 제너레이터](#113-대용량-데이터를-위한-배치-제너레이터)
- [요약](#요약)
- [점검 문제](#점검-문제)

---

<br>

## 1. 스크립트 파일과 함수

### 1.1 번지 점프: 인라인 스크립트

1장의 번지 점프 문제를 **이차 항력 모델(Quadratic Drag Model)** 로 다시 살펴본다. 지배 ODE는 다음과 같다:

```math
\frac{dv}{dt} = g - \frac{c_d}{m}v^2
```

여기서 $`c_d`$는 항력 계수(kg/m), $`m`$은 질량(kg), $`g`$는 중력 가속도(m/s$`^2`$)이다.

초기 조건 $`v(0) = 0`$에 대한 이 ODE의 **해석해(Analytical Solution)** 는 다음과 같다:

```math
v(t) = \sqrt{\frac{mg}{c_d}} \tanh\left(\sqrt{\frac{c_d g}{m}} \cdot t\right)
```

다음은 $`m = 68.1`$ kg, $`c_d = 0.25`$ kg/m인 점프 선수의 $`t = 4`$ s에서의 속도를 계산하는 간단한 인라인 스크립트이다:

```python
import math
import scipy.constants as pc

t, m, cd = 4, 68.1, 0.25
g = pc.g  # 9.80665 m/s^2 (standard gravity from scipy.constants)

v = math.sqrt(m * g / cd) * math.tanh(math.sqrt(cd * g / m) * t)
print(f'velocity = {v:.4f} m/s')  # velocity = 33.1022 m/s
```

> **참고:** `scipy.constants`는 완전한 정밀도의 물리 상수를 제공한다. `pc.g`는 표준 중력 가속도($`9.80665`$ m/s$`^2`$)를 제공하며, 이는 일반적으로 사용되는 근사값 $`9.81`$ m/s$`^2`$와 약간 다르다.

### 1.2 함수로 리팩토링

위의 인라인 스크립트는 재사용이 불가능하다. 계산을 **함수(Function)** 로 감싸면, 다른 매개변수로 반복적으로 호출할 수 있다:

```python
def freefall(t, m, cd):
    """Compute free-fall velocity of bungee jumper.

    Parameters:
        t: time (s)
        m: mass (kg)
        cd: drag coefficient (kg/m)
    Returns:
        velocity (m/s)
    """
    g = pc.g
    return math.sqrt(m * g / cd) * math.tanh(math.sqrt(cd * g / m) * t)
```

이제 이 함수를 어떤 매개변수 집합으로든 호출할 수 있다:

```python
print(freefall(4, 68.1, 0.25))    # 33.1022...
print(freefall(10, 80.0, 0.30))   # different jumper
print(freefall(0, 68.1, 0.25))    # 0.0 (at t=0, velocity is zero)
```

함수의 주요 이점은 다음과 같다:

1. **재사용성(Reusability)**: 다른 입력으로 여러 번 호출할 수 있다.
2. **가독성(Readability)**: 함수 이름이 코드가 무엇을 하는지 문서화한다.
3. **테스트 용이성(Testability)**: 알려진 값을 확인하여 정확성을 쉽게 검증할 수 있다.
4. **모듈성(Modularity)**: 복잡한 프로그램을 단순하고 잘 테스트된 함수들로 구축할 수 있다.

### 1.3 다중 반환값

Python 함수는 **튜플(Tuple)** 을 사용하여 여러 값을 반환할 수 있다. 이는 함수가 여러 관련 양을 계산하는 경우가 많은 과학 계산에서 매우 편리하다:

```python
def stats(x):
    """Return count, mean, and std of array x."""
    n = len(x)
    avg = sum(x) / n
    s = (sum((xi - avg) ** 2 for xi in x) / (n - 1)) ** 0.5
    return n, avg, s
```

호출자는 튜플을 별도의 변수로 언패킹한다:

```python
data = [2.3, 4.5, 1.2, 3.8, 5.1]
count, mean, std = stats(data)
print(f'n={count}, mean={mean:.2f}, std={std:.2f}')
# n=5, mean=3.38, std=1.61
```

> **참고:** 표준 편차 공식의 분모에 $`n-1`$을 사용하는 것(베셀 보정, Bessel's Correction)은 모집단 표준 편차가 아닌 **표본(Sample)** 표준 편차를 계산하기 때문이다. 실제 모집단 평균 대신 표본 평균을 사용하면 편차의 제곱합이 작게 나오는데, $`n-1`$로 나누면 분산 추정의 이 편향이 보정되어 표본 분산이 불편 추정량이 된다. 그 제곱근인 표본 표준 편차에는 작은 편향이 남는다.

---

<br>

## 2. 변수 스코프: LEGB 규칙

Python은 **LEGB 규칙**을 사용하여 변수 이름을 해석하며, 이는 스코프가 검색되는 순서를 정의한다:

| 스코프 | 설명 | 수명 |
|:------|:-----------|:---------|
| **L**ocal (지역) | 현재 함수 내에서 정의된 변수 | 함수 실행 중에만 존재 |
| **E**nclosing (둘러싸는) | 중첩 함수의 바깥(둘러싸는) 함수에 있는 변수 | 둘러싸는 함수가 실행되는 동안 존재 |
| **G**lobal (전역) | 모듈(파일) 수준에서 정의된 변수 | 전체 프로그램 동안 존재 |
| **B**uilt-in (내장) | Python의 내장 이름 (`print`, `len`, `range` 등) | 항상 사용 가능 |

*표 1. LEGB 규칙의 네 스코프*

Python이 변수 이름을 만나면, 이 스코프를 **순서대로** 검색한다: 지역 먼저, 그다음 둘러싸는, 전역, 내장 순이다. 처음 발견된 것을 사용하며, 네 스코프 어디에도 없으면 `NameError`가 발생한다.

```mermaid
flowchart LR
    L[Local] -->|없으면| E[Enclosing]
    E -->|없으면| G[Global]
    G -->|없으면| B[Built-in]
    B -->|없으면| X[NameError 발생]
```

*그림 1. LEGB 규칙의 이름 검색 순서*

### 2.1 지역 스코프와 전역 스코프

```python
x = 88  # global variable

def func():
    s = 10  # local to func
    print(x)  # can READ global x → 88
    print(s)  # local s → 10

func()
# print(s)  # NameError! s is local to func, not visible here
```

함수는 외부 스코프의 변수를 **읽을** 수 있지만, 명시적 선언 없이 그 이름에 값을 다시 대입하여 **수정할** 수는 없다. 함수 내에서 변수에 값을 대입하면, Python은 동일한 이름의 전역 변수가 존재하더라도 이를 **새로운 지역 변수**로 취급한다. 단, 전역 리스트에 `append`를 호출하는 것처럼 변수가 가리키는 가변 객체의 내용을 바꾸는 것은 대입이 아니므로 선언 없이도 가능하다:

```python
x = 88

def func():
    x = 10  # creates a NEW local x, does NOT modify the global x
    print(x)  # 10

func()
print(x)  # 88 — global x is unchanged
```

### 2.2 둘러싸는 스코프 (중첩 함수)

함수가 중첩될 때, 내부 함수는 둘러싸는(외부) 함수의 스코프에 있는 변수에 접근할 수 있다:

```python
def fun1():
    x = 'fun1_local'

    def fun2():
        print(x)  # accesses enclosing scope → 'fun1_local'

    fun2()

fun1()
# fun2()  # NameError — fun2 only exists inside fun1
```

둘러싸는 스코프는 **클로저(Closure)** 에서 중요하다. 클로저는 둘러싸는 함수가 실행을 마친 후에도 둘러싸는 스코프의 변수를 캡처하고 기억하는 함수이다.

### 2.3 global 키워드

함수 내에서 전역 변수를 수정하려면 `global` 키워드를 사용해야 한다:

```python
def modify_global():
    global x
    x = 99

x = 1
modify_global()
print(x)  # 99 — global x was modified
```

> **[프로그래밍언어]** **어휘적(정적) 스코핑(lexical scoping)** 은 이름이 가리키는 변수를 함수가 호출되는 위치가 아니라 소스 코드에서 함수가 정의된 위치로 정하는 규칙이며, LEGB 규칙은 Python에서 이 규칙을 구현한 것이다. 또한 Python은 함수 안의 대입을 지역 변수의 선언으로 해석하므로, 외부 스코프의 이름에 다시 대입하려면 `global`이나 `nonlocal`로 명시해야 한다.

마찬가지로, `nonlocal` 키워드는 둘러싸는(전역이 아닌) 스코프의 변수 수정을 허용한다:

```python
def outer():
    count = 0

    def increment():
        nonlocal count
        count += 1

    increment()
    increment()
    print(count)  # 2

outer()
```

> **참고:** `global`과 `nonlocal`의 과도한 사용은 코드를 이해하고 디버깅하기 어렵게 만든다. 값을 인수로 전달하고 결과를 반환하는 것이 좋다.

---

<br>

## 3. 인수

Python은 단순한 위치 매개변수를 넘어서는 유연한 인수 전달 메커니즘을 제공한다.

### 3.1 키워드 인수와 기본값

함수는 매개변수에 **기본값(Default Values)** 을 가질 수 있다. 기본값이 있는 매개변수는 함수 호출 시 선택 사항이다:

```python
def func_keyword(a, print_input=False):
    if print_input:
        print(f'Input: {a}')
    return a ** 2

# Call with just the required argument
print(func_keyword(5))                    # 25

# Call with the optional keyword argument
print(func_keyword(5, print_input=True))  # Input: 5 \n 25
```

키워드 인수는 모든 위치 인수 뒤에 오는 한, 순서에 관계없이 전달할 수 있다:

```python
def describe(name, age=25, city='Seoul'):
    return f'{name}, age {age}, from {city}'

print(describe('Alice', city='Sejong', age=22))
# Alice, age 22, from Sejong
```

### 3.2 가변 위치 인수 (*args)

`*args` 구문을 사용하면 함수가 임의 개수의 위치 인수를 받을 수 있다. 함수 내에서 `args`는 튜플이다:

```python
def mysum(*args):
    return sum(args)

print(mysum(1, 2, 3))       # 6
print(mysum(10, 20, 30, 40)) # 100
```

이는 미리 얼마나 많은 인수가 전달될지 알 수 없을 때 유용하다:

```python
def mean(*values):
    return sum(values) / len(values)

print(mean(1, 2, 3, 4, 5))  # 3.0
```

### 3.3 가변 키워드 인수 (**kwargs)

`**kwargs` 구문을 사용하면 함수가 임의 개수의 키워드 인수를 받을 수 있다. 함수 내에서 `kwargs`는 딕셔너리이다:

```python
def show_kwargs(**kwargs):
    for key, val in kwargs.items():
        print(f'{key} = {val}')

show_kwargs(name='Alice', age=30, major='CSE')
# name = Alice
# age = 30
# major = CSE
```

단일 함수 정의에서 모든 인수 타입을 결합할 수 있다. 순서는 위치 인수, `*args`, 기본값이 있는 키워드 인수, `**kwargs`여야 한다:

```python
def combined(a, b, *args, option=True, **kwargs):
    print(f'a={a}, b={b}')
    print(f'args={args}')
    print(f'option={option}')
    print(f'kwargs={kwargs}')

combined(1, 2, 3, 4, option=False, x=10, y=20)
# a=1, b=2
# args=(3, 4)
# option=False
# kwargs={'x': 10, 'y': 20}
```

---

<br>

## 4. NumPy를 이용한 파일 입출력

수치 계산에서는 파일에서 데이터를 읽고 결과를 저장하는 일이 빈번하다. NumPy는 이를 위한 효율적인 함수를 제공한다:

| 형식 | 저장 함수 | 불러오기 함수 | 사용 사례 |
|:-------|:-------------|:-------------|:---------|
| `.txt` / `.csv` | `np.savetxt` | `np.loadtxt` | 사람이 읽을 수 있는, 상호 운용 가능한 형식 |
| `.npy` (단일 배열) | `np.save` | `np.load` | 빠른 바이너리, 단일 배열 |
| `.npz` (다중 배열) | `np.savez` | `np.load` | 빠른 바이너리, 이름 있는 다중 배열 |

*표 2. 파일 형식별 NumPy 입출력 함수*

### 4.1 텍스트 파일 (CSV/TXT)

`np.loadtxt`는 텍스트 파일에서 열 형식의 수치 데이터를 읽는다. 헤더, 구분자, 열 선택을 처리한다:

```python
import numpy as np

# Load text data (e.g., KMA weather station data)
# File has 2 header rows; 3 columns: station_id, time, temperature
_, t, temp = np.loadtxt('./data/awsdata.txt', skiprows=2, unpack=True)
```

매개변수의 의미는 다음과 같다:

- `skiprows=2`: 처음 2행(헤더 행)을 건너뛴다.
- `unpack=True`: 데이터를 전치하여 각 열이 별도의 1D 배열이 되게 한다.
- `_` 변수는 버리려는 값(여기서는 관측소 ID 열)에 대한 Python 관례이다.

텍스트 파일로 데이터를 저장하려면 다음과 같이 한다:

```python
data = np.column_stack([t, temp])
np.savetxt('output.csv', data, delimiter=',', header='time,temperature', fmt='%.4f')
```

### 4.2 바이너리 파일 (NPY/NPZ)

바이너리 형식은 읽기/쓰기가 훨씬 빠르고 전체 부동 소수점 정밀도를 보존한다(서식 손실 없음):

```python
# Save a single array
np.save('data.npy', temp)
loaded = np.load('data.npy')

# Save multiple named arrays
np.savez('multi.npz', time=t, temperature=temp)
data = np.load('multi.npz')
print(data['time'])        # access by name
print(data['temperature'])
```

> **참고:** 바이너리 파일(`.npy`, `.npz`)은 사람이 읽을 수 없지만, 텍스트 파일보다 훨씬 빠르고 컴팩트하다. 다른 도구(Excel, R 등)와의 상호 운용이 필요할 때는 텍스트 형식을, Python 워크플로우 내의 중간 결과에는 바이너리 형식을 사용하라.

---

<br>

## 5. 구조적 프로그래밍

### 5.1 조건문

`if`/`elif`/`else` 구조는 분기 로직을 제공한다:

```python
def classify_temp(T):
    if T > 30:
        return 'Hot'
    elif T > 15:
        return 'Warm'
    else:
        return 'Cold'

print(classify_temp(35))  # Hot
print(classify_temp(20))  # Warm
print(classify_temp(5))   # Cold
```

### 5.2 For 및 While 반복문

**For 반복문**은 시퀀스를 순회한다:

```python
# Sum of integers 1 through 10
total = 0
for i in range(1, 11):
    total += i
print(total)  # 55
```

**`break`** 는 반복문을 즉시 종료하고, **`continue`** 는 남은 본문을 건너뛰고 다음 반복으로 넘어간다:

```python
# Find the first multiple of 7 greater than 0
for i in range(100):
    if i % 7 == 0 and i > 0:
        print(f'First multiple of 7: {i}')  # 7
        break
```

**While 반복문**은 조건이 `False`가 될 때까지 반복한다:

```python
count = 10
while count > 0:
    print(count, end=' ')
    count -= 1
# 10 9 8 7 6 5 4 3 2 1
```

조건이 만족될 때까지 반복하는 구조는 수치 해법에서 오차가 허용 범위 아래로 떨어질 때까지 반복하는 **반복 수렴(Iterative Convergence)** 에 특히 유용하다. 아래 예는 수렴하지 않을 때 무한 반복을 막기 위해 `while` 대신 최대 반복 횟수를 둔 `for` 반복문과 `break`를 사용한다:

```python
# Example: iterate until convergence
x = 1.0
tol = 1e-8
for _ in range(1000):  # safety limit
    x_new = 0.5 * (x + 2 / x)  # Babylonian method for sqrt(2)
    if abs(x_new - x) < tol:
        break
    x = x_new
print(f'sqrt(2) ≈ {x_new:.10f}')  # 1.4142135624
```

### 5.3 벡터화 vs. 반복문

이것은 과학 Python에서 가장 중요한 성능 개념이다. **벡터화(Vectorization)** 란 명시적인 Python 루프를 컴파일된 C 코드로 실행되는 NumPy 배열 연산으로 대체하는 것을 의미한다:

```python
import numpy as np
import time

t = np.linspace(0, 2 * np.pi, 1000000)

# SLOW: explicit Python for loop
start = time.time()
result_loop = np.zeros_like(t)
for i in range(len(t)):
    result_loop[i] = np.cos(10 * t[i])
loop_time = time.time() - start
print(f'Loop: {loop_time:.4f}s')

# FAST: vectorized NumPy operation
start = time.time()
result_vec = np.cos(10 * t)
vec_time = time.time() - start
print(f'Vectorized: {vec_time:.4f}s')

print(f'Speedup: {loop_time / vec_time:.0f}x')
# Vectorized is typically 10-100x faster!
```

> **[컴퓨터프로그래밍]** C 같은 컴파일 언어는 소스 코드를 실행 전에 기계어로 번역하므로 반복문의 한 반복이 몇 개의 기계 명령으로 끝나지만, 인터프리터 언어인 Python은 실행 중에 문장마다 해석하고 값의 타입을 확인한다. 벡터화는 반복을 NumPy 내부의 컴파일된 C 루프로 옮겨 이 해석 비용을 배열 전체에 대해 한 번만 치르게 하므로, 가능하면 명시적 루프 대신 NumPy 벡터화 연산을 쓴다.

루프가 왜 그렇게 느린가? Python은 **인터프리터 언어(Interpreted Language)** 이므로, 루프의 각 반복은 타입 검사, 함수 디스패치, 기타 오버헤드를 수반한다. NumPy의 벡터화 연산은 전체 배열에 대해 하나의 최적화된 C 루틴을 실행하여 이를 우회한다.

> **참고:** **경험 규칙.** 배열 요소를 하나씩 독립적으로 처리하는 `for` 반복문을 작성하고 있다면, 거의 확실히 훨씬 더 빠른 벡터화된 NumPy 대안이 존재한다. 단, 오일러 방법처럼 각 단계가 앞 단계의 결과에 의존하는 계산은 그대로 벡터화할 수 없다.

---

<br>

## 6. 일급 객체로서의 함수

Python에서 함수는 **일급 객체(First-Class Objects)** 이다. 즉, 함수를 변수에 대입하고, 인수로 전달하고, 다른 함수에서 반환하고, 데이터 구조에 저장할 수 있다. 이는 강력한 함수형 프로그래밍 패턴을 가능하게 한다.

### 6.1 람다 함수

**람다 함수(Lambda Functions)** 는 익명의 단일 표현식 함수이다:

```python
f1 = lambda x, y: x ** 2 + y ** 2
print(f1(3, 4))  # 25
```

람다는 다른 함수에 직접 전달되는 짧은 일회용 함수로 가장 유용하다:

```python
data = [(1, 'b'), (3, 'a'), (2, 'c')]
sorted_data = sorted(data, key=lambda item: item[1])
print(sorted_data)  # [(3, 'a'), (1, 'b'), (2, 'c')]
```

### 6.2 클로저와 고차 함수

**클로저(Closure)** 는 둘러싸는 스코프의 변수를 캡처하는 함수이다. **고차 함수(Higher-Order Function)** 는 다른 함수를 인수로 받거나 함수를 반환하는 함수이다:

```python
# Return a function (closure)
def make_power(n):
    return lambda x: x ** n

square = make_power(2)
cube = make_power(3)
print(square(5))  # 25
print(cube(5))    # 125
```

여기서 `make_power`는 둘러싸는 스코프에서 `n`의 값을 기억하는 람다를 반환하며, `make_power`의 실행이 끝난 후에도 그렇다. 이것이 클로저의 본질이다.

함수를 인수로 전달하는 예는 다음과 같다:

```python
def f_at_midpoint(func, a, b):
    """Evaluate func at the midpoint of [a, b]."""
    mid = (a + b) / 2
    return func(mid)

import numpy as np
print(f_at_midpoint(np.cos, 0, 10))  # cos(5) = 0.2837...
print(f_at_midpoint(np.sin, 0, 10))  # sin(5) = -0.9589...
```

이 패턴은 수치 해법의 기본이다. 근 찾기, 적분기, 최적화기처럼 많은 알고리즘이 함수를 입력으로 받는다.

### 6.3 functools.partial

`functools.partial`은 기존 함수의 일부 인수를 고정하여 새로운 함수를 만든다:

```python
from functools import partial

def power(base, exp):
    return base ** exp

square = partial(power, exp=2)
cube = partial(power, exp=3)

print(square(5))  # 25
print(cube(5))    # 125
```

이는 함수의 시그니처를 다른 함수가 기대하는 것과 맞추어야 할 때 유용하다. 예를 들어, ODE 풀이기가 `f(t, y)`를 기대하지만 함수에 추가 매개변수가 있는 경우, `partial`로 그 추가 매개변수를 고정할 수 있다.

### 6.4 zip, filter, 컴프리헨션

**`zip`** 은 여러 이터러블의 요소를 쌍으로 묶는다:

```python
names = ['Alice', 'Bob']
scores = [95, 87]

for name, score in zip(names, scores):
    print(f'{name}: {score}')
# Alice: 95
# Bob: 87
```

**`filter`** 는 조건을 만족하는 요소를 선택한다:

```python
def even(x):
    return x % 2 == 0

print(list(filter(even, range(10))))  # [0, 2, 4, 6, 8]
```

**컴프리헨션(Comprehensions)** 은 `map`과 `filter`의 간결하고 파이썬다운 대안이다:

```python
# List comprehension (replaces map + filter)
squares = [x ** 2 for x in range(10)]
# [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

evens = [x for x in range(10) if x % 2 == 0]
# [0, 2, 4, 6, 8]

# Dictionary comprehension
d = {name: score for name, score in zip(names, scores)}
# {'Alice': 95, 'Bob': 87}

# Set comprehension
unique_remainders = {x % 3 for x in range(10)}
# {0, 1, 2}
```

> **참고:** 현대 Python에서는 `map`과 `filter`보다 컴프리헨션이 일반적으로 선호되는데, 가독성이 더 높기 때문이다. 다만 `map(str, data)`처럼 내장 함수를 적용할 때는 `map`이 약간 더 빠를 수 있고, 람다를 넘기면 함수 호출 비용 때문에 대체로 컴프리헨션보다 느리다.

---

<br>

## 7. 예외 처리

실세계 데이터와 계산은 오류를 발생시킬 수 있다. Python의 `try`/`except` 메커니즘을 사용하면 프로그램을 중단시키지 않고 오류를 처리할 수 있다:

```python
def safe_divide(a, b):
    try:
        return a / b
    except ZeroDivisionError:
        print('Division by zero!')
        return None
    except TypeError:
        print('Invalid types!')
        return None

print(safe_divide(10, 3))    # 3.333...
print(safe_divide(10, 0))    # Division by zero! → None
print(safe_divide('a', 2))   # Invalid types! → None
```

실용적인 응용으로, 비수치 값이 포함될 수 있는 노이즈 데이터를 정제하는 예는 다음과 같다:

```python
raw = ['1.5', '2.3', 'N/A', '4.1', '-']
clean = []
for val in raw:
    try:
        clean.append(float(val))
    except ValueError:
        pass  # skip non-numeric values

print(clean)  # [1.5, 2.3, 4.1]
```

`try`/`except` 구조에는 `else`(예외가 발생하지 않으면 실행)와 `finally`(항상 실행, 정리에 유용)도 포함할 수 있다:

```python
try:
    result = 10 / 2
except ZeroDivisionError:
    print('Error!')
else:
    print(f'Success: {result}')  # runs because no exception
finally:
    print('Done.')               # always runs
```

---

<br>

## 8. 데코레이터

**데코레이터(Decorator)** 는 다른 함수를 받아 동작을 확장하고 수정된 함수를 반환하는 함수이다. 데코레이터는 `@` 구문을 사용한다:

```python
def simple_decorator(func):
    def wrapper(*args, **kwargs):
        print(f'Calling {func.__name__}...')
        result = func(*args, **kwargs)
        print(f'{func.__name__} finished.')
        return result
    return wrapper

@simple_decorator
def greet(name):
    print(f'Hello, {name}!')

greet('Alice')
# Output:
# Calling greet...
# Hello, Alice!
# greet finished.
```

`greet` 위의 `@simple_decorator` 구문은 다음과 동일하다:

```python
greet = simple_decorator(greet)
```

> **[프로그래밍언어]** **문법적 설탕(syntactic sugar)** 은 언어가 표현할 수 있는 것을 늘리지 않고 기존 구문을 더 짧게 쓰게 해 주는 구문이다. 데코레이터는 고차 함수에 함수를 넘기고 그 결과로 이름을 다시 묶는 `func = decorator(func)`를 함수 정의 위의 `@decorator` 한 줄로 쓰게 해 주므로, 웹 프레임워크(Flask 라우트), 테스팅, 로깅처럼 여러 함수에 같은 동작을 덧붙이는 곳에서 널리 사용된다.

실용적인 예시로, 실행 시간을 측정하는 데코레이터는 다음과 같다:

```python
import time

def timer(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        elapsed = time.time() - start
        print(f'{func.__name__} took {elapsed:.4f}s')
        return result
    return wrapper

@timer
def slow_function():
    time.sleep(1)
    return 'done'

slow_function()  # slow_function took 1.00XXs
```

---

<br>

## 9. OOP 기초

### 9.1 `__init__`, `__repr__`, `__str__`을 가진 클래스

객체 지향 프로그래밍(OOP)은 데이터(속성)와 동작(메서드)을 묶는 **객체(Object)** 중심으로 코드를 구성한다:

```python
class MyStudent:
    def __init__(self, name, scores):
        self.name = name
        self.scores = scores

    def __repr__(self):
        return f'MyStudent({self.name}, {self.scores})'

    def __str__(self):
        return f'{self.name}: avg={sum(self.scores)/len(self.scores):.1f}'

s = MyStudent('Alice', [95, 88, 92])
print(repr(s))  # MyStudent(Alice, [95, 88, 92])
print(str(s))   # Alice: avg=91.7
```

주요 특수 메서드(던더 메서드)는 다음과 같다:

| 메서드 | 목적 | 호출 시점 |
|:-------|:--------|:---------|
| `__init__` | 객체 초기화 | `MyStudent(...)` |
| `__repr__` | 명확한 문자열 표현 (개발자용) | `repr(obj)`, REPL |
| `__str__` | 사람이 읽기 쉬운 문자열 표현 | `print(obj)`, `str(obj)` |

*표 3. 주요 특수 메서드와 호출 시점*

### 9.2 데이터클래스

`@dataclass` 데코레이터(Python 3.7+)는 `__init__`, `__repr__`, 동등 비교 메서드 `__eq__`를 자동으로 생성하여 보일러플레이트를 제거한다. 대소 비교 메서드는 `@dataclass(order=True)`로 지정해야 생성된다:

```python
from dataclasses import dataclass

@dataclass
class Student:
    name: str
    scores: list

    @property
    def average(self):
        return sum(self.scores) / len(self.scores)

    @property
    def grade(self):
        avg = self.average
        if avg >= 90:
            return 'A'
        elif avg >= 80:
            return 'B'
        elif avg >= 70:
            return 'C'
        else:
            return 'F'
```

`@property` 데코레이터는 메서드를 읽기 전용 속성으로 변환하여, 괄호 없이 접근할 수 있게 한다:

```python
students = [
    Student('Alice', [95, 88, 92]),
    Student('Bob', [72, 85, 68]),
]

for s in students:
    print(f'{s.name}: avg={s.average:.1f}, grade={s.grade}')
# Alice: avg=91.7, grade=A
# Bob: avg=75.0, grade=C

# Sort by average (descending)
ranked = sorted(students, key=lambda s: s.average, reverse=True)
print(ranked)
# [Student(name='Alice', scores=[95, 88, 92]),
#  Student(name='Bob', scores=[72, 85, 68])]
```

> **참고:** 데이터클래스는 현대 Python에서 단순한 데이터 보유 클래스를 만드는 데 선호되는 방법이다. 수동으로 `__init__`을 작성하는 것에 비해 보일러플레이트를 줄이고 코드를 더 깔끔하고 오류가 적게 만든다.

---

<br>

## 10. 타입 힌트

**타입 힌트(Type Hints)**(PEP 484)는 함수 매개변수와 반환값의 예상 타입을 문서화한다. 실행 시점에 타입을 강제하지는 않지만 더 나은 IDE 지원, 문서화, 정적 분석을 가능하게 한다:

```python
from typing import List, Optional, Tuple

def train_test_split_idx(
    n: int,
    train_ratio: float = 0.8,
    seed: Optional[int] = None
) -> Tuple[List[int], List[int]]:
    """Split indices into train and test sets.

    Parameters:
        n: total number of samples
        train_ratio: fraction of data for training (default 0.8)
        seed: random seed for reproducibility (None = no seed)

    Returns:
        Tuple of (train_indices, test_indices)
    """
    import random
    if seed is not None:
        random.seed(seed)
    indices = list(range(n))
    random.shuffle(indices)
    split = int(n * train_ratio)
    return indices[:split], indices[split:]
```

사용법은 다음과 같다:

```python
train_idx, test_idx = train_test_split_idx(100, train_ratio=0.7, seed=42)
print(f'Train: {len(train_idx)} samples, Test: {len(test_idx)} samples')
# Train: 70 samples, Test: 30 samples
```

주요 타입 힌트는 다음과 같다:

| 힌트 | 의미 |
|:-----|:--------|
| `int`, `float`, `str` | 기본 타입 |
| `List[int]` | 정수 리스트 |
| `Tuple[int, float]` | 특정 요소 타입의 튜플 |
| `Optional[int]` | `int` 또는 `None` |
| `Dict[str, float]` | 문자열 키와 실수 값의 딕셔너리 |
| `Callable[[int], float]` | `int`를 받아 `float`를 반환하는 함수 |

*표 4. 주요 타입 힌트와 의미*

> **참고:** Python 3.9+부터는 `typing`에서 임포트하는 대신 내장 타입을 직접 사용할 수 있다(`list[int]`, `tuple[int, float]`). Python 3.10+부터는 `X | None`이 `Optional[X]`를 대체한다.

---

<br>

## 11. 제너레이터

### 11.1 메모리 효율성

제너레이터(Generator)는 모든 값을 메모리에 한꺼번에 저장하지 않고, 필요할 때 하나씩 **지연 평가(Lazy Evaluation)** 방식으로 값을 생성한다. 이는 대용량 데이터셋을 다룰 때 매우 중요하다:

```python
import sys

# List comprehension: stores ALL values in memory
list_comp = [x ** 2 for x in range(10000)]

# Generator expression: stores only the formula, computes on demand
gen_expr = (x ** 2 for x in range(10000))

print(sys.getsizeof(list_comp))  # ~85,000 bytes
print(sys.getsizeof(gen_expr))   # ~200 bytes (!)
```

차이는 극적이다: 리스트는 10,000개의 정수에 대한 참조를 저장하고, 제너레이터는 다음 값을 생성하는 데 필요한 상태만 저장한다. `sys.getsizeof`는 리스트가 가진 참조 배열의 크기만 재며, 정수 객체 자체가 차지하는 메모리는 포함하지 않는다. 구체적인 바이트 수는 Python 버전에 따라 조금씩 다르다.

### 11.2 yield를 사용한 제너레이터 함수

`yield` 키워드는 일반 함수를 **제너레이터 함수(Generator Function)** 로 변환한다. `yield`를 만날 때마다 함수는 실행을 일시 중지하고 값을 반환한다. 다음 반복에서 중단된 곳에서 다시 시작한다:

```python
def countdown(n):
    while n > 0:
        yield n
        n -= 1

for val in countdown(5):
    print(val, end=' ')  # 5 4 3 2 1
```

실행 흐름은 다음과 같다:
1. `countdown(5)`를 호출하면 제너레이터 객체가 반환되며, 본문은 아직 실행되지 않는다.
2. 첫 번째 `next()`는 `yield 5`까지 실행한 뒤 일시 중지하고 `5`를 반환한다.
3. 두 번째 `next()`는 `yield` 다음부터 재개하여 `n`을 `4`로 줄이고, `yield 4`에 도달하면 `4`를 반환한다.
4. 같은 과정이 `while` 조건이 `False`가 될 때까지 계속된다.
5. 함수가 (암묵적으로) 반환하면 `StopIteration`이 발생하여 루프가 종료된다.

### 11.3 대용량 데이터를 위한 배치 제너레이터

제너레이터는 대용량 데이터셋을 배치 단위로 처리할 때 특히 유용하다(기계학습과 데이터 분석에서 흔히 쓰인다):

```python
import numpy as np

def batch_generator(data, batch_size):
    """Yield successive batches from data."""
    for i in range(0, len(data), batch_size):
        yield data[i:i + batch_size]

data = np.arange(100)
for batch in batch_generator(data, 30):
    print(f'Batch size: {len(batch)}')
# Batch size: 30
# Batch size: 30
# Batch size: 30
# Batch size: 10
```

> **[프로그래밍언어]** **지연 평가(lazy evaluation)** 는 식의 값을 실제로 필요해질 때까지 계산하지 않는 평가 전략이다. 제너레이터는 이 전략을 반복에 적용하여 `yield`에서 실행을 멈추고 다음 값에 필요한 상태만 유지하므로, 사용 가능한 RAM보다 큰 데이터셋도 한 부분씩 처리할 수 있다.

또 다른 실용적인 예시로, 대용량 파일 전체를 메모리에 로드하지 않고 한 줄씩 읽는 제너레이터는 다음과 같다:

```python
def read_large_file(filepath):
    """Yield lines from a large file one at a time."""
    with open(filepath, 'r') as f:
        for line in f:
            yield line.strip()

# Process a 10GB file without running out of memory
# for line in read_large_file('huge_data.txt'):
#     process(line)
```

---

<br>

## 요약

| 개념 | 핵심 요약 |
|:------|:------------|
| 함수 | `def`로 정의하고 독스트링(Docstrings)으로 문서화하며, 다중 반환과 `*args`/`**kwargs`를 지원한다. |
| 스코프 | 이름은 LEGB 규칙으로 해석되며, `global`과 `nonlocal`로 바깥 스코프의 이름에 대입하고 클로저는 바깥 변수를 기억한다. |
| 파일 입출력 | 텍스트 파일은 `np.loadtxt`로, 바이너리 파일은 `np.save`/`np.load`와 `np.savez`로 다룬다. |
| 제어 흐름 | `if`/`elif`/`else`, `for`, `while`, `break`/`continue`로 분기와 반복을 제어한다. |
| 벡터화 | NumPy 벡터화 연산은 Python 루프보다 일반적으로 10~100배 빠르다. |
| 함수형 | `lambda`, `partial`, 고차 함수, 컴프리헨션으로 함수를 값처럼 다룬다. |
| 예외 | `try`/`except`/`else`/`finally`로 오류를 처리한다. |
| 데코레이터 | `@decorator`는 `func = decorator(func)`의 문법적 설탕이다. |
| OOP | `__init__`으로 객체를 초기화하고, `@dataclass`와 `@property`로 보일러플레이트를 줄인다. |
| 타입 힌트 | `int`, `List[int]`, `Optional[...]`, `Tuple[...]`로 예상 타입을 문서화하며, 실행 시점에 강제되지는 않는다. |
| 제너레이터 | `yield`를 사용한 지연 평가로 메모리 효율적인 반복을 제공한다. |

---

<br>

## 점검 문제

1. **LEGB 규칙:** 전역에 `x = 'global'`이 있고, 함수 `outer` 안에 `x = 'enclosing'`과 중첩 함수 `inner`가 있다. `inner`가 `print(x)`만 실행할 때 `outer()`를 호출하면 무엇이 출력되는가?

   <details>
   <summary>정답</summary>

   `'enclosing'`이 출력된다. `inner`에는 지역 변수 `x`가 없으므로 다음 순서인 둘러싸는 스코프(`outer`)에서 `x`를 찾고, 거기서 찾았으므로 전역 스코프까지 가지 않는다.

   </details>

2. **지역 변수와 전역 변수:** 전역에 `count = 0`이 있고 함수 `f`가 `global` 선언 없이 `count = 5`를 실행한다. `f()`를 호출한 뒤 `print(count)`의 출력은 무엇인가?

   <details>
   <summary>정답</summary>

   `0`이 출력된다. 함수 안에서 대입한 `count`는 새로운 지역 변수가 되어 전역 `count`를 바꾸지 않는다. 전역 변수를 바꾸려면 `global count`를 선언해야 한다.

   </details>

3. **다중 반환값과 베셀 보정:** 1.3절의 `stats` 함수에 `[2, 4, 6]`을 전달하면 반환되는 튜플은 무엇인가?

   <details>
   <summary>정답</summary>

   `(3, 4.0, 2.0)`이다. 평균은 12/3 = 4.0이고, 편차의 제곱합은 4 + 0 + 4 = 8이다. 이를 $`n - 1 = 2`$로 나누면 4이고, 그 제곱근인 표본 표준 편차는 2.0이다.

   </details>

4. **인수의 결합:** 3.3절의 `combined` 함수를 `combined(1, 2, 3, option=False, z=0)`로 호출하면 `args`, `option`, `kwargs`는 각각 무엇인가?

   <details>
   <summary>정답</summary>

   `a=1`, `b=2`를 채우고 남은 위치 인수 3이 `args`에 들어가므로 `args=(3,)`이다. `option=False`이고, 매개변수 이름과 맞지 않는 키워드 인수 `z=0`은 `kwargs={'z': 0}`에 들어간다.

   </details>

5. **배치 제너레이터:** 11.3절의 `batch_generator(np.arange(10), 4)`가 생성하는 배치의 크기를 순서대로 쓰라.

   <details>
   <summary>정답</summary>

   4, 4, 2이다. 시작 인덱스가 0, 4, 8이므로 `data[0:4]`, `data[4:8]`, `data[8:12]`가 생성되며, 마지막 슬라이스는 배열 끝에서 잘려 원소가 2개이다.

   </details>
