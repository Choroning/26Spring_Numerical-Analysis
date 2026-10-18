# Chapter 2 Lab — Python Programming

[Repository Guide](../../README.md) | [English](Concepts_Lab.md) | [한국어](Concepts_Lab.ko.md)

> **Last Updated:** 2026-10-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 2

> **Learning Objectives:**
> 1. Write structured programs for numerical computation
> 2. Apply modular programming principles
> 3. Debug and test numerical algorithms
> 4. Determine in which scope a variable name is resolved under the LEGB rule
> 5. Replace explicit loops with NumPy vectorized operations

---

## Table of Contents

- [1. Script Files and Functions](#1-script-files-and-functions)
  - [1.1 Bungee Jumper: Inline Script](#11-bungee-jumper-inline-script)
  - [1.2 Refactoring as a Function](#12-refactoring-as-a-function)
  - [1.3 Multiple Return Values](#13-multiple-return-values)
- [2. Variable Scope: LEGB Rule](#2-variable-scope-legb-rule)
  - [2.1 Local and Global Scope](#21-local-and-global-scope)
  - [2.2 Enclosing Scope (Nested Functions)](#22-enclosing-scope-nested-functions)
  - [2.3 The global Keyword](#23-the-global-keyword)
- [3. Arguments](#3-arguments)
  - [3.1 Keyword Arguments and Defaults](#31-keyword-arguments-and-defaults)
  - [3.2 Variable Positional Arguments (*args)](#32-variable-positional-arguments-args)
  - [3.3 Variable Keyword Arguments (**kwargs)](#33-variable-keyword-arguments-kwargs)
- [4. File I/O with NumPy](#4-file-io-with-numpy)
  - [4.1 Text Files (CSV/TXT)](#41-text-files-csvtxt)
  - [4.2 Binary Files (NPY/NPZ)](#42-binary-files-npynpz)
- [5. Structured Programming](#5-structured-programming)
  - [5.1 Conditional Statements](#51-conditional-statements)
  - [5.2 For and While Loops](#52-for-and-while-loops)
  - [5.3 Vectorization vs. Loops](#53-vectorization-vs-loops)
- [6. Functions as First-Class Objects](#6-functions-as-first-class-objects)
  - [6.1 Lambda Functions](#61-lambda-functions)
  - [6.2 Closures and Higher-Order Functions](#62-closures-and-higher-order-functions)
  - [6.3 functools.partial](#63-functoolspartial)
  - [6.4 zip, filter, and Comprehensions](#64-zip-filter-and-comprehensions)
- [7. Exception Handling](#7-exception-handling)
- [8. Decorators](#8-decorators)
- [9. OOP Basics](#9-oop-basics)
  - [9.1 Classes with `__init__`, `__repr__`, `__str__`](#91-classes-with-__init__-__repr__-__str__)
  - [9.2 Dataclasses](#92-dataclasses)
- [10. Type Hints](#10-type-hints)
- [11. Generators](#11-generators)
  - [11.1 Memory Efficiency](#111-memory-efficiency)
  - [11.2 Generator Functions with yield](#112-generator-functions-with-yield)
  - [11.3 Batch Generator for Large Data](#113-batch-generator-for-large-data)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Script Files and Functions

### 1.1 Bungee Jumper: Inline Script

We revisit the bungee jumper problem from Chapter 1, now with a **quadratic drag model**. The governing ODE is:

```math
\frac{dv}{dt} = g - \frac{c_d}{m}v^2
```

where $`c_d`$ is the drag coefficient (kg/m), $`m`$ is mass (kg), and $`g`$ is gravitational acceleration (m/s$`^2`$).

The **analytical solution** for this ODE with $`v(0) = 0`$ is:

```math
v(t) = \sqrt{\frac{mg}{c_d}} \tanh\left(\sqrt{\frac{c_d g}{m}} \cdot t\right)
```

Here is a simple inline script that computes the velocity at $`t = 4`$ s for a jumper with $`m = 68.1`$ kg and $`c_d = 0.25`$ kg/m:

```python
import math
import scipy.constants as pc

t, m, cd = 4, 68.1, 0.25
g = pc.g  # 9.80665 m/s^2 (standard gravity from scipy.constants)

v = math.sqrt(m * g / cd) * math.tanh(math.sqrt(cd * g / m) * t)
print(f'velocity = {v:.4f} m/s')  # velocity = 33.1022 m/s
```

> **[Computer Programming]** `scipy.constants` provides physical constants with full precision. `pc.g` gives the standard acceleration due to gravity ($`9.80665`$ m/s$`^2`$), which is slightly different from the commonly used approximation $`9.81`$ m/s$`^2`$.

### 1.2 Refactoring as a Function

The inline script above is not reusable. By wrapping the computation in a **function**, we can call it repeatedly with different parameters:

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

Now the function can be called with any set of parameters:

```python
print(freefall(4, 68.1, 0.25))    # 33.1022...
print(freefall(10, 80.0, 0.30))   # different jumper
print(freefall(0, 68.1, 0.25))    # 0.0 (at t=0, velocity is zero)
```

The key benefits of functions are:

1. **Reusability**: the function can be called many times with different inputs.
2. **Readability**: the function name documents what the code does.
3. **Testability**: correctness is easy to verify by checking known values.
4. **Modularity**: complex programs are built from simple, well-tested functions.

### 1.3 Multiple Return Values

Python functions can return multiple values using **tuples**. This is extremely convenient for scientific computing, where a function often computes several related quantities:

```python
def stats(x):
    """Return count, mean, and std of array x."""
    n = len(x)
    avg = sum(x) / n
    s = (sum((xi - avg) ** 2 for xi in x) / (n - 1)) ** 0.5
    return n, avg, s
```

The caller unpacks the tuple into separate variables:

```python
data = [2.3, 4.5, 1.2, 3.8, 5.1]
count, mean, std = stats(data)
print(f'n={count}, mean={mean:.2f}, std={std:.2f}')
# n=5, mean=3.38, std=1.61
```

> **Note:** The standard deviation formula uses $`n-1`$ in the denominator (Bessel's correction) because we are computing the **sample** standard deviation, not the population standard deviation. Using the sample mean instead of the true population mean makes the sum of squared deviations smaller, and dividing by $`n-1`$ corrects this bias in the variance estimate so that the sample variance becomes an unbiased estimator. Its square root, the sample standard deviation, still carries a small bias.

---

<br>

## 2. Variable Scope: LEGB Rule

Python resolves variable names using the **LEGB rule**, which defines the order in which scopes are searched:

| Scope | Description | Lifetime |
|:------|:-----------|:---------|
| **L**ocal | Variables defined inside the current function | Exists only during function execution |
| **E**nclosing | Variables in the enclosing (outer) function for nested functions | Exists while the enclosing function runs |
| **G**lobal | Variables defined at the module (file) level | Exists for the entire program |
| **B**uilt-in | Python's built-in names (`print`, `len`, `range`, etc.) | Always available |

*Table 1. The four scopes of the LEGB rule*

When Python encounters a variable name, it searches these scopes **in order**: Local first, then Enclosing, then Global, then Built-in. It uses the first match found, and raises a `NameError` if the name is in none of the four scopes.

```mermaid
flowchart LR
    L[Local] -->|not found| E[Enclosing]
    E -->|not found| G[Global]
    G -->|not found| B[Built-in]
    B -->|not found| X[NameError raised]
```

*Figure 1. Name lookup order under the LEGB rule*

### 2.1 Local and Global Scope

```python
x = 88  # global variable

def func():
    s = 10  # local to func
    print(x)  # can READ global x → 88
    print(s)  # local s → 10

func()
# print(s)  # NameError! s is local to func, not visible here
```

A function can **read** variables from outer scopes, but it **cannot modify** them by reassigning the name without explicit declaration. If you assign to a variable inside a function, Python treats it as a **new local variable**, even if a global variable with the same name exists. Changing the contents of a mutable object that the variable refers to, such as calling `append` on a global list, is not an assignment and works without any declaration:

```python
x = 88

def func():
    x = 10  # creates a NEW local x, does NOT modify the global x
    print(x)  # 10

func()
print(x)  # 88 — global x is unchanged
```

### 2.2 Enclosing Scope (Nested Functions)

When functions are nested, the inner function can access variables from the enclosing (outer) function's scope:

```python
def fun1():
    x = 'fun1_local'

    def fun2():
        print(x)  # accesses enclosing scope → 'fun1_local'

    fun2()

fun1()
# fun2()  # NameError — fun2 only exists inside fun1
```

The enclosing scope is important for **closures**. A closure is a function that captures and remembers variables from its enclosing scope even after the enclosing function has finished executing.

### 2.3 The global Keyword

To modify a global variable from inside a function, you must use the `global` keyword:

```python
def modify_global():
    global x
    x = 99

x = 1
modify_global()
print(x)  # 99 — global x was modified
```

> **[Programming Languages]** The LEGB rule is Python's implementation of lexical (static) scoping. Most modern languages use similar scoping rules. The key point is that a function can read variables from outer scopes, but it needs an explicit `global` or `nonlocal` keyword to reassign those names.

Similarly, the `nonlocal` keyword allows modification of variables in the enclosing (but not global) scope:

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

> **Note:** Excessive use of `global` and `nonlocal` makes code harder to reason about and debug. Prefer passing values as arguments and returning results instead.

---

<br>

## 3. Arguments

Python offers flexible argument-passing mechanisms that go well beyond simple positional parameters.

### 3.1 Keyword Arguments and Defaults

Functions can have **default values** for parameters. Parameters with defaults are optional when calling the function:

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

Keyword arguments can be passed in any order, as long as they come after all positional arguments:

```python
def describe(name, age=25, city='Seoul'):
    return f'{name}, age {age}, from {city}'

print(describe('Alice', city='Sejong', age=22))
# Alice, age 22, from Sejong
```

### 3.2 Variable Positional Arguments (*args)

The `*args` syntax allows a function to accept any number of positional arguments. Inside the function, `args` is a tuple:

```python
def mysum(*args):
    return sum(args)

print(mysum(1, 2, 3))       # 6
print(mysum(10, 20, 30, 40)) # 100
```

This is useful when you do not know in advance how many arguments will be passed:

```python
def mean(*values):
    return sum(values) / len(values)

print(mean(1, 2, 3, 4, 5))  # 3.0
```

### 3.3 Variable Keyword Arguments (**kwargs)

The `**kwargs` syntax allows a function to accept any number of keyword arguments. Inside the function, `kwargs` is a dictionary:

```python
def show_kwargs(**kwargs):
    for key, val in kwargs.items():
        print(f'{key} = {val}')

show_kwargs(name='Alice', age=30, major='CSE')
# name = Alice
# age = 30
# major = CSE
```

You can combine all argument types in a single function definition. The order must be positional arguments, `*args`, keyword arguments with defaults, and `**kwargs`:

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

## 4. File I/O with NumPy

Numerical computing frequently involves reading data from files and saving results. NumPy provides efficient functions for this:

| Format | Save Function | Load Function | Use Case |
|:-------|:-------------|:-------------|:---------|
| `.txt` / `.csv` | `np.savetxt` | `np.loadtxt` | Human-readable, interoperable |
| `.npy` (single array) | `np.save` | `np.load` | Fast binary, single array |
| `.npz` (multiple arrays) | `np.savez` | `np.load` | Fast binary, multiple named arrays |

*Table 2. NumPy file I/O functions by format*

### 4.1 Text Files (CSV/TXT)

`np.loadtxt` reads columnar numeric data from text files. It handles headers, delimiters, and column selection:

```python
import numpy as np

# Load text data (e.g., KMA weather station data)
# File has 2 header rows; 3 columns: station_id, time, temperature
_, t, temp = np.loadtxt('./data/awsdata.txt', skiprows=2, unpack=True)
```

The parameters mean the following:

- `skiprows=2`: skips the first 2 rows (header lines).
- `unpack=True`: transposes the data so each column becomes a separate 1D array.
- The `_` variable is a Python convention for values we intend to discard (here, the station ID column).

To save data to a text file:

```python
data = np.column_stack([t, temp])
np.savetxt('output.csv', data, delimiter=',', header='time,temperature', fmt='%.4f')
```

### 4.2 Binary Files (NPY/NPZ)

Binary formats are much faster to read/write and preserve full floating-point precision (no formatting loss):

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

> **Note:** Binary files (`.npy`, `.npz`) are not human-readable, but they are significantly faster and more compact than text files. Use text formats when you need interoperability with other tools (Excel, R, etc.) and binary formats for intermediate results within a Python workflow.

---

<br>

## 5. Structured Programming

### 5.1 Conditional Statements

The `if`/`elif`/`else` structure provides branching logic:

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

### 5.2 For and While Loops

**For loops** iterate over a sequence:

```python
# Sum of integers 1 through 10
total = 0
for i in range(1, 11):
    total += i
print(total)  # 55
```

**`break`** exits the loop immediately, and **`continue`** skips the rest of the body and moves to the next iteration:

```python
# Find the first multiple of 7 greater than 0
for i in range(100):
    if i % 7 == 0 and i > 0:
        print(f'First multiple of 7: {i}')  # 7
        break
```

**While loops** repeat until a condition becomes `False`:

```python
count = 10
while count > 0:
    print(count, end=' ')
    count -= 1
# 10 9 8 7 6 5 4 3 2 1
```

Repeating until a condition is met is particularly useful in numerical methods for **iterative convergence**, where you repeat until the error drops below a tolerance. The example below uses a `for` loop with a maximum iteration count and `break` instead of `while`, so that it cannot loop forever if it fails to converge:

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

### 5.3 Vectorization vs. Loops

This is the single most important performance concept in scientific Python. **Vectorization** means replacing explicit Python loops with NumPy array operations, which execute in compiled C code:

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

> **[Computer Programming]** Vectorization is the single most important optimization technique in scientific Python. Instead of looping over elements in Python (slow interpreter), vectorized operations push the loop down to C (fast compiled code). Whenever possible, prefer NumPy vectorized operations over explicit loops.

Why is the loop so much slower? Python is an **interpreted language**, so each iteration of the loop involves type checking, function dispatch, and other overhead. NumPy's vectorized operations bypass this by executing a single optimized C routine on the entire array at once.

> **Note:** **Rule of thumb.** If you find yourself writing a `for` loop that processes array elements independently one at a time, there is almost certainly a vectorized NumPy equivalent that will be dramatically faster. A computation in which each step depends on the result of the previous step, such as Euler's method, cannot be vectorized directly.

---

<br>

## 6. Functions as First-Class Objects

In Python, functions are **first-class objects**. That is, they can be assigned to variables, passed as arguments, returned from other functions, and stored in data structures. This enables powerful functional programming patterns.

### 6.1 Lambda Functions

**Lambda functions** are anonymous, single-expression functions:

```python
f1 = lambda x, y: x ** 2 + y ** 2
print(f1(3, 4))  # 25
```

Lambdas are most useful as short, throwaway functions passed directly to other functions:

```python
data = [(1, 'b'), (3, 'a'), (2, 'c')]
sorted_data = sorted(data, key=lambda item: item[1])
print(sorted_data)  # [(3, 'a'), (1, 'b'), (2, 'c')]
```

### 6.2 Closures and Higher-Order Functions

A **closure** is a function that captures variables from its enclosing scope. A **higher-order function** is a function that takes another function as an argument or returns a function:

```python
# Return a function (closure)
def make_power(n):
    return lambda x: x ** n

square = make_power(2)
cube = make_power(3)
print(square(5))  # 25
print(cube(5))    # 125
```

Here, `make_power` returns a lambda that remembers the value of `n` from the enclosing scope, even after `make_power` has finished executing. This is the essence of a closure.

Passing a function as an argument looks like this:

```python
def f_at_midpoint(func, a, b):
    """Evaluate func at the midpoint of [a, b]."""
    mid = (a + b) / 2
    return func(mid)

import numpy as np
print(f_at_midpoint(np.cos, 0, 10))  # cos(5) = 0.2837...
print(f_at_midpoint(np.sin, 0, 10))  # sin(5) = -0.9589...
```

This pattern is fundamental to numerical methods. Many algorithms, such as root finders, integrators, and optimizers, take a function as input.

### 6.3 functools.partial

`functools.partial` creates a new function by fixing some arguments of an existing function:

```python
from functools import partial

def power(base, exp):
    return base ** exp

square = partial(power, exp=2)
cube = partial(power, exp=3)

print(square(5))  # 25
print(cube(5))    # 125
```

This is useful when you need to adapt a function's signature to match what another function expects. For example, an ODE solver might expect `f(t, y)`, but your function has extra parameters; `partial` can fix those extra parameters.

### 6.4 zip, filter, and Comprehensions

**`zip`** pairs elements from multiple iterables:

```python
names = ['Alice', 'Bob']
scores = [95, 87]

for name, score in zip(names, scores):
    print(f'{name}: {score}')
# Alice: 95
# Bob: 87
```

**`filter`** selects elements that satisfy a predicate:

```python
def even(x):
    return x % 2 == 0

print(list(filter(even, range(10))))  # [0, 2, 4, 6, 8]
```

**Comprehensions** are concise, Pythonic alternatives to `map` and `filter`:

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

> **Note:** Comprehensions are generally preferred over `map` and `filter` in modern Python because they are more readable. However, `map` can be marginally faster when it applies a built-in function, as in `map(str, data)`, while passing a lambda usually makes it slower than a comprehension because of the function call overhead.

---

<br>

## 7. Exception Handling

Real-world data and computations can produce errors. Python's `try`/`except` mechanism allows you to handle errors gracefully instead of crashing the program:

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

As a practical application, the following example cleans noisy data that may contain non-numeric values:

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

The `try`/`except` structure can also include `else` (runs if no exception occurred) and `finally` (always runs, useful for cleanup):

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

## 8. Decorators

A **decorator** is a function that takes another function, extends its behavior, and returns the modified function. Decorators use the `@` syntax:

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

The `@simple_decorator` syntax above `greet` is equivalent to writing:

```python
greet = simple_decorator(greet)
```

> **[Programming Languages]** Decorators are syntactic sugar for higher-order functions. `@decorator` above a function definition is equivalent to `func = decorator(func)`. They are widely used in web frameworks (Flask routes), testing, and logging.

As a practical example, a decorator that measures execution time looks like this:

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

## 9. OOP Basics

### 9.1 Classes with `__init__`, `__repr__`, `__str__`

Object-oriented programming organizes code around **objects** that bundle data (attributes) and behavior (methods):

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

The key special methods (dunder methods) are:

| Method | Purpose | Called by |
|:-------|:--------|:---------|
| `__init__` | Initialize the object | `MyStudent(...)` |
| `__repr__` | Unambiguous string representation (for developers) | `repr(obj)`, REPL |
| `__str__` | Human-readable string representation | `print(obj)`, `str(obj)` |

*Table 3. Key special methods and when they are called*

### 9.2 Dataclasses

The `@dataclass` decorator (Python 3.7+) eliminates boilerplate by automatically generating `__init__`, `__repr__`, and the equality method `__eq__`. Ordering methods are generated only when requested with `@dataclass(order=True)`:

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

The `@property` decorator turns a method into a read-only attribute, so you can access it without parentheses:

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

> **Note:** Dataclasses are the preferred way to create simple data-holding classes in modern Python. They reduce boilerplate and make the code cleaner and less error-prone compared to writing `__init__` manually.

---

<br>

## 10. Type Hints

**Type hints** (PEP 484) document the expected types of function parameters and return values. They do not enforce types at runtime but enable better IDE support, documentation, and static analysis:

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

Usage:

```python
train_idx, test_idx = train_test_split_idx(100, train_ratio=0.7, seed=42)
print(f'Train: {len(train_idx)} samples, Test: {len(test_idx)} samples')
# Train: 70 samples, Test: 30 samples
```

Common type hints:

| Hint | Meaning |
|:-----|:--------|
| `int`, `float`, `str` | Basic types |
| `List[int]` | List of integers |
| `Tuple[int, float]` | Tuple with specific element types |
| `Optional[int]` | `int` or `None` |
| `Dict[str, float]` | Dictionary with string keys and float values |
| `Callable[[int], float]` | Function taking `int`, returning `float` |

*Table 4. Common type hints and their meanings*

> **Note:** From Python 3.9+, you can use built-in types directly (`list[int]`, `tuple[int, float]`) instead of importing from `typing`. From Python 3.10+, `X | None` replaces `Optional[X]`.

---

<br>

## 11. Generators

### 11.1 Memory Efficiency

Generators produce values one at a time on demand through **lazy evaluation**, instead of storing all values in memory at once. This is critical when working with large datasets:

```python
import sys

# List comprehension: stores ALL values in memory
list_comp = [x ** 2 for x in range(10000)]

# Generator expression: stores only the formula, computes on demand
gen_expr = (x ** 2 for x in range(10000))

print(sys.getsizeof(list_comp))  # ~85,000 bytes
print(sys.getsizeof(gen_expr))   # ~200 bytes (!)
```

The difference is dramatic: the list stores references to 10,000 integers, while the generator stores only the state needed to produce the next value. `sys.getsizeof` measures only the list's array of references and does not include the memory taken by the integer objects themselves. The exact byte counts vary slightly between Python versions.

### 11.2 Generator Functions with yield

The `yield` keyword turns a regular function into a **generator function**. Each time `yield` is encountered, the function suspends execution and returns a value. On the next iteration, it resumes from where it left off:

```python
def countdown(n):
    while n > 0:
        yield n
        n -= 1

for val in countdown(5):
    print(val, end=' ')  # 5 4 3 2 1
```

Execution flow:
1. Calling `countdown(5)` returns a generator object without executing the body yet.
2. The first `next()` executes until `yield 5`, suspends, and returns `5`.
3. The second `next()` resumes after `yield`, decreases `n` to `4`, and returns `4` when it reaches `yield 4`.
4. The same process continues until the `while` condition is `False`.
5. When the function returns (implicitly), `StopIteration` is raised, ending the loop.

### 11.3 Batch Generator for Large Data

Generators are especially useful for processing large datasets in batches (common in machine learning and data analysis):

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

> **[Computer Programming]** Generators use lazy evaluation, so values are computed on demand rather than stored in memory. This is essential when working with datasets larger than available RAM. The `yield` keyword suspends function execution and resumes where it left off on the next call.

Another practical example is a generator that reads a large file line by line without loading it all into memory:

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

## Summary

| Concept | Key Summary |
|:------|:------------|
| Functions | Functions are defined with `def` and documented with docstrings, and they support multiple returns and `*args`/`**kwargs`. |
| Scope | Names are resolved by the LEGB rule, `global` and `nonlocal` allow assignment to outer names, and closures remember outer variables. |
| File I/O | Text files are handled with `np.loadtxt`, and binary files with `np.save`/`np.load` and `np.savez`. |
| Control Flow | `if`/`elif`/`else`, `for`, `while`, and `break`/`continue` control branching and repetition. |
| Vectorization | NumPy vectorized operations are typically 10-100x faster than Python loops. |
| Functional | `lambda`, `partial`, higher-order functions, and comprehensions treat functions as values. |
| Exceptions | `try`/`except`/`else`/`finally` handles errors. |
| Decorators | `@decorator` is syntactic sugar for `func = decorator(func)`. |
| OOP | `__init__` initializes an object, and `@dataclass` and `@property` reduce boilerplate. |
| Type Hints | `int`, `List[int]`, `Optional[...]`, and `Tuple[...]` document expected types, which are not enforced at runtime. |
| Generators | `yield` provides memory-efficient iteration through lazy evaluation. |

---

<br>

## Self-Check Questions

1. **LEGB rule:** A global `x = 'global'` exists, and the function `outer` contains `x = 'enclosing'` and a nested function `inner`. If `inner` only executes `print(x)`, what is printed when `outer()` is called?

   <details>
   <summary>Answer</summary>

   `'enclosing'` is printed. `inner` has no local `x`, so Python looks in the next scope, the enclosing scope of `outer`, finds `x` there, and never reaches the global scope.

   </details>

2. **Local and global variables:** A global `count = 0` exists, and the function `f` executes `count = 5` without a `global` declaration. What does `print(count)` show after calling `f()`?

   <details>
   <summary>Answer</summary>

   It prints `0`. The `count` assigned inside the function becomes a new local variable and does not change the global `count`. To change the global variable, `global count` must be declared.

   </details>

3. **Multiple return values and Bessel's correction:** What tuple does the `stats` function of Section 1.3 return for `[2, 4, 6]`?

   <details>
   <summary>Answer</summary>

   It returns `(3, 4.0, 2.0)`. The mean is 12/3 = 4.0, and the sum of squared deviations is 4 + 0 + 4 = 8. Dividing by $`n - 1 = 2`$ gives 4, whose square root, the sample standard deviation, is 2.0.

   </details>

4. **Combining arguments:** If the `combined` function of Section 3.3 is called as `combined(1, 2, 3, option=False, z=0)`, what are `args`, `option`, and `kwargs`?

   <details>
   <summary>Answer</summary>

   After `a=1` and `b=2` are filled, the remaining positional argument 3 goes into `args`, so `args=(3,)`. `option=False`, and the keyword argument `z=0`, which matches no parameter name, goes into `kwargs={'z': 0}`.

   </details>

5. **Batch generator:** List, in order, the sizes of the batches produced by `batch_generator(np.arange(10), 4)` from Section 11.3.

   <details>
   <summary>Answer</summary>

   The sizes are 4, 4, and 2. The start indices are 0, 4, and 8, so `data[0:4]`, `data[4:8]`, and `data[8:12]` are produced, and the last slice is cut off at the end of the array and holds 2 elements.

   </details>
