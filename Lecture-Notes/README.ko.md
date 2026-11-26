# 강의 노트

[레포지토리 안내](../README.ko.md) | [English](README.md) | [한국어](README.ko.md)

> **최종 수정일:** 2026-11-27

## 1부 모델링, 컴퓨터, 오차 분석

| 단원 | 주제 | 문서 |
|:-----|:-----|:-----|
| [1장 수학적 모델링과 공학 문제 해결](Ch01_Mathematical-Modeling/) | 수학적 모델, Newton의 제2법칙과 번지 점프 예제, 보존 법칙과 정상 상태, 이차 항력 모델의 해석해, Euler 방법과 그 국소 오차와 전역 오차, 수치적으로 푸는 공학 문제의 유형, 과도 계산과 정상 상태 계산, 낙하 거리, 위로 출발하는 경우의 해석해와 Euler 해, 수치 오차와 모델 오차 | [강의 노트](Ch01_Mathematical-Modeling/Concepts_Lecture.ko.md) |
| [2장 Python 기초](Ch02_Python-Fundamentals/) | 개발 환경과 노트북, Python 기본 문법, 대입과 별칭, 내장 타입과 문자열, 수치 자료형의 성질, NumPy 배열 연산과 브로드캐스팅, 벡터화와 범용 함수, NumPy와 SciPy 함수, Matplotlib 그래프와 로그 축, NumPy의 뷰와 복사본, 연산자 우선순위와 음수의 거듭제곱, np.matrix와 ndarray, math 함수와 atan2, 히스토그램, 탐색적 자료 분석 | [강의 노트](Ch02_Python-Fundamentals/Concepts_Lecture.ko.md), [실습 노트](Ch02_Python-Fundamentals/Concepts_Lab.ko.md) |
| [3장 Python 프로그래밍](Ch03_Python-Programming/) | 구조적 프로그래밍과 의사코드, 제어 구조와 반복 계산의 종료 조건, 위치 인수, 키워드 인수, 가변 인수, 함수를 인수로 전달하기와 map, 변수 스코프와 LEGB 규칙, 모듈과 독스트링, 번지 점프 예제를 푸는 일반화된 Euler 함수, 파일 입출력, 예외 처리, 수치 프로그램의 검증, 데코레이터, 객체 지향 기초, 타입 힌트, 제너레이터, 논리 연산자의 우선순위와 조건의 순서, 사전 검사와 사후 검사 반복과 else 절, 스텝 수와 시각 계산의 함정, 벡터화할 수 있는 계산과 없는 계산 | [강의 노트](Ch03_Python-Programming/Concepts_Lecture.ko.md), [실습 노트](Ch03_Python-Programming/Concepts_Lab.ko.md) |
| [4장 반올림 오차와 절단 오차](Ch04_Roundoff-and-Truncation-Errors/) | 유효 숫자, 정확도와 정밀도, 편향과 불확실성, 오차의 정의와 Scarborough 기준, 오차의 원인, 정수 표현과 2의 보수, 부동 소수점 표현과 편향 지수, 비정규화 수와 단위 반올림, 부동 소수점 산술의 문제점, Taylor 급수와 나머지 항, 근사의 차수, 유한 차분 근사, 총 수치 오차와 최적 스텝 크기, 안정성과 조건수, 수치 오차의 관리 | [강의 노트](Ch04_Roundoff-and-Truncation-Errors/Concepts_Lecture.ko.md), [실습 노트](Ch04_Roundoff-and-Truncation-Errors/Concepts_Lab.ko.md) |

## 2부 근과 최적화

| 단원 | 주제 | 문서 |
|:-----|:-----|:-----|
| [5장 근: 구간법](Ch05_Roots-Bracketing-Methods/) | 근 문제와 그래프 방법, 시행착오, 구간법과 개방법의 비교, 증분 탐색과 그 구현, 이분법과 그 수렴, 수렴 차수, 근사 상대 오차와 중지 기준, 가위치법과 수정 가위치법, 가위치법의 수렴 속도, SciPy의 bisect와 brentq | [강의 노트](Ch05_Roots-Bracketing-Methods/Concepts_Lecture.ko.md) |
| [6장 근: 개방법](Ch06_Roots-Open-Methods/) | 고정점 반복의 그래프 해석과 오차 추정, Wegstein 방법, Newton-Raphson법과 이차 수렴, Newton-Raphson법의 실패, 중근과 수정 Newton-Raphson법, 연립 비선형 방정식으로의 확장, 할선법과 그 수렴 차수, 역 이차 보간, Brent 방법, SciPy의 근 찾기 함수, 다항식의 근과 동반 행렬 | [강의 노트](Ch06_Roots-Open-Methods/Concepts_Lecture.ko.md) |
| [7장 최적화](Ch07_Optimization/) | 국소 최적성과 전역 최적성, 최적해의 존재, 1차원 최적화, 황금 분할 탐색의 반복 횟수와 정밀도 한계, 포물선 보간과 Brent 방법, 등고선과 기울기, 볼록 최적화 문제의 1차와 2차 특성화, 최급 하강법과 그 수렴 속도, 최적화의 Newton 방법과 정류점의 판정, SciPy 최적화와 제약 최적화 | [강의 노트](Ch07_Optimization/Concepts_Lecture.ko.md) |

## 3부 선형계

| 단원 | 주제 | 문서 |
|:-----|:-----|:-----|
| [8장 연립 일차 대수방정식의 풀이](Ch08_Linear-Algebraic-Equations/) | 행렬 표기와 특수 행렬, NumPy의 행렬 표기, 행렬 연산 규칙, 연립방정식의 행렬 형태, 정방 시스템의 해의 존재와 유일성, 소규모 시스템의 그래프 방법과 열 관점, 행렬식과 Cramer 법칙, 단순 Gauss 소거법과 계산 비용, 단순 Gauss 소거법의 함정, 악조건 시스템과 잔차, 질량-스프링 시스템과 회로의 메시 해석 | [강의 노트](Ch08_Linear-Algebraic-Equations/Concepts_Lecture.ko.md), [실습 노트](Ch08_Linear-Algebraic-Equations/Concepts_Lab.ko.md) |
| [9장 Gauss 소거법: 개선된 방법](Ch09_Gauss-Elimination/) | 부분 피봇팅의 단계별 추적, 성장 인자, 스케일링된 부분 피봇팅, Gauss-Jordan 방법과 그 연산 횟수, 삼중대각 시스템과 Thomas 알고리즘, Thomas 알고리즘의 정규화 형태와 안정성 조건, 열전도 막대의 삼중대각 시스템, Gauss 소거법을 이용한 행렬식 계산, 연산 횟수의 유도, 조건수 | [강의 노트](Ch09_Gauss-Elimination/Concepts_Lecture.ko.md), [실습 노트](Ch09_Gauss-Elimination/Concepts_Lab.ko.md) |
| [10장 LU 분해](Ch10_LU-Factorization/) | Gauss 소거법과 LU 분해의 관계, 기본 행렬, LU 분해의 존재 조건, 전진 대입과 후진 대입, 피봇팅을 적용한 LU 분해와 승수 교환, Cholesky 분해와 공식의 유도, Doolittle 알고리즘과 LU 분해의 연산 횟수, scipy.linalg.lu의 규약, lu_factor와 lu_solve, 다중 우변 문제와 열저항 네트워크 | [강의 노트](Ch10_LU-Factorization/Concepts_Lecture.ko.md), [실습 노트](Ch10_LU-Factorization/Concepts_Lab.ko.md) |
| [11장 역행렬과 조건수](Ch11_Matrix-Inverse-and-Condition/) | 역행렬과 LU 분해를 이용한 계산, 역행렬 계산의 비용, 질량-스프링 시스템의 물리적 해석, 자극과 응답 해석, 악조건의 판별, 벡터 노름과 행렬 노름, 노름의 동등성과 행렬 노름의 부등식, 행렬 조건수와 오차 한계 정리의 유도, 특이값과 조건수, 잔차와 오차, Hilbert 행렬, 중첩과 비례성, 계수 변화의 1차 추정, 스펙트럼 반지름, 조건수의 성질과 행 스케일링, 한계가 등호가 되는 섭동, Vandermonde 행렬과 정밀도 규칙 | [강의 노트](Ch11_Matrix-Inverse-and-Condition/Concepts_Lecture.ko.md), [실습 노트](Ch11_Matrix-Inverse-and-Condition/Concepts_Lab.ko.md) |
| [12장 선형계의 반복법](Ch12_Iterative-Methods/) | Gauss-Seidel 방법, 대각 우위와 수렴 정리, Jacobi 반복법, 완화법과 최적 완화 계수, 반복법의 수렴 속도, 전처리, 비선형 연립방정식의 연속 대입과 Newton-Raphson법, 비선형 반복의 수렴 판정, scipy.optimize.root | [강의 노트](Ch12_Iterative-Methods/Concepts_Lecture.ko.md) |
| [13장 고유값 방법: 대칭 행렬](Ch13_Eigenvalues-Symmetric-Matrices/) | 고유값과 고유벡터의 계산, 고유값과 상미분방정식의 안정성, 결합 2차 ODE 시스템, 질량-스프링 시스템, 거듭제곱법과 그 수렴 원리, Rayleigh 몫, 역거듭제곱법과 편이 역거듭제곱법, 수축법, 대칭 행렬의 실수 고유값과 직교 고유벡터, 스펙트럼 분해, 진동 모드와 고유벡터의 정규화, 고차 특성 다항식과 근의 민감도, 복소 고유값과 감쇠 진동, 질량 행렬과 강성 행렬, 대칭화와 M-직교성, 3층 건물의 고유 진동 | [강의 노트](Ch13_Eigenvalues-Symmetric-Matrices/Concepts_Lecture.ko.md), [실습 노트](Ch13_Eigenvalues-Symmetric-Matrices/Concepts_Lab.ko.md) |

## 4부 곡선 적합

| 단원 | 주제 | 문서 |
|:-----|:-----|:-----|
| [14장 선형 회귀](Ch14_Linear-Regression/) | 통계 복습과 정규 분포, 난수와 Monte Carlo 시뮬레이션, 직선 최소제곱 회귀와 단계별 예제, 정규 방정식의 행렬 형태, 최대 우도 해석과 회귀의 가정, 적합도의 정량화, 잔차 분석, 계수의 신뢰 구간, 거듭제곱과 포화 성장률 모델의 선형화, np.polyfit, MAD 기반 표준 편차, 분산 계산 공식의 뺄셈 소거, Python 기술 통계와 ddof, 정규 확률 규칙과 평균의 신뢰 구간, 히스토그램, 이상값과 지렛대, 역변환의 치우침, 효소 반응 속도 모델의 비교, 균일 분포와 오차 전파 | [강의 노트](Ch14_Linear-Regression/Concepts_Lecture.ko.md) |
| [15장 비선형 회귀](Ch15_Nonlinear-Regression/) | 표준 오차와 결정 계수, 비선형 관계의 선형화, 다항 회귀, 다중 선형 회귀와 예제, 임의 기저 함수의 일반 선형 최소제곱, 정규 방정식의 조건과 QR 분해, 수정 결정 계수와 예측 결정 계수, Gauss-Newton 비선형 회귀와 예제, 선형화 적합과 비선형 적합의 비교, curve_fit과 least_squares, minimize를 이용한 비선형 회귀, 공선성과 표준화, 다변수 거듭제곱 모델, 단계적 회귀, 매개변수의 상관, curve_fit의 경계, 선형화 근사와 Monte Carlo 불확실성 | [강의 노트](Ch15_Nonlinear-Regression/Concepts_Lecture.ko.md) |
| [16장 Fourier 해석](Ch16_Fourier-Analysis/) | 우함수와 기함수, Fourier 급수와 Gibbs 현상, 복소 Fourier 계수, Fourier 적분, Fourier 변환과 직사각형 펄스, 정현파의 최소제곱 적합과 예제, 주파수 영역과 시간 영역, 이산 Fourier 변환, Nyquist 주파수와 앨리어싱, 고속 Fourier 변환의 분할과 연산 횟수, 파워 스펙트럼, 스펙트럼 누설과 Parseval 정리, 잡음 속 주기 검출, 삼각파의 Fourier 급수, Fourier 적분의 계수 함수, 선 스펙트럼과 위상, 급수와 변환과 DFT의 관계, DFT 칸의 해석, 앨리어싱의 위상 반전, 추세 제거 | [강의 노트](Ch16_Fourier-Analysis/Concepts_Lecture.ko.md) |
| [17장 다항식 보간](Ch17_Interpolation/) | 다항식 보간과 그 유일성, Horner 방법, Newton 보간 다항식, Lagrange 보간 다항식, 보간 오차와 오차 추정, 역보간, 외삽과 진동, Chebyshev 노드, SciPy의 다항식 보간 | [강의 노트](Ch17_Interpolation/Concepts_Lecture.ko.md) |
| [18장 스플라인과 구간별 보간](Ch18_Splines/) | k차 스플라인과 연속성 조건, 선형 스플라인, 이차 스플라인, 삼차 스플라인, 고정 끝 조건과 not-a-knot 끝 조건, Thomas 알고리즘을 이용한 구현, CubicSpline과 수렴 차수, PCHIP, 평활 스플라인과 LOESS, 쌍선형 보간 | [강의 노트](Ch18_Splines/Concepts_Lecture.ko.md) |

## 5부 적분과 미분

| 단원 | 주제 | 문서 |
|:-----|:-----|:-----|
| [19장 수치 적분 공식](Ch19_Numerical-Integration-Formulas/) | 정확도 차수와 미정 계수법, 사다리꼴 법칙, Simpson의 1/3 법칙과 3/8 법칙, 오차 공식의 유도, 복합 적분 공식과 수렴 차수의 수치 확인, 적분 법칙의 결합, 비등간격 데이터, 누적 적분, 이중 적분과 dblquad, Boole 법칙과 열린 공식, 연속 함수의 평균과 총량, 닫힌 공식과 열린 공식의 가중치와 오차 상수, 이상 적분과 복합 중점 법칙, 사례 연구: 변하는 힘이 한 일 | [강의 노트](Ch19_Numerical-Integration-Formulas/Concepts_Lecture.ko.md) |
| [20장 함수의 수치 적분](Ch20_Numerical-Integration-of-Functions/) | Romberg 적분과 Newton-Cotes 공식의 관계, Gauss 구적법, Legendre 다항식과 Gauss 절점, 고차 Gauss-Legendre 절점과 가중치, Gauss 구적법의 오차, 적응적 구적법, SciPy의 quad, 함숫값을 재사용하는 Romberg 알고리즘, quad의 인수와 출력과 실패 진단, 사례 연구: 실효 전류 | [강의 노트](Ch20_Numerical-Integration-of-Functions/Concepts_Lecture.ko.md) |
| [21장 수치 미분](Ch21_Numerical-Differentiation/) | 1계 도함수와 2계 도함수의 유한 차분 근사와 오차, 최적 스텝 크기, 고차 전진 차분과 5점 중심 차분, Richardson 외삽법, 비등간격 데이터의 도함수, 잡음이 있는 데이터의 미분, 수치 편도함수, 미정 계수법과 비중심 스텐실, NumPy의 diff와 gradient, 미분 가능 프로그래밍 | [강의 노트](Ch21_Numerical-Differentiation/Concepts_Lecture.ko.md) |

## 6부 상미분방정식

| 단원 | 주제 | 문서 |
|:-----|:-----|:-----|
| [22장 상미분방정식의 초기값 문제](Ch22_Initial-Value-Problems/) | Euler 방법, 국소 절단 오차와 전역 절단 오차의 추적, 2차 Runge-Kutta 방법과 Ralston 방법, 고전적 4차 Runge-Kutta 방법, 방법별 수렴 차수와 비용, 연립 ODE와 그 손 계산, 수치적 안정성과 안정 영역, SciPy의 solve_ivp | [강의 노트](Ch22_Initial-Value-Problems/Concepts_Lecture.ko.md) |
| [23장 강성 ODE와 적응적 방법](Ch23_Stiff-ODEs/) | 적응적 Runge-Kutta 방법, 내장 RK 방법 RKF45와 오차 추정, 스텝 크기 제어, 다단계 방법, 강성 시스템과 고유값의 강성비, 연립 후진 Euler 방법과 Newton 방법, 후진 차분 공식(BDF), A-안정성과 음적 방법의 안정 영역, solve_ivp의 강성 풀이기, 사건 감지, 스텝 반감법, 좁은 펄스와 max_step, 교정자 반복, 사례 연구: 줄이 있는 번지 점프와 간헐 분수 | [강의 노트](Ch23_Stiff-ODEs/Concepts_Lecture.ko.md) |
| [24장 경계값 문제](Ch24_Boundary-Value-Problems/) | 경계값 문제와 해석해, 선형 사격법과 비선형 사격법, 도함수 경계 조건의 사격법, 유한 차분법, 도함수 경계 조건과 가상 노드, 비선형 유한 차분법, 사격법과 유한 차분법의 비교, SciPy의 solve_bvp, 열 수지로부터의 지배 방정식, 대류(Robin) 경계 조건의 사격법과 유한 차분법, 비선형 차분 방정식의 연속 대입 | [강의 노트](Ch24_Boundary-Value-Problems/Concepts_Lecture.ko.md) |
