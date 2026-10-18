# [2026학년도 봄학기] 수치해석

![Last Commit](https://img.shields.io/github/last-commit/Choroning/26Spring_Numerical-Analysis)
![Languages](https://img.shields.io/github/languages/top/Choroning/26Spring_Numerical-Analysis)

이 레포지토리는 수치해석 과목의 이중 언어 강의 노트와 수업 결과물을 정리한 아카이브입니다.

*작성자: 박철원 (고려대학교(세종), 컴퓨터소프트웨어학과) - 2026년 기준 3학년*
<br><br>

## 📑 목차

- [레포지토리 소개](#about-this-repository)
- [강의 정보](#course-information)
- [사전 요구사항](#prerequisites)
- [레포지토리 구조](#repository-structure)
- [라이선스](#license)

---


<br><a name="about-this-repository"></a>
## 📝 레포지토리 소개

이 레포지토리에는 Chapra와 Clough의 교재를 따라 과학과 공학 문제를 푸는 기본 수치 방법을 다루는 수치해석 과목의 이중 언어 강의 노트가 담겨 있습니다.

- `Lecture-Notes/`에는 교재의 장마다 강의를 정리한 강의 노트(`Concepts_Lecture`)가 있고, 1\~3장과 8\~13장에는 실습 내용을 정리한 실습 노트(`Concepts_Lab`)가 함께 있습니다.
- 강의자료의 장 번호가 교재의 장 번호와 같아 단원도 장 번호(`Ch##`)로 나눴으며, 교재의 6개 부(Part)에 걸친 전체 커리큘럼을 다룹니다. 부 구분은 `Lecture-Notes` 색인에 있습니다.
- 모든 강의 노트와 실습 노트는 영어 파일(`.md`)과 한국어 파일(`.ko.md`)이 같은 구성으로 짝을 이룹니다.
- 강의 노트에는 교재의 수치 알고리즘을 Python으로 구현한 코드가 함께 실려 있습니다.
- `Coursework/`에는 주간 과제 12개로 완성해 제출한 실습 노트북(Jupyter)이 있습니다.

> **🤖 AI 에이전트 활용**
> 이 과목은 AI 에이전트 사용을 권장합니다.
> 수업 전반에 걸쳐 [Claude Code](https://claude.ai/download)와 [Gemini CLI](https://github.com/google-gemini/gemini-cli)를 코딩 어시스턴트로 활용했고, 강의 노트의 정리와 번역, 레포지토리 구조 정리에도 Claude Code를 활용했습니다.

<br><a name="course-information"></a>
## 📚 강의 정보

- **학기:** 2026학년도 봄학기 (3월 - 6월)
- **소속:** 고려대학교(세종)

|학수번호      |강의명    |이수구분|교수자|개설학과|
|:----------:|:-------|:----:|:------:|:----------------|
|`DCSS305-00`|수치해석(영강)|전공선택|강신후 교수|컴퓨터소프트웨어학과|

### 과목 개요

과학과 공학 문제를 푸는 기본 수치 방법을 개관하는 영어 강의 전공선택 과목으로, Chapra와 Clough의 교재를 따라 각 방법을 Python으로 구현한다. 수치 방법은 데이터 과학, 생명 과학, 사회 과학, 대기 과학, 금융, 항공우주공학, 기계공학, 의공학 등 과학과 공학의 여러 응용에 널리 쓰인다. 오차 해석, 선형 시스템과 비선형 방정식의 풀이, 고유값, 최적화, 회귀, 푸리에 해석, 보간, 수치 미분과 수치 적분, 상미분방정식을 다루며, 기본 수치 방법과 계산 응용 능력을 결합해 과학과 공학의 복잡한 문제를 정식화하고 푸는 능력을 기르는 것을 목표로 한다. 교수자, 수업 일정과 형식, 평가, 수업 운영 규정, 주차별 계획 및 진도, 참고 자료는 [강의 정보 문서](Course-Information.ko.md)에 정리했다.

<br><a name="prerequisites"></a>
## ✅ 사전 요구사항

- 선형대수와 이산수학을 이해하고 있어야 합니다. 강의계획서는 Python도 권장 선수 지식으로 꼽습니다.
- Python 인터프리터가 설치되어 있어야 합니다.
- 과학 계산 라이브러리에 익숙해야 합니다.

- **💻 개발 환경**

| 도구 | 회사 |  운영체제  | 비고 |
|:-----|:-------:|:----:|:------|
|Python 3|Python Software Foundation|macOS|    |
|NumPy|NumFOCUS|macOS|수치 연산|
|SciPy|NumFOCUS|macOS|과학 연산|
|Matplotlib|NumFOCUS|macOS|시각화|
|JupyterLab|Project Jupyter|macOS|대화형 노트북|

<br><a name="repository-structure"></a>
## 🗂 레포지토리 구조

```plaintext
26Spring_Numerical-Analysis
├── Lecture-Notes               교재 장별 강의 노트와 실습 노트
├── Coursework                  주간 과제 노트북
├── Course-Information.ko.md    수업 운영 정보
├── Course-Information.md
├── LICENSE                     문서 라이선스 (CC BY-NC-SA 4.0)
├── LICENSE-CODE                코드 라이선스 (PolyForm Noncommercial 1.0.0)
├── README.ko.md
└── README.md
```

단원 목록과 교재의 부(Part) 구분은 [`Lecture-Notes` 색인](Lecture-Notes/README.ko.md)에 있습니다. 과제 목록은 [`Coursework` 색인](Coursework/README.ko.md)에 있습니다.

<br><a name="license"></a>
## 🤝 라이선스

이 레포지토리의 문서(강의 노트, 해설, 요약지)는 [CC BY-NC-SA 4.0](LICENSE), 코드는 [PolyForm Noncommercial 1.0.0](LICENSE-CODE)에 따라 배포됩니다. 강의자료와 교재에서 가져온 그림과 예제 코드, 강의가 배포한 과제 원본과 템플릿 코드는 원저작자에게 권리가 있으며 이 라이선스의 대상이 아닙니다.

---
