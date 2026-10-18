# [Spring 2026] Numerical Analysis

![Last Commit](https://img.shields.io/github/last-commit/Choroning/26Spring_Numerical-Analysis)
![Languages](https://img.shields.io/github/languages/top/Choroning/26Spring_Numerical-Analysis)

This repository archives bilingual lecture notes and coursework for the Numerical Analysis course.

*Author: Cheolwon Park (Korea University Sejong, CSE) - Year 3 (Junior) as of 2026*
<br><br>

## 📑 Table of Contents

- [About This Repository](#about-this-repository)
- [Course Information](#course-information)
- [Prerequisites](#prerequisites)
- [Repository Structure](#repository-structure)
- [License](#license)

---


<br><a name="about-this-repository"></a>
## 📝 About This Repository

This repository contains the bilingual lecture notes of a Numerical Analysis course that covers the fundamental numerical methods for solving scientific and engineering problems, following the textbook by Chapra and Clough.

- `Lecture-Notes/` holds a lecture note (`Concepts_Lecture`) for each textbook chapter, and Chapters 1-3 and 8-13 also have a lab note (`Concepts_Lab`) that summarizes the lab content.
- The lecture slides follow the textbook chapter numbers, so the units are also numbered by chapter (`Ch##`), and they cover the full curriculum across the six parts of the textbook. The part grouping is in the `Lecture-Notes` index.
- Every lecture note and lab note comes as an English file (`.md`) and a Korean file (`.ko.md`) with the same structure.
- The lecture notes include Python implementations of the numerical algorithms from the textbook.
- `Coursework/` holds the lab notebooks (Jupyter) completed and submitted as the 12 weekly homework assignments.

> **🤖 AI Agent Usage**
> This course encourages the use of AI agents.
> [Claude Code](https://claude.ai/download) and [Gemini CLI](https://github.com/google-gemini/gemini-cli) were used as coding assistants throughout the course, and Claude Code was also used to organize and translate the lecture notes and to restructure the repository.

<br><a name="course-information"></a>
## 📚 Course Information

- **Semester:** Spring 2026 (March - June)
- **Affiliation:** Korea University Sejong

| Course&nbsp;Code| Course            | Type          | Instructor      | Department                              |
|:----------:|:------------------|:-------------:|:---------------:|:----------------------------------------|
|`DCSS305-00`|NUMERICAL ANALYSIS (English)|Major Elective|Prof. Shinhoo&nbsp;Kang|Department of Computer Science and Software Engineering|

### Course Overview

An English-taught major elective that surveys the fundamental numerical methods for solving scientific and engineering problems and implements each method in Python following the textbook by Chapra and Clough. Numerical methods are widely used in science and engineering applications, including data science, life science, social science, atmospheric science, finance, aerospace engineering, mechanical engineering, and biomedical engineering. The course covers error analysis, solving linear systems and nonlinear equations, eigenvalues, optimization, regression, Fourier analysis, interpolation, numerical differentiation and integration, and ordinary differential equations, and aims to build the ability to formulate and solve complex scientific and engineering problems by integrating fundamental numerical methods with computational application skills. The instructor, schedule and class format, assessment, course policies, weekly plan and progress, and references are in the [course information document](Course-Information.md).

<br><a name="prerequisites"></a>
## ✅ Prerequisites

- An understanding of linear algebra and discrete mathematics is required. The syllabus also lists Python as a recommended prerequisite.
- A Python interpreter must be installed.
- Familiarity with scientific computing libraries is required.

- **💻 Development Environment**

| Tool | Company |  OS  | Notes |
|:-----|:-------:|:----:|:------|
|Python 3|Python Software Foundation|macOS|    |
|NumPy|NumFOCUS|macOS|Numerical computing|
|SciPy|NumFOCUS|macOS|Scientific computing|
|Matplotlib|NumFOCUS|macOS|Visualization|
|JupyterLab|Project Jupyter|macOS|Interactive notebooks|

<br><a name="repository-structure"></a>
## 🗂 Repository Structure

```plaintext
26Spring_Numerical-Analysis
├── Lecture-Notes               Lecture notes and lab notes by textbook chapter
├── Coursework                  Weekly homework notebooks
├── Course-Information.ko.md    Course operations
├── Course-Information.md
├── LICENSE                     Documentation license (CC BY-NC-SA 4.0)
├── LICENSE-CODE                Code license (PolyForm Noncommercial 1.0.0)
├── README.ko.md
└── README.md
```

The units and the part grouping of the textbook are listed in the [`Lecture-Notes` index](Lecture-Notes/README.md). The homework is listed in the [`Coursework` index](Coursework/README.md).

<br><a name="license"></a>
## 🤝 License

The documents in this repository (lecture notes, explanations, and summary sheets) are licensed under [CC BY-NC-SA 4.0](LICENSE), and the code under [PolyForm Noncommercial 1.0.0](LICENSE-CODE). Figures and example code taken from lecture slides and textbooks, as well as assignment handouts and template code distributed by the course, belong to their original authors and are not covered by these licenses.

---
