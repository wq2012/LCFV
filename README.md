# Label Consistent Fisher Vectors (LCFV)

[![View Label Consistent Fisher Vectors (LCFV) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/47730-label-consistent-fisher-vectors-lcfv)
[![Octave application](https://github.com/wq2012/LCFV/actions/workflows/octave.yml/badge.svg)](https://github.com/wq2012/LCFV/actions/workflows/octave.yml)

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
  - [Prepare Data](#prepare-data)
  - [LCFV1](#lcfv1)
  - [LCFV2](#lcfv2)
- [Copyright and Citation](#copyright-and-citation)

## Overview

LCFV is a method for adding supervised information to Fisher vectors. This package allows you to compute a transformation matrix to be applied to Fisher vectors, taking the original Fisher vectors and class labels as input.

**Note**: This package does **not** provide code for computing Fisher vectors. You must compute them yourself before using this package (e.g., using [INRIA's Fisher vector implementation](http://lear.inrialpes.fr/src/inria_fisher/)).

![LCFV Logo](resources/LCFV_logo.png)

## Installation

You can clone this repository or download it from MATLAB File Exchange.

```bash
git clone https://github.com/wangquan/LCFV.git
```

## Usage

Please check `code/run_demo.m` for a complete example.

### Prepare Data

You need:
- `fv`: Fisher vectors matrix (dimension `N x D` or transposed).
- `labels`: Class labels.

### LCFV1

```matlab
% G: Fisher vectors (D x N)
% C: Label comparison matrix
% alpha: Parameter (tune this!)

[M1, W1] = solve_LCFV1(G, C, alpha);
LCFV1 = M1 * G;
```

### LCFV2

```matlab
M2 = solve_LCFV2(G, C, alpha);
LCFV2 = M2 * G;
```

## Copyright and Citation

```
Copyright (C) 2012 Quan Wang <wangq10@rpi.edu>,
Signal Analysis and Machine Perception Laboratory,
Department of Electrical, Computer, and Systems Engineering,
Rensselaer Polytechnic Institute, Troy, NY 12180, USA
```

This software was developed as part of the following research. If you use this software in your research, please cite:

**Plain Text:**

> Quan Wang, Xin Shen, Meng Wang, and Kim L. Boyer.
"Label consistent fisher vectors for supervised feature aggregation."
In 2014 22nd International Conference on Pattern Recognition, pp. 3588-3593. IEEE, 2014.

> Quan Wang.
Exploiting Geometric and Spatial Constraints for Vision and Lighting Applications.
Ph.D. dissertation, Rensselaer Polytechnic Institute, 2014.

**BibTeX:**

```bibtex
@inproceedings{wang2014label,
  title={Label consistent Fisher vectors for supervised feature aggregation},
  author={Wang, Quan and Shen, Xin and Wang, Meng and Boyer, Kim L},
  booktitle={Pattern Recognition (ICPR), 2014 22nd International Conference on},
  pages={2507--2512},
  year={2014},
  organization={IEEE}
}

@phdthesis{wang2014exploiting,
  title={Exploiting Geometric and Spatial Constraints for Vision and Lighting Applications},
  author={Quan Wang},
  year={2014},
  school={Rensselaer Polytechnic Institute},
}
```

This library is also available at MathWorks:
https://www.mathworks.com/matlabcentral/fileexchange/47730-label-consistent-fisher-vectors-lcfv

