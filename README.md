# Label Consistent Fisher Vectors (LCFV)

## Overview

LCFV is a method to add supervised information to Fisher vectors. With this package, you can compute a transformation matrix to be applied to Fisher vectors. LCFV takes the original Fisher vectors and class labels as input.
In this package, we do not provide the code for computing Fisher vectors. You need to compute the Fisher vectors yourself before using our package. For example, you can use the INRIA's Fisher vector implementation: http://lear.inrialpes.fr/src/inria_fisher/

The method is described in our ICPR 2014 paper:

```
Quan Wang, Xin Shen, Meng Wang, Kim L. Boyer, "Label Consistent Fisher Vectors for Supervised Feature Aggregation", 22nd International Conference on Pattern Recognition (ICPR), 2014.
```

This library is also available at MathWorks:
* https://www.mathworks.com/matlabcentral/fileexchange/47730-label-consistent-fisher-vectors-lcfv