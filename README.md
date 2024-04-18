# Construction-of-the-best-approximation-of-a-function-through-basis-functions-of-finite-elements
>>>A program written in C++ for constructing the best approximation of a function through basis functions of finite elements. There is work with the matrix, which is stored in tape form. And four methods are used for solving systems of linear algebraic equations (including iterative ones)
___
Since in this case a strip matrix is ​​obtained, it is much more efficient to store only its non-zero elements. But since it is also symmetrical, we will store only the upper part of the tape. This makes the algorithm significantly more memory efficient.
![image](https://github.com/chay-s-penkoy/Construction-of-the-best-approximation-of-a-function-through-basis-functions-of-finite-elements/assets/129063012/4a3a3bed-3293-4a2f-9613-4658f42f4cc9)

___
***1.*** A function is introduced to calculate *y = f(x)*.

***2.*** The coordinates of the ends of the segment [a,b] are specified.

***3.*** The number K of finite elements (intervals) is specified.

***4.*** The number of nodes N of the finite element is specified.

***5.*** The total number of grid nodes M on the segment [a,b] is calculated and
a uniform grid with step h is constructed.

***6.*** the best approximation is constructed - a function of class *C^0[a,b]* ,
using basis functions of finite element nodes

***7.*** Using a random number generator,
are specified on each finite element L> N internal “random
points."

***8.*** A system of linear algebraic equations is formed to determine the best fit element coefficients based on the values ​​of the function *y = f(x)* at “random points”,
generated in the previous paragraph.

***9.*** The system is solved in five ways (Conjugate gradient method, upper relaxation method, Cholesky method, LU, Gauss) written by hand and another from the library.

***10.*** The absolute and relative errors of the best approximations at “random points” using the norms of vectors l1,l2,linf

***11.*** At points of the segment [a,b] with a step h /100, the absolute and relative error of the best approximation using norms of vectors l1,l2,linf
