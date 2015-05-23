Reproducing our results
-----------------------

Please run myDriver.m to reproduce the results.

Note
----

In an earlier version of our code, we had achieved the results you can find in the Images/old-values/ folder. Unfortunately, we lost the code that produced those results - it came from initializing the mean estimate in a particular way. The latest version of our code will reproduce the results you can find in the Images/ folder.

Basically, the algorithm now is converging towards a different extremum than before.

Images produced
---------------

"initial-pointset.png" shows the pointsets before running the alignment algorithm.
"converged-pointset.png" shows the aligned pointsets after convergence.
"computed-shape-mean.png" shows the computed shape mean after convergence.
"eigenvalues.png" is a bar plot of the eigenvalues.

The images named "mode<number>.png" show the variation of the mean shape along the <number> largest eigenvalue. For example, "mode1.png" shows the mean shape shifted along the direction of the eigenvector corresponding to the largest eigenvalue, with the coefficient of the eigenvector being set to -2, -1, 0, 1 and 2 times the square root of the largest eigenvalue respectively.

Team
----

Ashray Malhotra
Rohan Prinja