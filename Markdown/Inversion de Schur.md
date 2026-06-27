## Principe :

- [[Block Matrix Formulas.pdf]]
- on souhaite obtenir l'inverse d'une matrice par blocs inversible $$
M:=\begin{bmatrix}A&B\\C&D\end{bmatrix}
$$ sous la forme d'une matrice par blocs $$
M^{-1}=\begin{bmatrix}\bar{A}&\bar{B}\\\bar{C}&\bar{D}\end{bmatrix}
$$
- il existe 2 formules d'inversion par blocs basées sur les compléments de Schur $M/A$ et $M/D$

## Inversion selon $M/A$ :

- en réécrivant $MM^{-1}=I$, on trouve $$
\left\{\begin{align}
A\bar{A}+B\bar{C}&=I\\
A\bar{B}+B\bar{D}&=0\\
C\bar{A}+D\bar{C}&=0\\
C\bar{B}+D\bar{D}&=I\\
\end{align}\right.
$$
- en particulier, en supposant $A$ inversible, il vient de $(1)$ que $$
\bar{A}=A^{-1}(I-B\bar{C})
$$
- en remplaçant $\bar{A}$ par son expression dans $(3)$, on a alors $$
\begin{align}
&C\bar{A}+D\bar{C}=0\\
\iff\quad&CA^{-1}(I-B\bar{C})+D\bar{C}=0\\
\iff\quad&CA^{-1}+(D-CA^{-1}B)\bar{C}=0\\
\end{align}
$$
- posons $S:=D-CA^{-1}B$ le complément de Schur de $M$ par rapport à $A$ (aussi noté $M/A$), et supposons son inversibilité
- on peut alors réécrire $\bar{C}$ sous la forme $$
\bar{C}=-S^{-1}CA^{-1}
$$ puis $$
\bar{A}=A^{-1}(I-B\bar{C})=A^{-1}+A^{-1}BS^{-1}CA^{-1}
$$
- à partir de $(2)$, on déduit également que $$
\bar{B}=-A^{-1}B\bar{D}
$$
- en remplaçant $\bar{B}$ par son expression dans $(4)$, il vient $$
\begin{align}
&C\bar{B}+D\bar{D}=I\\
\iff\quad&-CA^{-1}B\bar{D}+D\bar{D}=I\\
\iff\quad&(D-CA^{-1}B)\bar{D}=I\\
\iff\quad&\bar{D}=S^{-1}\\
\end{align}
$$
- finalement $$
\bar{B}=-A^{-1}BS^{-1}
$$
- on a donc obtenu les expressions des blocs de $M^{-1}$ en fonction de ceux de $M$, qu'on résume ici $$
\left\{\begin{align}
S=&M/A=D-CA^{-1}B\\
\bar{A}=&A^{-1}+A^{-1}BS^{-1}CA^{-1}\\
\bar{B}=&-A^{-1}BS^{-1}\\
\bar{C}=&-S^{-1}CA^{-1}\\
\bar{D}=&S^{-1}\\
\end{align}\right.
$$

## Inversion selon $M/D$ :

- si on suppose que $D$ est inversible à la place de $A$, on peut refaire les calculs précédents de manière similaire
- la démonstration de l'inversion avec $M/D$ en connaissant celle avec $M/A$ est en fait équivalente à celle des [[Lemmes de Woodbury]]
- cette inversion s'écrit $$
\left\{\begin{align}
S=&M/D=A-BD^{-1}C\\
\bar{A}=&S^{-1}\\
\bar{B}=&-S^{-1}BD^{-1}\\
\bar{C}=&-D^{-1}CS^{-1}\\
\bar{D}=&D^{-1}+D^{-1}CS^{-1}BD^{-1}\\
\end{align}\right.
$$
- en particulier, si $A$ et $D$ sont inversibles, on peut combiner les 2 formes bloc par bloc