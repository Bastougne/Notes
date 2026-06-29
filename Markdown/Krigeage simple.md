## But :

- méthode de régression par [[Processus Gaussiens]] (_Gaussian process regression_ ou GPR en anglais)
- c'est une méthode qui garantit de trouver un estimateur linéaire non-biaisé à minimum de variance pour l'image d'une fonction bruitée et partiellement connue sous réserve qu'elle suive un processus gaussien
![[Étapes du krigeage.pdf]]

## Références :

- [[Gaussian Processes for Machine Learning.pdf]]
- [[Le Krigeage Revue de la Théorie et Application à l'Interpolation Spatiale de Données de Précipitations.pdf]]
- supports de cours : [[Krigeage.pdf]] [[Introduction aux Statistiques Spatiales le Krigeage.pdf]]
- [[Gaussian Processes for Regression a Quick Introduction.pdf]] ajoute un terme de bruit dans la définition du noyau de covariance pour simplifier les notations : $\sigma_n^2 \: \delta(x'-x)$
- ce n'est possible que si les incertitudes sur $z$ et $\tilde z$ sont les mêmes (donc les matrices de covariance associées aux capteurs ayant permit l'acquisition des deux jeux de données)

## Notations :

- pour $f\sim\mathcal{GP}(m,k):\mathbb{R}^d\rightarrow\mathbb{R}$, on souhaite estimer $z:=f(x)+\varepsilon$ connaissant des mesures connues $\tilde{Z}:=\begin{bmatrix}f(\tilde{x}_1)+\tilde\varepsilon_1&\dots&f(\tilde{x}_n)+\tilde\varepsilon_n\end{bmatrix}^T$, bruitées par des bruits blancs gaussiens de variances connues $$
\varepsilon\sim\mathcal{N}(0,R)\quad\text{et}\quad\forall i\in[\![1,n]\!],\:\tilde\varepsilon_i\sim\mathcal{N}(0,\tilde{R}_i)
$$
- on suppose ici que la fonction de moyenne $m$ du processus gaussien est une constante connue
- en posant $f_{centré}=f-m$, on se ramène donc à un processus gaussien de moyenne nulle
- l'estimateur souhaité est de la forme $$
z^*=\sum_{i=1}^{\tilde{n}}\lambda_i\tilde{z}_i=\Lambda^T\tilde{Z}
$$ qui est bien linéaire et sans biais car $$
\mathbb{E}(z^*-z)=0-\sum_{i=1}^{\tilde{n}}\lambda_i\times0
$$

## Coefficients du krigeage simple :

- on souhaite minimiser la variance $$
\begin{align}
\text{Cov}(z^*-z)=&\sum_{i=1}^{\tilde{n}}\sum_{j=1}^{\tilde{n}}\lambda_i\lambda_jk(\tilde{x}_i,\tilde{x}_j)+\sum_{i=1}^{\tilde{n}}\lambda_i^2\tilde{R}-2\sum_{i=1}^{\tilde{n}}\lambda_ik(\tilde{x}_i,x)+k(x,x)+R\\
=&\Lambda^T\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)\Lambda-2\Lambda^Tk(\tilde{X},x)+k(x,x)+R
\end{align}
$$
- en dérivant par rapport à $\Lambda$, il vient $$
\begin{align}
&\frac{d}{d\Lambda}\left(\Lambda^T\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)\Lambda-2\Lambda^Tk(\tilde{X},x)+k(x,x)+R\right)=0\\
\iff\quad&2\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)\Lambda-2k(\tilde{X},x)=0\\
\iff\quad&\Lambda=\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)^{-1}k(\tilde{X},x)\\
\end{align}
$$

## Estimateur du krigeage simple :

- l'estimateur du krigeage simple s'écrit donc $$
z^*_s(x,\tilde{X},\tilde{Z})=\Lambda^T\tilde{Z}=k(x,\tilde{X})\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)^{-1}\tilde{Z}
$$
- dans le cas d'une moyenne connue $m$, on a alors $$
z^*_s(x,\tilde{X},\tilde{Z},m)=m+k(x,\tilde{X})\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)^{-1}(\tilde{Z}-m\textbf{1})
$$
- la matrice de covariance du krigeage simple est la covariance minimisée $$
\begin{align}
R^*_s(x,\tilde{X})=&\Lambda^T\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)\Lambda-2\Lambda^Tk(\tilde{X},x)+k(x,x)+R\\
=&k(x,\tilde{X})\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)^{-1}\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)k(\tilde{X},x)-\\
&2k(x,\tilde{X})\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)^{-1}k(\tilde{X},x)+k(x,x)+R\\
=&k(x,x)+R-k(x,\tilde{X})\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)^{-1}k(\tilde{X},x)\\
\end{align}
$$

## Espérance _a posteriori_ :

- puisque $z$ et les éléments de $\tilde{Z}$ suivent tous le $f$ bruité, leur loi jointe suit une loi normale multivariée $$
\begin{bmatrix}
z\\\tilde{Z}
\end{bmatrix}\sim\mathcal{N}\left(\begin{bmatrix}
0\\0\end{bmatrix},\begin{bmatrix}
k(x,x)+RI_n&k(x,\tilde{X})\\k(x,\tilde{X})^T&k(\tilde{x},\tilde{x})+\tilde{R}I_{\tilde{n}}
\end{bmatrix}\right)
$$
- on va inverser la matrice de covariance de la loi jointe en utilisant l'[[Inversion de Schur]] du bloc inférieur droit, et en posant $$
\left\{\begin{align}
A:=&k(x,x)+RI_n\\
B:=&k(x,\tilde{X})\\
C:=&k(x,\tilde{X})\\
D:=&k(\tilde{x},\tilde{x})+\tilde{R}I_{\tilde{n}}\\
S:=&A-BD^{-1}C\\
\end{align}\right.
$$
- en remarquant que $D=D^T$ et $C=B^T$ dans la matrice de covariance de la loi jointe, posons $$
\begin{align}
\alpha:=&\begin{bmatrix}z\\\tilde{Z}\end{bmatrix}^T\begin{bmatrix}A&B\\C&D\end{bmatrix}^{-1}\begin{bmatrix}z\\\tilde{Z}\end{bmatrix}\\
=&\begin{bmatrix}z\\\tilde{Z}\end{bmatrix}^T\begin{bmatrix}\bar{A}&\bar{B}\\\bar{C}&\bar{D}\end{bmatrix}\begin{bmatrix}z\\\tilde{Z}\end{bmatrix}\\
=&z^T\bar{A}z+z^T\bar{B}\tilde{Z}+\tilde{Z}^T\bar{C}z+\tilde{Z}^T\bar{D}\tilde{Z}\\
=&z^TS^{-1}z+z^TS^{-1}(-BD^{-1}\tilde{Z})+(-BD^{-1}\tilde{Z})^TS^{-1}z+\\
&(-BD^{-1}\tilde{Z})^TS^{-1}(-BD^{-1}\tilde{Z})+\tilde{Z}^TD^{-1}\tilde{Z}\\
=&(z-BD^{-1}\tilde{Z})^TS^{-1}(z-BD^{-1}\tilde{Z})+\tilde{Z}^TD^{-1}\tilde{Z}
\end{align}
$$
- la loi conditionnelle $p(z|\tilde{Z})$ s'écrit alors $$
\begin{align}
p(z|\tilde{Z})=&\frac{p(z,\tilde{Z})}{p(\tilde{Z})}\\
\propto&\frac{\text{exp}\left(-\frac{1}{2}\alpha\right)}{\text{exp}\left(-\frac{1}{2}\tilde{Z}^TD^{-1}\tilde{Z}\right)}\\
\propto&-\frac{1}{2}\left((z-BD^{-1}\tilde{Z})^TS^{-1}(z-BD^{-1}\tilde{Z})\right)
\end{align}
$$
- l'estimateur du krigeage simple est donc également l'espérance _a posteriori_ $$
z^*=k(x,\tilde{X})\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)^{-1}\tilde{Z}=\mathbb{E}(z|\tilde{Z})
$$ et leur matrice de covariance sont égales (et égales au complément de Schur inférieur droit de la matrice de covariance de la loi jointe)

## Estimations simultanées :

- on peut estimer simultanément $Z:=\begin{bmatrix}f(x_1)+\varepsilon_1&\dots&f(x_n)+\varepsilon_n\end{bmatrix}^T$ en connaissant les mêmes échantillons bruités $\tilde{Z}$
- le vecteur des poids $\Lambda$ devient alors une matrice, et l'estimateur du krigeage devient un vecteur colonne correspondant à la concaténation des estimateurs individuels en chaque points
- la matrice de covariance associée est telle que ses éléments diagonaux sont exactement les covariances obtenues pour chaque estimation ponctuelle

## Estimation de $\mathbb{R}^d\rightarrow\mathbb{R}^{d'}$ :

- les entrées du processus Gaussien peuvent être de dimension finie quelconque $d$, mais ses sorties sont scalaires
- $k$ et $m$ sont donc des fonctions de $\mathbb{R}^d\rightarrow\mathbb{R}$
- il faut donc effectuer une régression pour chaque dimension $i\in[\![1,d']\!]$ des sorties souhaitées
- il semble difficile d'introduire des relations entre les dimensions des sorties préalablement au krigeage
