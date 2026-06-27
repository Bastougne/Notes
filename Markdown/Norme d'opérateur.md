## Définition :
- notons $||\cdot||$ la norme matricielle d'opérateur induite par la norme euclidienne $$
||M||:=||M||_{2,2}=\sup_{x\neq0}\left(\frac{||Mx||_2}{||x||_2}\right)
$$ qui vérifie donc $$
\forall x\in\mathbb{R^n},\:||Mx||_2\leq||M||\:||x||_2
$$

## Propriétés :
- pour $x\in\mathbb{R}^n$, $M_1\in M_{m\times p}(\mathbb{R})$ et $M_2\in M_{p\times n}(\mathbb{R})$, on a $$
||M_1M_2x||_2\leq||M_1||\:||M_2x||_2\leq||M_1||\:||M_2||\:||x||_2
$$
- la norme de d'opérateur est donc sous-multiplicative, _i.e._ $$
||M_1M_1||\leq||M_1||\:||M_2||
$$
- de même $$
||(M_1+M_2)x||_2^2\leq||M_1x||_2^2+||M_2x||_2^2=(||M_1||^2+||M_2||^2)\:||x||_2^2
$$
- la norme de d'opérateur vérifie l'inégalité triangulaire, _i.e._ $$
||M_1+M_1||\leq\sqrt{||M_1||^2+||M_2||^2}\leq||M_1||+||M_2||
$$
- enfin, décomposons $M\in M_{m\times n}(\mathbb{R})$ et $x$ par blocs $$
M=\begin{bmatrix}M_{1,1}&\cdots&M_{1,\tilde n}\\\vdots&\ddots&\vdots\\M_{\tilde m,1}&\cdots&M_{\tilde m,\tilde n}\end{bmatrix}\quad\text{et}\quad x=\begin{bmatrix}x_1\\\vdots\\x_{\tilde n}\end{bmatrix}
$$
- on a $$
\begin{align}
||Mx||_2^2=&\left|\left|\begin{bmatrix}M_{1,1}&\cdots&M_{1,\tilde n}\\\vdots&\ddots&\vdots\\M_{\tilde m,1}&\cdots&M_{\tilde m,\tilde n}\end{bmatrix}\begin{bmatrix}x_1\\\vdots\\x_{\tilde n}\end{bmatrix}\right|\right|_2^2\\
=&\left|\left|\sum_{i=1}^{\tilde m}\sum_{j=1}^{\tilde n}M_{i,j}x_j\right|\right|_2\\
\leq&\sum_{i=1}^{\tilde m}\sum_{j=1}^{\tilde n}||M_{i,j}x_j||_2^2\\
\leq&\left(\sum_{i=1}^{\tilde m}\sum_{j=1}^{\tilde n}||M_{i,j}||^2\right)\max_{j\in[\![1,\tilde n]\!]}(||x_j||_2^2)\\
\leq&\left(\sum_{i=1}^{\tilde m}\sum_{j=1}^{\tilde n}||M_{i,j}||^2\right)||x||_2^2\\
\end{align}
$$
- et donc $$
||M||\leq\sqrt{\sum_{i=1}^{\tilde m}\sum_{j=1}^{\tilde n}||M_{i,j}||^2}\leq\sum_{i=1}^{\tilde m}\sum_{j=1}^{\tilde n}||M_{i,j}||
$$

## Majorations utiles :
- pour $x\in\mathbb{R}^n$ et $M\in M_{m\times n}(\mathbb{R})$, on a $$
\begin{align}
||Mx||^2_2=&\left|\left|\begin{bmatrix}
\displaystyle\sum_{j=1}^nM_{1,j}x_j&\dots&\displaystyle\sum_{j=1}^nM_{m,j}x_j
\end{bmatrix}^T\right|\right|^2_2\\
=&\sum_{i=1}^m\left(\sum_{j=1}^nM_{i,j}x_j\right)^2\\
\leq&\sum_{i=1}^m\left(\sum_{j=1}^nM_{i,j}^2\right)\left(\sum_{j=1}^nx_j^2\right)\\
\leq&||M||^2_F\:||x||^2_2\\
\end{align}
$$ par inégalité de Cauchy-Schwarz
- d'où la norme de Frobenius domine la norme d'opérateur $$
||M||\leq||M||_F
$$
- en particulier, si $M$ est bornée termes à termes par une matrice $A$, alors $$
||M||\leq||M||_F\leq||A||_F
$$
