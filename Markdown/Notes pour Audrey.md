## Version 1

- $$
e(x):=\tilde{f}_{pl}(x)-\tilde{f}_p(x)=k(x,\tilde{X})\times g\times\tilde{Z}
$$
- $$\begin{align}
g:=&{\underbrace{\begin{bmatrix}A&B\\B^T&D\end{bmatrix}}_{G:=}}^{-1}-\begin{bmatrix}A^{-1}&0\\0&0\end{bmatrix}\\
=&\begin{bmatrix}A^{-1}BS^{-1}B^TA^{-1}&A^{-1}BS^{-1}\\-S^{-1}B^TA^{-1}&S^{-1}\end{bmatrix},&&S:=D-B^TA^{-1}B\\
=&\underbrace{\begin{bmatrix}-A^{-1}B\\I\end{bmatrix}}_{M:=}S^{-1}\begin{bmatrix}-A^{-1}B\\I\end{bmatrix}^T&&\phantom{S\:}=D-(A^{-1}B)^TA(A^{-1}B)
\end{align}$$
 - d'où :
 $$\begin{align}
|e(x)|=&|k(x,\tilde{X})\times MS^{-1}M^T\times\tilde{Z}|\\
\leq&||k(x,\tilde{X})||\times||M||^2\:||S^{-1}||\times||\tilde{Z}||\\
\leq&\tilde{N}\sigma_{max}\times(||A^{-1}B||^2+1)\:||S^{-1}||\times||\tilde{Z}||
\end{align}$$

## Version 2

- $$\begin{align}
\tilde{f}_{pl}(x):=&\sum_ik(x,\tilde{x}_i)\times(\alpha_{pl})_i=k(x,\tilde{X})\times\alpha_{pl}\\
\tilde{f}_p(x):=&\sum_ik(x,\tilde{x}_i)\times(\alpha_p)_i=k(x,\tilde{X})\times\alpha_p
\end{align}$$
- $$\begin{align}
\alpha_{pl}=&\begin{bmatrix}A^{-1}+A^{-1}BS^{-1}B^TA^{-1}&A^{-1}BS^{-1}\\-S^{-1}B^TA^{-1}&S^{-1}\end{bmatrix}\begin{bmatrix}\tilde{Z}_p\\\tilde{Z}_l\end{bmatrix}\\
=&\begin{bmatrix}\alpha_p\\0\end{bmatrix}+\underbrace{MS^{-1}M^T}_{\text{terme de correction ?}}\tilde{Z}
\end{align}$$

## Version 3

- écriture avec les RKHS fausse ?
- en tous cas pas directe :
$$\begin{align}
|f(x)-\tilde{f}(x)|=&|\langle f-\tilde{f}|k(x,\cdot)\rangle|\\
=&|\langle f-\tilde{f}|P_nk(x,\cdot)+(I-P_n)k(x,\cdot)\rangle|\\
=&|0+\langle f-\tilde{f}|(I-P_n)k(x,\cdot)\rangle|&\text{par orthogonalité}\\
\leq&||f-\tilde{f}||\times||k(x,\cdot)-P_nk(x,\cdot)||\\
\leq&||f-\tilde{f}||\times\sqrt{||k(x,\cdot)||^2-||P_nk(x,\cdot)||^2}\\
\leq&||f-\tilde{f}||\times\sqrt{k(x,x)-k(x,\tilde{X})G^{-1}k(\tilde{X,x})}\\
\end{align}$$
- le $k(x,x)$ ne vient pas de $f$ mais de $P_n$

## Pistes

- $S\preceq D\preceq G$
- pour $A\in SDP$, $||A||_{op}=\lambda_{max}(A)$ et $||A^{-1}||_{op}=1/\lambda_{min}(A)$
- [Inégalités de Weyl](https://fr.wikipedia.org/wiki/In%C3%A9galit%C3%A9s_de_Weyl) et [Théorème de Rayleigh](https://major-prepa.com/mathematiques/quotient-encadrement-rayleigh/)
- itérer en retirant 1 point par 1 point ($D$ scalaire, $B$ vecteur)
- [LMI](https://fr.wikipedia.org/wiki/In%C3%A9galit%C3%A9_matricielle_lin%C3%A9aire)
- formuler avec la matrice d'information de Fisher ?
- passer par un RKHS qui en contient un autre ?
