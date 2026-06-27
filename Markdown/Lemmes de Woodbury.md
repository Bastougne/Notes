## But :
- [[Matrix Inversion Identities.pdf]]
- on souhaite inverser des matrices s'écrivant sous la forme $A+BDC$
- si $D$ est de dimension inférieure à $A$, et $A$ est facilement inversible, l'inversion peut être plus rapide que par un pivot de Gauss usuel
- notamment utile dans des cas de réduction de rang ou d'inversion par blocs

## Lemme 1 :

- pour $P\in\mathcal{M}_n(\mathbb{R})$ et $I$ la matrice identité de taille $n$, on a $$
\begin{align}
(I+P)^{-1}=&(I+P)^{-1}(I+P-P)\\
=&(I+P)^{-1}(I+P)-(I+P)^{-1}P\\
=&I-(I+P)^{-1}P
\end{align}
$$

## Lemme 2 :

- pour $Q\in\mathcal{M}_n(\mathbb{R})$, on a $$
P+PQP=P(I+QP)=(I+PQ)P
$$
- en inversant $I+QP$ et $I+PQ$ dans la dernière égalité, on trouve alors $$
\begin{align}
P(I+QP)=&(I+PQ)P\\
P=&(I+PQ)P(I+QP)^{-1}\\
(I+PQ)^{-1}P=&P(I+QP)^{-1}\\
\end{align}
$$

## Lemme de Woodbury classique :

- soient $A,B,C,D$ des matrices telles que $A+BDC$ ait un sens
- on suppose $A$ et $D$ inversibles, d'où $$
\begin{align}
(A+BDC)^{-1}=&\left(A(I+A^{-1}BDC)\right)^{-1}\\
=&(I+A^{-1}BDC)^{-1}A^{-1}\\
=&(I-(I+A^{-1}BDC)^{-1}A^{-1}BDC)A^{-1}&\text{d'après (1)}\\
=&A^{-1}-(I+(A^{-1}B)DC)^{-1}(A^{-1}B)DCA^{-1}\\
=&A^{-1}-A^{-1}B(I+DCA^{-1}B)^{-1}DCA^{-1}&\text{d'après (2)}\\
=&A^{-1}-A^{-1}B(D^{-1}+CA^{-1}B)^{-1}CA^{-1}\\
\end{align}
$$

## Lemmes de Woodbury factorisés :

- avec les mêmes hypothèses, on écrit $$
\begin{align}
(A+BDC)^{-1}BD=&((I+BDCA^{-1})A)^{-1}BD\\
=&A^{-1}(I+BDCA^{-1})^{-1}BD\\
=&A^{-1}B(I+DCA^{-1}B)^{-1}D&\text{d'après (2)}\\
=&A^{-1}B(D^{-1}+CA^{-1}B)^{-1}\\
\end{align}
$$
- et de même $$
\begin{align}
DC(A+BDC)^{-1}=&DC(A(I+A^{-1}BDC))^{-1}\\
=&DC(I+A^{-1}BDC)^{-1}A^{-1}\\
=&D(I+CA^{-1}BD)^{-1}CA^{-1}&\text{d'après (2)}\\
=&(D^{-1}+CA^{-1}B)^{-1}CA^{-1}\\
\end{align}
$$
