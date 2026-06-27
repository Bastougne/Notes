- [[Computing a Nearest Symmetric Positive Semidefinite Matrix.pdf]]

## Problématique :

- erreurs numériques lors du calcul des matrices de covariance
- elles doivent être réelles, symétriques et définies positives
- pour $M\in M_n(\mathbb{R})$, on cherche la matrice RSDP $\tilde M\in S^{++}_n(\mathbb{R})$ la plus proche de $M$ au sens de la norme de Frobenius $$
A\in M_{n\times m}(\mathbb{R})\mapsto||A||_F:=\sqrt{\sum_{i=1}^m \sum_{j=1}^n{A_{i,j}}^2}
$$
- on rappelle que la norme de Frobenius est la norme associée au produit scalaire canonique de $M_n(\mathbb{R})$ $$
(A,B)\in M_n(\mathbb{R})^2\mapsto \langle A,B\rangle:= \text{tr}(A^TB)$$ $$\forall A\in M_n(\mathbb{R}), ||A||_F=\sqrt{\langle A,A\rangle}
$$

## Principe :

- on commence par trouver la matrice (_semi_-définie) positive $\hat M\in S^{+}_n(\mathbb{R})$ la plus proche de $M$
- en dimension 1, le nombre positif le plus proche de $m\in \mathbb{R}$ est $\hat m\in\mathbb{R}_+:=\frac{1}{2}(m+|m|)$
- en dimension supérieure, la valeur absolue est remplacée par le terme symétrique positif d'une décomposition polaire
- on augmente légèrement les valeurs propres de $\hat M$ pour obtenir la matrice $\tilde M$ souhaitée

## Algorithme :

- symétrisation : $B:=\frac{1}{2}(M+M^T)$
- décomposition en valeurs singulières (toujours possible) : $B=UDV^T$
- calcul du terme symétrique positif d'une décomposition polaire : $P=VDV^T$
	- un terme orthogonal associé est $Q=UV^T$
	- démonstration : les 2 matrices sont dans les bons ensembles et :$$Q\times P=UV^T\times VDV^T=UDV^T=B$$
	- remarque : décomposition possible pour toute matricé réelle carrée, unicité de $P$ en général, unicité de $Q$ si $B$ est inversible
- calcul de $\hat M:=\frac{1}{2}(B+P)$
- calcul de $s:=\min(\text{sp}(\hat M))$
	- normalement, $s=0$, mais on le calcule quand-même en cas d'erreurs numériques
- correction des valeurs propres : $$
\varepsilon:=\left\{\begin{matrix}
-s+a\quad\text{avec}\quad a<<|s|&\text{si }s\leq0\\
0&\text{sinon}
\end{matrix}\right.
$$
- calcul de $\tilde M:=\hat M+\varepsilon I_n$
	- $\varepsilon I_n$ et $\hat M$ sont simultanément diagonalisables car elles commutent
	- $\text{sp}(\tilde M)=\{\lambda+\varepsilon|\lambda\in\text{sp}(\hat M)\}\subset\mathbb{R}_+^*$
- pour $a$ petit positif ou nul, $\tilde M$ est la matrice recherchée
