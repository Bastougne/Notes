## Estimateur de l'espérance a posteriori

- plaçons-nous dans l'espace $\mathcal{L}^2(\mathbb{R}^p)$ muni de la distance associée au produit scalaire canonique des variables aléatoire $$
(A,B)\in\mathcal{L}^2(\mathbb{R}^p)\mapsto\langle A,B\rangle:=\text{tr}\left(\mathbb{E}(A^TB)\right)$$
- dans $(\mathcal{L}^2(\mathbb{R}^p),\langle\cdot|\cdot\rangle)$, le problème de minimisation revient à résoudre $$
\underset{\psi\in\mathcal{L}^2(\mathbb{R}^p)\:|\:\mathbb{E}(|\psi(Y)|^2)<+\infty}{\text{argmin}}(\langle\psi(Y)-\varphi(X)|\psi(Y)-\varphi(X)\rangle^2)
$$ _i.e._ $$
\underset{\psi\in\mathcal{L}^2(\mathbb{R}^p)\:|\:\mathbb{E}(|\psi(Y)|^2)<+\infty}{\text{argmin}}(\mathbb{E}(\psi(Y)-\varphi(X))^2)
$$
- un tel estimateur est appelé estimateur MMSE (de l'anglais _minimum mean squared errror_)
- posons $\hat\psi(y):=\mathbb{E}(\varphi(X)|Y=y)$
- on a, pour $\psi\in\mathcal{L}^2(\mathbb{R}^p)$ tel que $\mathbb{E}(|\psi(Y)|^2)<+\infty$ $$
\begin{align}
\text{tr}\left(\mathbb{E}\left((\psi(Y)-\varphi(X))^T(\psi(Y)-\varphi(X))\right)\right)=&\text{tr}\left(\mathbb{E}\left((\psi(Y)-\hat\psi(Y)+\hat\psi(Y)-\varphi(X))^T\right.\right.\\
&\left.\left.\phantom{\text{tr}(\:\:\mathbb{E}(\:\:}(\psi(Y)-\hat\psi(Y)+\hat\psi(Y)-\varphi(X))\right)\right)\\
\\
=&\text{tr}\left(\mathbb{E}\left((\psi(Y)-\hat\psi(Y))^T(\psi(Y)-\hat\psi(Y))\right)\right)+\\
&\text{tr}\left(\mathbb{E}\left((\psi(Y)-\hat\psi(Y))^T(\hat\psi(Y)-\varphi(X))\right)\right)+\\
&\text{tr}\left(\mathbb{E}\left((\hat\psi(Y)-\varphi(X))^T(\psi(Y)-\hat\psi(Y))\right)\right)+\\
&\text{tr}\left(\mathbb{E}\left((\hat\psi(Y)-\varphi(X))^T(\hat\psi(Y)-\varphi(X))\right)\right)\\
\end{align}
$$
- or $$
\begin{align}
\text{tr}\left(\mathbb{E}\left((\psi(Y)-\hat\psi(Y))^T(\hat\psi(Y)-\varphi(X))\right)\right)&=\int\int(\psi(y)-\hat\psi(y))^T(\hat\psi(y)-\varphi(x))\:p(dx,dy)\\
&=\int\int(\psi(y)-\hat\psi(y))^T(\hat\psi(y)-\varphi(x))\:p(dx|y)\:p(dy)\\
&=\int(\psi(y)-\hat\psi(y))^T\left(\int(\hat\psi(y)-\varphi(x))\:p(dx|y)\right)\:p(dy)\\
&=\int(\psi(y)-\hat\psi(y))^T\times0\:p(dy)\\
&=0
\end{align}
$$ par définition de $\hat\psi(Y)$, et de même pour $\text{tr}\left(\mathbb{E}\left((\hat\psi(Y)-\varphi(X))^T(\psi(Y)-\hat\psi(Y))\right)\right)$
- finalement $$
\begin{align}
\text{tr}\left(\mathbb{E}\left((\psi(Y)-\varphi(X))^T(\psi(Y)-\varphi(X))\right)\right)=&\text{tr}\left(\mathbb{E}\left((\psi(Y)-\hat\psi(Y))^T(\psi(Y)-\hat\psi(Y))\right)\right)+\\
&\text{tr}\left(\mathbb{E}\left((\hat\psi(Y)-\varphi(X))^T(\hat\psi(Y)-\varphi(X))\right)\right)\\
\geq&\text{tr}\left(\mathbb{E}\left((\hat\psi(Y)-\varphi(X))^T(\hat\psi(Y)-\varphi(X))\right)\right)
\end{align}
$$
- la moyenne _a posteriori_ $\hat\psi(Y)=\mathbb{E}(\varphi(X)|Y)$ est donc l'estimateur MMSE de $\varphi(X)$ connaissant $Y$
- $\hat\psi(Y)$ est en fait la projection au sens de $\langle\cdot|\cdot\rangle$ de $\varphi(X)$ sur l'espace des variables aléatoires qui dépendent de $Y$

## Estimateur de la médiane a posteriori

- l'estimateur MMSE est très sensible aux valeurs extrêmes
- il n'est donc pas adapté aux densités multimodales ou à queues lourde
- plaçons-nous dans $\mathcal{L}^2(\mathbb{R}^p)$ muni de la distance (non-euclidienne) $$
(A,B)\in\mathcal{L}^1(\mathbb{R}^p)\mapsto d(A,B):=\mathbb{E}(|A-B|)
$$
- le problème de minimisation devient $$
\underset{\psi\in\mathcal{L}^1(\mathbb{R}^p)\:|\:\mathbb{E}(|\psi(Y)|)<+\infty}{\text{argmin}}(\mathbb{E}(|\psi(Y)-\varphi(X)|))
$$ qu'on réécrit en $$
\underset{m\in\mathbb{R}}{\text{argmin}}(\mathbb{E}(|\varphi(X)-m|))=\underset{m\in\mathbb{R}}{\text{argmin}}(\int|\varphi(x)-m|\:p(x)\:dx))
$$
- en dérivant cette expression, il vient que $m$ doit satisfaire $$
\begin{align}
&\frac{d}{dm}\left(\int|\varphi(x)-m|\:p(x)\:dx\right)=0\\
\iff\quad&\frac{d}{dm}\left(\int\left((\varphi(x)-m)\:\mathbb{1}_{\varphi(x)\geq m}(x)-(\varphi(x)-m)\:\mathbb{1}_{\varphi(x)<m}(x)\right)\:p(x)\:dx\right)=0\\
\iff\quad&\int\left(-\mathbb{1}_{\varphi(x)\geq m}(x)+\mathbb{1}_{\varphi(x)<m}(x)\right)\:p(x)\:dx=0\\
\iff\quad&\mathbb{E}\left(\mathbb{1}_{\varphi(X)\geq m}\right)=\mathbb{E}\left(\mathbb{1}_{\varphi(X)<m}\right)\\
\iff\quad&\mathbb{P}(\varphi(X)\geq m)=\mathbb{P}(\varphi(X)<m)\\
\end{align}
$$
- il s'agit de la médiane _a posteriori_
- cet estimateur est moins sensible aux valeurs extrêmes que le MMSE

## Estimateur du maximum a posteriori

- définit par $$
\underset{x\in\mathbb{R}^p}{\text{argmax}}(\mathbb{P}(X=x|Y=y))=\underset{x\in\mathbb{R}^p}{\text{argmax}}(\mathbb{P}(Y=y|X=x)\:\mathbb{P}(X=x))
$$
- très sensible aux multimodalités

## Estimateur du maximum de vraisemblance

- définit par $$
\underset{x\in\mathbb{R}^p}{\text{argmax}}(\mathbb{P}(Y=y|X=x))
$$
- cet estimateur n'est pas bayésien car il ne prend pas en compte l'information a priori