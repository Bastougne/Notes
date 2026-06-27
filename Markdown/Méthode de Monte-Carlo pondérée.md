## But :
- amélioration de la [[Méthode de Monte-Carlo]]
- on cherche toujours à calculer l'espérance $$
\mathbb{E}_p(\varphi(X))=\int\varphi(x)\:p(x)\:dx=\langle\mu,\varphi\rangle
$$ d'une variable aléatoire de densité $p$ vérifiant $\mu(dx)=p(x)\:dx$
- on ne peut pas appliquer une méthode de Mont-Carlo usuelle car l'échantillonnage depuis $p$ est difficile

## Mesure empirique pondérée :
- supposons qu'on connaisse une loi $q$ (dite d'importance, ou _proposal_ en anglais) de densité plus simple à échantillonner, et qui vérifie $\nu(dx)=q(x)\:dx$
- on aurait alors $$
p(x)=\frac{p(x)}{q(x)}\:q(x)
$$ et donc $$
\mathbb{E}_p(\varphi(X))=\int\left(\varphi(x)\:\frac{p(x)}{q(x)}\right)\:q(x)\:dx=\mathbb{E}_q\left(\varphi(X)\:\frac{p(X)}{q(X)}\right)
$$
- en notation de théorie de la mesure, cela revient à écrire $$
\langle\mu,\varphi\rangle=\langle\nu,\varphi\frac{d\mu}{d\nu}\rangle
$$
- en construisant $X^{(1)},...,X^{(N)}$ des variables aléatoires indépendantes et identiquement distribuées selon la mesure $\nu$, on approche $\mu$ par $$
S^N_\delta(\mu,\nu):=\frac{d\mu}{d\nu}S^N_\delta(\nu)=\frac{1}N\sum_{i=1}^N\frac{d\mu\left(X^{(i)}\right)}{d\nu\left(X^{(i)}\right)}\delta_{X^{(i)}}
$$ la distribution empirique pondérée par $\displaystyle{\frac{d\mu}{d\nu}}$

## Estimateur de Monte-Carlo pondéré :
- l'estimateur de Monte-Carlo pondéré est donné par $$
\langle S^N_\delta(\mu,\nu),\varphi\rangle=\langle S^N_\delta(\nu),\frac{d\mu}{d\nu}\varphi\rangle=\frac{1}N\sum_{i=1}^N\frac{d\mu\left(X^{(i)}\right)}{d\nu\left(X^{(i)}\right)}\varphi\left(X^{(i)}\right)
$$
- par méthode de Monte-Carlo, on a que $\langle S^N_\delta(\mu,\nu),\varphi\rangle$ est non-biaisé et que $$
\mathbb{E}\left(|\langle S^N_\delta(\mu,\nu)-\mu,\varphi\rangle|^2\right)=\frac{1}{N}\mathbb{V}_\mu\left(\frac{d\mu}{d\nu}\varphi\right)
$$

## Lois proportionnelles :
- cet estimateur n'est pas normalisé en général, car $$
\langle S^N_\delta(\mu,\nu),1\rangle=\frac{1}N\sum_{i=1}^N\frac{d\mu\left(X^{(i)}\right)}{d\nu\left(X^{(i)}\right)}
$$
- en outre, on ne connait parfois que $$
\tilde\mu=Z_\mu\:\mu\quad\text{et}\quad\tilde\nu=Z_\nu\:\nu
$$ ou de manière équivalente $$
\tilde{p}=Z_\mu\:p\quad\text{et}\quad\tilde{q}=Z_\nu\:q
$$
- on a alors $$
\langle S^N_\delta(\mu,\nu),1\rangle=\frac{1}N\sum_{i=1}^N\frac{Z_\mu}{Z_\nu}\frac{d\tilde\mu\left(X^{(i)}\right)}{d\tilde\nu\left(X^{(i)}\right)}
$$
- on peut donc réécrire l'estimateur de Monte-Carlo pondéré sous la forme $$
\langle S^N_\delta(\mu,\nu),\varphi\rangle=\sum_{i=1}^N w^{(i)}\varphi\left(X^{(i)}\right)
$$ avec les poids pondérés $w^{(i)}$ et non pondérés $r^{(i)}$ respectivement définis par $$
w^{(i)}:=\frac{r^{(i)}}{\langle S^N_\delta(\mu),1\rangle}=\frac{r^{(i)}}{\displaystyle\sum_{j=1}^N r^{(j)}}
$$ et $$
r^{(i)}:=\frac{d\tilde\mu\left(X^{(i)}\right)}{d\tilde\nu\left(X^{(i)}\right)}=\frac{\tilde{p}\left(X^{(i)}\right)}{\tilde{q}\left(X^{(i)}\right)}
$$
- la distribution empirique pondérée s'écrit donc $$
S^N_\delta(\mu,\nu)=\sum_{i=1}^N w^{(i)}\delta_{X^{(i)}}
$$
