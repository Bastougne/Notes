## But :
- [[Filtrage Bayésien et Approximation Particulaire Méthodes de Monte Carlo.pdf]]
- on souhaite calculer l'espérance $$
\mathbb{E}_p(\varphi(X))=\int\varphi(x)\:p(x)\:dx
$$ d'une variable aléatoire $X$ de densité de probabilité $p$ de variance finie
- en notation de théorie de la mesure, pour une mesure $\mu$ telle que $$
\mu(dx)=p(x)\:dx
$$ on cherche à calculer $$
\langle\mu,\varphi\rangle:=\int\varphi(x)\:\mu(dx)
$$

## Mesure empirique :
- la mesure $\mu$ étant continue, le calcul de l'intégrale est difficile en toute généralité
- pour $N\in\mathbb{N}^*$, on construit $X^{(1)},...,X^{(N)}$ des variables aléatoires indépendantes et identiquement distribuées selon la mesure $\mu$
- on approche alors $\mu$ par $$
S^N_\delta(\mu):=\frac{1}N\sum_{i=1}^N\delta_{X^{(i)}}
$$ la distribution empirique associée aux $\left(X^{(i)}\right)_{i\in[\![1,N]\!]}$
- ainsi, pour toute fonction mesurable bornée $\phi$, on a $$
\langle\mu,\phi\rangle\approx\langle S^N_\delta(\mu),\phi\rangle=\frac{1}N\sum_{i=1}^N\phi\left(X^{(i)}\right)
$$
- c'est-à-dire $$
\mu\approx S^N_\delta(\mu)=\frac{1}N\sum_{i=1}^N\delta_{X^{(i)}}
$$

## Estimateur de Monte-Carlo :
- l'estimateur de Monte-Carlo de $\mathbb{E}_p(\varphi(X))$ est alors $$
\langle S^N_\delta(\mu),\varphi\rangle=\frac{1}N\sum_{i=1}^N\varphi\left(X^{(i)}\right)
$$
- c'est un estimateur non-biaisé car $$
\begin{align}
\text{Biais}(\langle S^N_\delta(\mu),\varphi\rangle)&=\mathbb{E}(\langle S^N_\delta(\mu),\varphi\rangle)-\langle\mu,\varphi\rangle\\
&=\mathbb{E}\left(\frac{1}N\sum_{i=1}^N\varphi\left(X^{(i)}\right)\right)-\langle\mu,\varphi\rangle\\
&=\frac{1}N\sum_{i=1}^N\mathbb{E}\left(\varphi\left(X^{(i)}\right)\right)-\langle\mu,\varphi\rangle\\
&=\frac{1}N\sum_{i=1}^N\left(\int\varphi\left(X^{(i)}(\omega)\right)\:\mathbb{P}(d\omega)\right)-\langle\mu,\varphi\rangle\\
&=\frac{1}N\sum_{i=1}^N\left(\int\varphi(x)\:\mathbb{P}_{X^{(i)}}(dx)\right)-\langle\mu,\varphi\rangle\\
&=\frac{1}N\sum_{i=1}^N\left(\int\varphi(x)\:\mu(dx)\right)-\langle\mu,\varphi\rangle\\
&=\frac{1}N\sum_{i=1}^N\left(\langle\mu,\varphi\rangle\right)-\langle\mu,\varphi\rangle\\
&=0
\end{align}
$$ par linéarité de l'espérance et théorème de transfert
- notons $$
\mathbb{V}_\mu(\phi):=\langle\mu,|\phi-\langle\mu,\phi\rangle|^2\rangle=\langle\mu,|\phi|^2\rangle-|\langle\mu,\phi\rangle|^2
$$
- on a alors $$
\begin{align}
\mathbb{E}\left(|\langle S^N_\delta(\mu)-\mu,\varphi\rangle|^2\right)&=\mathbb{E}\left(|\langle S^N_\delta(\mu),\varphi\rangle-\langle\mu,\varphi\rangle|^2\right)\\
&=\mathbb{E}\left(\left|\frac{1}N\sum_{i=1}^N\left(\varphi\left(X^{(i)}\right)-\langle\mu,\varphi\rangle\right)\right|^2\right)\\
&=\frac{1}{N^2}\sum_{i=1}^N\mathbb{E}\left(\left|\varphi\left(X^{(i)}\right)-\langle\mu,\varphi\rangle\right|^2\right)\\
&=\frac{1}{N^2}\sum_{i=1}^N\langle\mu,|\varphi-\langle\mu,\varphi\rangle|^2\rangle\\
&=\frac{1}{N}\mathbb{V}_\mu(\varphi)
\end{align}
$$ par indépendance des $X^{(i)}$

## Théorème central limite et loi forte des grands nombres :
- par théorème central limite $$
\langle S^N_\delta(\mu),\varphi\rangle\xrightarrow{\mathcal{L}}\langle\mu,\varphi\rangle
$$ et $$\sqrt{N}\langle S^N_\delta(\mu)-\mu,\varphi\rangle\xrightarrow{\mathcal{L}}\mathcal{N}(0,\mathbb{V}_\mu(\varphi))
$$
- de plus, par loi forte des grands nombres $$
\langle S^N_\delta(\mu),\varphi\rangle\xrightarrow{\text{p.s.}}\langle\mu,\varphi\rangle
$$
- l'estimateur de Monte-Carlo converge donc presque sûrement vers l'espérance théorique à une vitesse de $\sqrt{N}$