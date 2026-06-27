## But :

- [[Random Fileds and Geometry.pdf]]
- on souhaite modéliser les incertitudes liées aux observations que l'on fait lors d'un déplacement
- à chaque instant et chaque endroit, l'observation peut être vue comme une variable aléatoire

## Principe :

- un processus stochastique généralise la notion de v.a. en tant que "fonction
	aléatoire" _ie_ une famille $X$ de v.a. de $\mathbb{R}^d$ indexées par un
	ensemble $U$, (le temps $\mathbb{R}_+$, l'espace $\mathbb{R}^2$ ou
	$\mathbb{R}^3$, etc)
- un processus gaussien est un processus stochastique tel que toute combinaison
	linéaire d'éléments de $X$ suit une loi gaussienne (éventuellement
	multidimensionnelle) : $$
	\forall\:a\in\mathbb{R}^{(U)},\exists\:(m,\Sigma)\in \mathbb{R}^d\times\mathcal{M}_{d}(\mathbb{R})\:|\:\sum_{u\in U}a_uX(u)\sim\mathcal{N}(m,\Sigma)
	$$
- rappel : pour tout $X\sim\mathcal{N}(m,\Sigma)$ et $X'\sim\mathcal{N}(m',\Sigma')$, on a $X+X'\sim\mathcal{N}(m+m',\Sigma+\Sigma')$
- pour une combinaison linéaire d'un seul élément, on note alors $X(u) \sim\mathcal{N}(m(u),\Sigma(u,u))$
- on parlera de noyau de covariance $k$ pour désigner une fonction qui associe $\Sigma$ à un ensemble d'entrées $u$

## Notations :

- on note donc un processus gaussien comme $$X(\cdot)\sim\mathcal{GP}(m(\cdot),k(\cdot, \cdot\cdot))$$

- avec ces définitions, on a :
	- $m(u)=\mathbb{E}(X(u))$
	- $k(u,v)=\mathbb{E}\left((X(u)-m(u))(X(v)-m(v))^T\right)$

## Interprétation :

- la fonction de moyenne peut être vue comme une approximation peu précise du modèle d'observation
- son noyau de covariance comme une mesure du bruit d'observation (par exemple pour un bruit additif gaussien) ou du manque de précision de la moyenne par rapport à la réalité de terrain

## Remarques :

- les processus gaussiens sont une généralisation de la loi gaussienne multidimensionnelle car $$
\begin{bmatrix}
X(u_1)\\\vdots\\X(u_n)
\end{bmatrix}\sim\mathcal{N}\left(\begin{bmatrix}
m(u_1)\\\vdots\\m(u_n)
\end{bmatrix},\begin{bmatrix}
k(u_1,u_1)&\dots&k(u_1,u_n)\\\vdots&\ddots&\vdots\\k(u_n,u_1)&\dots&k(u_n,u_n)
\end{bmatrix}\right)
$$
- en dimensions supérieures à 1, $k$ n'est plus une fonction scalaire
- la matrice de covariance est donc une matrice par blocs