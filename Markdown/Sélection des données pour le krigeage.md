## But :

- [[Computationally Efficient Data Selection in Gaussian Process Regression-Based Particle Filters for Autonomous Navigation.pdf]]
- les multiplications matricielles et l'inversion de la matrice de Gram $$G=k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}$$ sont des lentes si $\tilde{n}$ est grand
- supposons que le noyau de covariance du krigeage est isotrope, _i.e._ $$
\exists f:\mathbb{R}_+\rightarrow\mathbb{R}\:|\:\forall x,x'\in(\mathbb{R}^n)^2,k(x,x')=f(||x-x'||)
$$ avec $f$ décroissante
- les échantillons suffisamment éloignés des particules ralentissent donc les calculs sans améliorer significativement la précision de la régression, et on peut ne pas les prendre en compte pour gagner en rapidité dans les calculs
- en notant $\{\tilde{x}_i,i\in[\![1,\tilde{n}]\!]\}$ l'ensemble des données d'entraînement, on va déterminer des ensembles $S$ tels que les données sélectionnées pour le krigeage seront les $\{\tilde{x}_i,i\in[\![1,\tilde{n}]\!]\}\cap S$ 

## Bornes d'erreur :

- [[Bornes d'erreur du krigeage adaptatif]]
- [[Notes pour Audrey]]
- [[Bornes de krigeage dans L2]]

## Fenêtrage autour de l'état prédit :

- présenté dans [[Adaptive Kriging Particle Filter and its Application to Terrain Aided Navigation.pdf]] et [[Présentation Fusion 2024.pdf]]
- une heuristique simple consiste à sélectionner uniquement les échantillons à l'intérieur d'une hypersphère
- puisque l'étape de krigeage intervient entre la prédiction et la correction, le centre de cette hypersphère est l'état prédit du filtre $\hat{x}_{k|k-1}$, _i.e._ le barycentre des particules à kriger
- le rayon de l'hypersphère va dépendre de l'ellipsoïde de confiance à $3\sigma$ car, pour une variable aléatoire $X\sim\mathcal{N}(\mu,\sigma)$, $$
\mathbb{P}(|X-\mu|\leq3\sigma)\approx99.7\%
$$ et donc la plupart des particules se trouvent en général dans cet ensemble
- considérer l'hypersphère circonscrite à l'ellipsoïde prédite, dont le rayon est $$
\max_{\lambda\in\text{sp}(\hat{P}_{k|k-1})}\left(3\sqrt{\lambda}\right)
$$ permet de sélectionner des échantillons uniformément dans l'espace, ce qui ne serait pas le cas si on sélectionnait selon la matrice de covariance d'une distribution très resserrée suivant certaines dimensions par rapport à d'autres
- on dilate ce rayon d'un facteur $\alpha$ afin de s'assurer que les particules les plus proches de l'extérieur soient interpolées et non extrapolées
- la sélection se fait donc comme illustré ci-dessous ![[Fenêtre de krigeage.png]]
- finalement, l'ensemble de sélection associé à cette méthode est $$
S_\text{fenêtrage}:=\mathcal{B}\left(\hat{x}_{k|k-1},\alpha\times\max_{\lambda\in\text{sp}(\hat{P}_{k|k-1})}\left(3\sqrt{\lambda}\right)\right),\quad\quad\alpha\geq1
$$

## Partitionnement des particules :

**TODO : partitionnement adaptatif et simultané**

- [[A Clustered Approach to Adaptive Kriging Particle Filter for Magnetic Field Aided Navigation.pdf]]
- si le fenêtrage est une heuristique simple et rapide à mettre en œuvre, elle présente 3 problèmes principaux :
	- il n'est pas certain que toutes les particules soient dans l'ellipsoïde à $3\sigma$, donc le krigeage va potentiellement extrapoler les observations de certaines particules
	- le rayon de la fenêtre dépend du nuage de particules, ce qui signifie que le nombre d'échantillons (et donc la durée de la régression) varie grandement avec la distribution des particules dans l'espace
	- en particulier, si le nuage est fortement multimodal, la fenêtre peut prendre en compte des échantillons loin de toutes les particules, ce qui n'apporte pas d'information utile
- le [[Partitionnement]] (_clustering_ en anglais) des particules permet de passer d'un nuage de forme quelconque à une union de sous-nuages assimilables à des ellipsoïdes $$
\left\{x_k^{(i)},i\in[\![1,N]\!]\right\}=\bigcup_{c_k\in\mathcal{C}_k}c_k
$$
- pour chaque _cluster_ $c_k\in\mathcal{C}_k$, on peut définir un état estimé $\hat{x}_{k|k-1}^{(c_k)}$ et une matrice de covariance $\hat{P}_{k|k-1}^{(c_k)}$ comme on l'aurait fait pour le nuage total
- l'ensemble associé à cette méthode est donc $$
S_\text{partitionnement}:=\bigcup_{c_k\in\mathcal{C}_k}\mathcal{B}\left(\hat{x}_{k|k-1}^{(c_k)},\alpha\times\max_{\displaystyle\small\lambda\in\text{sp}\left(\hat{P}_{k|k-1}^{(c_k)}\right)}\left(3\sqrt{\lambda}\right)\right),\quad\quad\alpha\geq1
$$

## Dilatation de Minkowski :

- dans le cas où toutes les particules sont associées à un même _cluster_, les ensembles considérés par les 2 méthodes précédentes sont confondus
- si à l'inverse, on considère que chaque particule appartient à son propre _cluster_, on ne peut plus définir les matrices $\hat{P}_{k|k-1}^{(c_k)}$, mais on peut tout de même généraliser la méthode
- pour $A$ et $B$ des ensembles, on définit la somme de Minkowski $$
A\oplus B:=\{a+b,a\in A,b\in B\}
$$
- en supposant que chaque _cluster_ $c_k$ de la méthode précédente a la même matrice $P$ comme matrice de covariance, on trouve alors $$
S_\text{partitionnement}=\left\{\hat{x}_{k|k-1}^{(c_k)},c_k\in\mathcal{C}_k\right\}\oplus\mathcal{B}\left(0,\alpha\times\max_{\displaystyle\small\lambda\in\text{sp}(P)}\left(3\sqrt{\lambda}\right)\right)
$$
- à la limite d'une particule par _cluster_, la méthode de dilatation de Minkowski propose donc d'introduire un rayon $r_k$ commun à toutes les particules, et de considérer l'ensemble de sélection $$
S_\text{dilatation}=\left\{x_k^{(i)},i\in[\![1,N]\!]\right\}\oplus\mathcal{B}\left(0,r_k\right)
$$
