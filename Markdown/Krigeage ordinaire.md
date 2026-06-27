## Principe :

- généralisation du [[Krigeage simple]] au cas où la fonction de moyenne est inconnue mais constante
- on travaille toujours avec un processus Gaussien de $\mathbb{R}^d\rightarrow\mathbb{R}$
- comme précédemment, on cherche à trouver le meilleur estimateur linéaire non-biaisé
- on cherche donc  toujours un estimateur de la forme $$
z^*=\sum_{i=1}^{\tilde{n}}\lambda_i\tilde{z}_i
$$
- on souhaite également minimiser $$
\text{Cov}(z^*-z)=\Lambda^T\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)\Lambda-2\Lambda^Tk(\tilde{X},x)+k(x,x)+R
$$ sous la contrainte que $$
\begin{align}
&\mathbb{E}(z^*-z)=0\\
\iff\quad&\sum_{i=1}^{\tilde{n}}\lambda_i\mathbb{E}(\tilde{z}_i)-\mathbb{E}(z)=0\\
\iff\quad&\left(\sum_{i=1}^{\tilde{n}}\lambda_i-1\right)m=0\\
\iff\quad&\sum_{i=1}^{\tilde{n}}\lambda_i-1=0&\text{en toutes généralités}
\end{align}
$$

## Méthode des multiplicateurs de Lagrange :

- on pose donc le [Lagrangien](https://fr.wikipedia.org/wiki/Multiplicateur_de_Lagrange) $$
\mathcal{L}(\Lambda,\mu):=f(\Lambda)+2\langle\mu|g(\Lambda)\rangle=f(\Lambda)+2\mu g(\Lambda)
$$ avec $$
f(\Lambda):=\text{Cov}(z^*-z)\quad\text{et}\quad g(\Lambda):=\sum_{i=1}^{\tilde{n}}\lambda_i-1=\mathbf{1}^T\Lambda-1
$$
- cherchons un extremum de $\mathcal{L}(\Lambda,\mu)$ en dérivant par rapport à $\Lambda$ $$
\begin{align}
&\frac{d}{d\Lambda}\mathcal{L}(\Lambda,\mu)=0\\
\iff\quad&2\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)\Lambda-2k(\tilde{X},x)+2\mu\mathbf{1}=0\\
\iff\quad&\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)\Lambda+\mu\mathbf{1}=k(\tilde{X},x)
\end{align}
$$
- on a donc le système $$
\begin{align}
&\left\{\begin{matrix}\left(k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\right)\Lambda+\mu\mathbf{1}&=&k(\tilde{X},x)\\\mathbf{1}^T\Lambda&=&1\end{matrix}\right.\\
\iff\quad&\begin{bmatrix}k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}&\mathbf{1}\\\mathbf{1}^T&0\end{bmatrix}\begin{bmatrix}\Lambda\\\mu\end{bmatrix}=\begin{bmatrix}k(\tilde{X},x)\\1\end{bmatrix}\\
\iff\quad&\begin{bmatrix}\Lambda\\\mu\end{bmatrix}=\begin{bmatrix}k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}&\mathbf{1}\\\mathbf{1}^T&0\end{bmatrix}^{-1}\begin{bmatrix}k(\tilde{X},x)\\1\end{bmatrix}
\end{align}
$$
- sous cette forme, on voit une ressemblance forte avec l'équation du krigeage simple

## Coefficients du krigeage ordinaire :

- puisque la matrice à inverser est par blocs, on va utiliser l'inversion de Schur sur le bloc supérieur gauche (l'inférieur droit n'étant pas inversible)
- on a donc $$
\begin{bmatrix}\Lambda\\\mu\end{bmatrix}=\begin{bmatrix}\bar{A}&\bar{B}\\\bar{C}&\bar{D}\end{bmatrix}\begin{bmatrix}k(\tilde{X},x)\\1\end{bmatrix}
$$ avec $$
\left\{\begin{align}
G:=&k(\tilde{X},\tilde{X})+\tilde{R}I_{\tilde{n}}\\
S:=&0-\mathbf{1}^TG^{-1}\mathbf{1}&\in\mathbb{R}\\
\bar{A}=&G^{-1}+G^{-1}\mathbf{1}S^{-1}\mathbf{1}^TG^{-1}\\
\bar{B}=&-G^{-1}\mathbf{1}S^{-1}\\
\bar{C}=&-S^{-1}\mathbf{1}^TG^{-1}\\
\bar{D}=&S^{-1}
\end{align}\right.
$$
- les coefficients sont donc donnés par $$
\begin{align}
\Lambda=&\bar{A}\:k(\tilde{X},x)+\bar{B}\\
=&G^{-1}k(\tilde{X},x)+\frac{(G^{-1}\mathbf{1})(\mathbf{1}^TG^{-1})}{-\mathbf{1}^TG^{-1}\mathbf{1}}k(\tilde{X},x)-\frac{G^{-1}\mathbf{1}}{-\mathbf{1}^TG^{-1}\mathbf{1}}\\
=&G^{-1}k(\tilde{X},x)+G^{-1}\mathbf{1}\underbrace{\frac{1-\mathbf{1}^TG^{-1}k(\tilde{X},x)}{\mathbf{1}^TG^{-1}\mathbf{1}}}_{\in\mathbb{R}}
\end{align}
$$

## Estimateur du krigeage ordinaire :

- l'estimateur du krigeage ordinaire est donc $$
z^*_o(x,\tilde{X},\tilde{Z})=\Lambda^T\tilde{Z}=k(x,\tilde{X})G^{-1}\tilde{Z}+\frac{1-k(x,\tilde{X})G^{-1}\mathbf{1}}{\mathbf{1}^TG^{-1}\mathbf{1}}\mathbf{1}^TG^{-1}\tilde{Z}
$$
- on peut réécrire cet estimateur pour faire apparaître un terme de moyenne krigée $$
z^*_o(x,\tilde{X},\tilde{Z})=\frac{\mathbf{1}^TG^{-1}\tilde{Z}}{\mathbf{1}^TG^{-1}\mathbf{1}}+k(x,\tilde{X})G^{-1}\left(\tilde{Z}-\mathbf{1}\frac{\mathbf{1}^TG^{-1}\tilde{Z}}{\mathbf{1}^TG^{-1}\mathbf{1}}\right)
$$
- simplifions $\Lambda^TG\Lambda\in\mathbb{R}$ $$
\begin{align}
\Lambda^TG\Lambda=&\left(k(x,\tilde{X})G^{-1}+\frac{1-k(x,\tilde{X})G^{-1}\mathbf{1}}{\mathbf{1}^TG^{-1}\mathbf{1}}\mathbf{1}^TG^{-1}\right)G\Lambda\\
=&k(x,\tilde{X})G^{-1}G\Lambda+\frac{1-k(x,\tilde{X})G^{-1}\mathbf{1}}{\mathbf{1}^TG^{-1}\mathbf{1}}\mathbf{1}^TG^{-1}G\Lambda\\
\\
=&k(x,\tilde{X})\left(G^{-1}k(\tilde{X},x)+G^{-1}\mathbf{1}\frac{1-\mathbf{1}^TG^{-1}k(\tilde{X},x)}{\mathbf{1}^TG^{-1}\mathbf{1}}\right)+\\
&\frac{1-k(x,\tilde{X})G^{-1}\mathbf{1}}{\mathbf{1}^TG^{-1}\mathbf{1}}\mathbf{1}^T\left(G^{-1}k(\tilde{X},x)+G^{-1}\mathbf{1}\frac{1-\mathbf{1}^TG^{-1}k(\tilde{X},x)}{\mathbf{1}^TG^{-1}\mathbf{1}}\right)\\
\\
=&k(x,\tilde{X})G^{-1}k(\tilde{X},x)+\frac{1-k(x,\tilde{X})G^{-1}\mathbf{1}}{\mathbf{1}^TG^{-1}\mathbf{1}}\mathbf{1}^TG^{-1}\mathbf{1}\frac{1-\mathbf{1}^TG^{-1}k(\tilde{X},x)}{\mathbf{1}^TG^{-1}\mathbf{1}}+\\
&\underbrace{k(x,\tilde{X})G^{-1}\mathbf{1}}_{\in\mathbb{R}}\frac{1-\mathbf{1}^TG^{-1}k(\tilde{X},x)}{\mathbf{1}^TG^{-1}\mathbf{1}}+\frac{1-k(x,\tilde{X})G^{-1}\mathbf{1}}{\mathbf{1}^TG^{-1}\mathbf{1}}\underbrace{\mathbf{1}^TG^{-1}k(\tilde{X},x)}_{\in\mathbb{R}}\\
\\
=&k(x,\tilde{X})G^{-1}k(\tilde{X},x)+\frac{\left(1-k(x,\tilde{X})G^{-1}\mathbf{1}\right)^2}{\mathbf{1}^TG^{-1}\mathbf{1}}+2\frac{1-k(x,\tilde{X})G^{-1}\mathbf{1}}{\mathbf{1}^TG^{-1}\mathbf{1}}\mathbf{1}^TG^{-1}k(\tilde{X},x)\\
\end{align}
$$
- la matrice de covariance du krigeage ordinaire est la variance minimisée $$
\begin{align}
R^*_o(x,\tilde{X})=&\Lambda^TG\Lambda-2\Lambda^Tk(\tilde{X},x)+k(x,x)+R\\
=&k(x,\tilde{X})G^{-1}k(\tilde{X},x)+\frac{\left(1-k(x,\tilde{X})G^{-1}\mathbf{1}\right)^2}{\mathbf{1}^TG^{-1}\mathbf{1}}+\\
&2\frac{1-k(x,\tilde{X})G^{-1}\mathbf{1}}{\mathbf{1}^TG^{-1}\mathbf{1}}\mathbf{1}^TG^{-1}k(\tilde{X},x)-2k(x,\tilde{X})G^{-1}k(\tilde{X},x)-\\
&2\frac{1-k(x,\tilde{X})G^{-1}\mathbf{1}}{\mathbf{1}^TG^{-1}\mathbf{1}}\mathbf{1}^TG^{-1}k(\tilde{X},x)+k(x,x)+R\\
\\
=&k(x,x)+R-k(x,\tilde{X})G^{-1}k(\tilde{X},x)+\frac{\left(1-k(x,\tilde{X})G^{-1}\mathbf{1}\right)^2}{\mathbf{1}^TG^{-1}\mathbf{1}}
\end{align}
$$

## Lien avec le krigeage simple :

- en posant $$
m_o(\tilde{X},\tilde{Z}):=\frac{\mathbf{1}^TG^{-1}\tilde{Z}}{\mathbf{1}^TG^{-1}\mathbf{1}}
$$ on remarque que $$
z^*_o(x,\tilde{X},\tilde{Z})=z^*_s(x,\tilde{X},\tilde{Z},m_o(\tilde{X},\tilde{Z}))
$$
- de plus, en posant $$
R^*_{m_o}(x,\tilde{X}):=\frac{\left(1-k(x,\tilde{X})G^{-1}\mathbf{1}\right)^2}{\mathbf{1}^TG^{-1}\mathbf{1}}
$$ il vient que $$
R^*_o(x,\tilde{X})=R^*_s(x,\tilde{X})+R^*_{m_o}(x,\tilde{X})
$$

## Résultats :

- Krigeage 2D : ![[GPR 2d MATLAB.png]]
- Krigeage 3D : ![[GPR 3d MATLAB.png]]
- Krigeage sur carte : ![[GPR carte MATLAB.png]]
- RPF avec krigeage (construction de carte par données bruitées sans correction en cours de mission) : ![[RPF avec krigeage.png]]

## Remarques :

- [[Choix du noyau de covariance du krigeage]]
- [[Estimation des hyper-paramètres du krigeage]]
- la [[Sélection des données pour le krigeage]] et la [[Réduction de rang pour le krigeage]] visent à réduire le temps de calcul nécessaire au krigeage sans trop dégrader sa qualité
- le krigeage peut être combiné avec un filtre particulaire dans un paradigme SLAM comme le montre [[Terrain-Aided SLAM with Limited-Size Reference Maps using Gaussian Processes.pdf]]
- [[Lien entre le krigeage et les espaces de Hilbert reproduisants]]
