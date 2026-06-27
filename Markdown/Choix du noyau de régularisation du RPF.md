## Noyau d'Epanechnikov :

- introduit par [[Non-Parametric Estimation of a Multivariate Probability Density.pdf]] :
- forme du noyau : $$
K_{Epa}(x):=\frac{d+2}{2\:c_d}\left(1-||x||^2\right)\times\mathbb{1}_{\{||x||<1\}}(x)
$$ avec $c_d$ le volume de la boule unité de $\mathbb{R}^d$
- facteur de dilatation : $$
h_{Epa}:=A(K_{Epa})\:N^{-\frac{1}{d+4}}\quad\text{et}\quad
A(K_{Epa}):=\left(\frac{8(d+4)\left(2\sqrt{\pi}\right)^d}{c_d}\right)^{\frac{1}{d+4}}
$$
- minimise le MISE (_mean integrated squared error_ en anglais, erreur quadratique intégrée moyenne) si la loi a priori $p$ est gaussienne et que sa matrice de covariance est l'identité (on va donc appliquer une transformation linéaire, dite de _whitening_ pour que la matrice de covariance des particules soit $I_N$)
- en pratique, on multiplie $h_{Epa}$ par un coefficient $\mu\in[0,1;0,7]$ de réglage pour éviter le sur-lissage des multimodalités de la densité a posteriori
- $1/\mu$ correspond approximativement au nombre de modes de la densité a posteriori

## Noyau Gaussien

- forme du noyau : $$
K_{Gauss}(x):=(2\pi)^{-\frac{d}{2}}\:e^{-\frac{x^2}{2}}
$$
- facteur de dilatation : $$
h_{Epa}:=A(K_{Gauss})\:N^{-\frac{1}{d+4}}\quad\text{et}\quad
A(K_{Gauss}):=\left(\frac{4}{d+2}\right)^{\frac{1}{d+4}}
$$
- sous-optimal par rapport au noyau d'Epanechnikov, mais proche et plus simple à calculer
- permet aussi de combiner avec des corrections de Kalman
- de même, on multiplie $h_{Gauss}$ par $\mu\in[0,1;0,7]$ pour éviter le sur-lissage

## Noyau boîte :

- Nicolas a parlé d'un _Box Regularised Particle Filter_ dans sa thèse
- il utilise des noyaux dont le support est un produit cartésien d'intervalles
- en pratique, il utilise des noyaux uniformes sur ces boîtes pour simplifier les équations
