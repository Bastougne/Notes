## Principe :

- les hypothèse de linéarités sur les fonctions de modèle et de mesure sont parfois trop restrictives
- on peut linéariser au premier ordre ces fonctions autour du point d'estimation pour utiliser les équations de Kalman

## Hypothèses :

- $x_k=f_k(x_{k-1},u_k)+w_k$
- $z_k=h_k(x_k)+n_k$
- $f_k$ et $h_k$ sont connues et différentiables en les $x_k$
- $w_k$ et $n_k$ sont des bruits blancs gaussiens
- $x_0$ est gaussien

## Matrices jacobiennes :

- la matrice jacobienne d'une fonction $\varphi:\mathbb{R}^n\to\mathbb{R}^m$ différentiable en un point $x\in\mathbb{R}^n$ est donnée par : $$
\frac{\partial\varphi}{\partial x}=\begin{bmatrix}
\displaystyle{\frac{\partial\varphi_1}{\partial x_1}}&\dots&\dots&\displaystyle{\frac{\partial\varphi_1}{\partial x_n}}\\
\vdots&\ddots&\ddots&\vdots\\
\displaystyle{\frac{\partial\varphi_m}{\partial x_1}}&\dots&\dots&\displaystyle{\frac{\partial\varphi_m}{\partial x_n}}
\end{bmatrix}=\begin{bmatrix}
\displaystyle{\frac{\partial\varphi}{\partial x_1}}&\dots&\displaystyle{\frac{\partial \varphi}{\partial x_n}}
\end{bmatrix} 
$$
- on a le développement limité suivant $$
\varphi(x)=\varphi(x_0)+\frac{\partial\varphi}{\partial x}(x_0)\:\overrightarrow{xx_0}+o({||x-x_0||}_2)
$$
- en particulier, si $\varphi:x\mapsto F\:x+b$, alors $$
\forall x\in\mathbb{R}^n,\:\frac{\partial\varphi}{\partial x}(x)=F
$$
- on va donc remplacer $F_k$ et $H_k$ respectivement par $$
\frac{\partial f_k}{\partial x_{k-1}}(x_{k-1},u_k)\quad\text{et}\quad\frac{\partial h_k}{\partial x_k}(x_k)
$$

## Algorithme :

- posons $$
\tilde{F}_k=\frac{\partial f_k}{\partial x_{k-1}}(\hat{x}_{k-1|k-1},u_k)\quad\text{et}\quad\tilde{H}_k=\frac{\partial h_k}{\partial x_k}(\hat{x}_{k|k-1})
$$
- prédiction
	- $\hat{x}_{k|k-1}=f_k(\hat{x}_{k-1|k-1},u_k)$
	- $P_{k|k-1}=\tilde{F}_k\:P_{k-1|k-1}\:\tilde{F}_k^T+Q_k$
- correction
	- $y_k=z_k-h_k(\hat{x}_{k|k-1})$
	- $S_k=\tilde{H}_k\:P_{k|k-1}\:\tilde{H}_k^T+R_k$
	- $K_k=P_{k|k-1}\:\tilde{H}_k^T\:S_k^{-1}$
	- $\hat{x}_{k|k}=\hat{x}_{k|k-1}+K_k\:y_k$
	- $P_{k|k}=(I-K_k\:\tilde{H}_k)\:P_{k|k-1}$

## Remarque :
- l'EKF est une approximation analytique
- non-applicable si $f_k$ ou $h_k$ sont discontinues
- reste un filtre de Kalman donc suppose une loi _a posteriori_ gaussienne
- échoue si les non-linéarités sont telles que la loi réelle est fortement multimodale ou à queue lourde
- on peut rajouter des ordres supplémentaires dans la décomposition en série de Taylor (_higher-order EKF_, cf [[Smoothing Filtering and Prediction Estimating the Past Present and Future Nonlinear Prediction Filtering and Smoothing.pdf]]), mais pas beaucoup mieux pour un coût calculatoire important
- on peut aussi itérer l'étape de correction afin de linéariser au plus proche de l'état corrigé au lieu de l'état prédit (_iterative EKF_ ou _IEKF_, cf [[Performance Evaluation of Iterated Extended Kalman Filter with Variable Step-Length.pdf]]), mais fonctionne surtout quand on observe tout l'état (d'après Ristic)
