## But :

- [[Resampling in Particle Filtering Comparison.pdf]]
- [[Resampling Methods for Particle Filtering Classification Implementation and Strategies.pdf]]
- basé sur la fonction de répartition des poids (_cumulative distribution function_ ou CDF en anglais) :
![[CDF resampling.png]]
- certains algorithmes de filtrage particulaire se servent de l'indice de la particule parent de la nouvelle particule $i^{(j)}$, on la calcule donc ici directement

## Rééchantillonnage multinomial

-  $$\left(\left\{x_k^{(j)*},w_k^{(j)*},i^{(j)}\right\}_{j=1}^N\right)=Resampling\left(\left\{x_k^{(i)},w_k^{(i)}\right\}_{i=1}^N\right)$$
	- calcul de la CDF : $$\left\{\begin{align}cdf_1=&0&\\cdf_i=&cdf_{i-1}+w_k^{(i)},&i\in[\![2,N]\!]\end{align}\right.$$
	- tirage : $$\forall j\in[\![1,N]\!],u_j\sim\mathbb{U}([0,1])$$
	- POUR $j=1:N$
		- recherche (dichotomique car croissance de $(cdf_i)_{i=1}^N$) du parent : $$i\text{ tel que }cdf_i>u_j>cdf_{i-1}$$
		- création des nouvelles particules : $$\left\{\begin{align}x_k^{(j)*}=&x_k^{(i)}\\w_k^{(j)*}=&\frac{1}{N}\\i^{(j)}=&i\end{align}\right.$$
	- FIN POUR

## Rééchantillonnage de Kitagawa :

-  $$\left(\left\{x_k^{(j)*},w_k^{(j)*},i^{(j)}\right\}_{j=1}^N\right)=Resampling\left(\left\{x_k^{(i)},w_k^{(i)}\right\}_{i=1}^N\right)$$
	- calcul de la CDF : $$\left\{\begin{align}cdf_1=&0&\\cdf_i=&cdf_{i-1}+w_k^{(i)},&i\in[\![2,N]\!]\end{align}\right.$$
	- tirage : $$\left\{\begin{align}u_1\sim&\mathbb{U}([0,1/N])&\\u_j=&u_1+\frac{j-m}{N},&j\in[\![2,N]\!]&\quad\text{avec }m\text{ fixé dans }[0,1]\end{align}\right.$$
	- $i=1$
	- POUR $j=1:N$
		- TANT QUE $u_j>cdf_i$
			- $i=i+1$
		- FIN TANT QUE
		- création des nouvelles particules : $$\left\{\begin{align}x_k^{(j)*}=&x_k^{(i)}\\w_k^{(j)*}=&\frac{1}{N}\\i^{(j)}=&i\end{align}\right.$$
	- FIN POUR

## Rééchantillonnage déterministe :

-  $$\left(\left\{x_k^{(j)*},w_k^{(j)*},i^{(j)}\right\}_{j=1}^N\right)=Resampling\left(\left\{x_k^{(i)},w_k^{(i)}\right\}_{i=1}^N\right)$$
	- calcul de la CDF : $$\left\{\begin{align}cdf_1=&0&\\cdf_i=&cdf_{i-1}+w_k^{(i)},&i\in[\![2,N]\!]\end{align}\right.$$
	- tirage : $$\forall j\in[\![1,N]\!],u_j=\frac{j-m}{N}\quad\text{avec }m\text{ fixé dans }[0,1]$$
	- $i=1$
	- POUR $j=1:N$
		- TANT QUE $cdf_i<u_j$
			- $i=i+1$
		- FIN TANT QUE
		- création des nouvelles particules : $$\left\{\begin{align}x_k^{(j)*}=&x_k^{(i)}\\w_k^{(j)*}=&\frac{1}{N}\\i^{(j)}=&i\end{align}\right.$$
	- FIN POUR
