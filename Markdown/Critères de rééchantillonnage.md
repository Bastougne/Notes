## Taille effective :

- la "taille effective" $N_{eff}$ de l'échantillon est défini par : $$
N_{eff} := \frac{N}{1 +
\text{Var}\left(\frac{p\left(\textbf{x}_{0:k}^{(i)}|\textbf{z}_{1:k}\right)}
{q\left(\textbf{x}_{0:k}^{(i)}|\textbf{z}_{1:k}\right)}\right)}
$$
- on ne peut pas évaluer $N_{eff}$, mais on peut en calculer un estimateur $\widehat{N_{eff}}$ : $$
\widehat{N_{eff}} := \frac{1}
{\displaystyle{\sum_{i=1}^N\left(w_k^{(i)}\right)}} \leq N
$$
- $\widehat{N_{eff}}$ faible indique une forte dégénérescence
- on ne rééchantillonne que quand $\widehat{N_{eff}} < N_{seuil}$

## Entropie

**TODO ?**
