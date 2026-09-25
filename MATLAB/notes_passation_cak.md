# CA-OK réuni ou distinct — ce qui a été fait dans la conversation du 2026-09-22

Note de passation, écrite à la demande de Bastien, qui a mené ce sujet dans deux conversations en parallèle. Elle ne dit que ce qui vient de celle-ci. Le détail chiffré est dans `notes_localite_multimodalite.md`.

## Vocabulaire

- **CA-OK distincts** : la version du manuscrit. Une fenêtre, un ajustement d'hyperparamètres et une résolution par mode du mean-shift.
- **CA-OK réunis** : les modes ne servent qu'à *choisir* les points. On réunit leurs fenêtres, on ajuste une fois et on krige une fois, donc toutes les particules sont notées sous la même loi.

## Dans le code

`krigeage.m`, `query_cak` : option **`par.krig.cak_union`**, ajoutée le 2026-09-22. À vrai, elle rend la version réunie ; absente ou fausse, rien ne change, et les lots antérieurs restent reproductibles.

## Les lots, tous à 100 essais sur 10 relevés

Ce sont exactement les cent premiers essais de `20260912_030842`, donc tout se compare essai par essai avec l'AK, le CA-OK et la NNAK du manuscrit.

| lot | réglages | modèles |
|---|---|---|
| `20260917_111800` | $\alpha_{\mathrm{cov}} = 1{,}2$, $\tilde{N}_{\min} = 10$, bande 1 km | AK, CA-OK |
| `20260917_113709` | idem, bande 4 km | CA-OK |
| `20260922_120645` | réglages du manuscrit, `refit='global'` | AK, CA-OK |
| `20260922_124202` | réglages du manuscrit, `cak_union=true` | CA-OK |

## Ce qu'ils disent

1. **Le CA-OK perd en RMSE à l'acquisition par sa queue de distribution, pas par sa médiane.** À 7,5 km, au pas 20, 30 % de ses essais convergés sont encore au-delà de 5 km, contre 10 % pour l'AK. Le premier virage les remet tous en place, et les erreurs finales sont identiques.
2. **La cause est l'ajustement des hyperparamètres par cluster.** Avec un seul jeu partagé, l'écart du CA-OK à l'AK tombe de +215 m à +7 m au pas 20 à 7,5 km, et disparaît à toutes les mailles.
3. **Le CA-OK réuni fait ce qu'on attend de lui** — plus d'écart à l'AK nulle part, queue ramenée de 30 à 15 % — **mais il rend l'avantage de coût** : moitié de l'AK aux mailles fines, autant que lui aux mailles lâches, là où l'union a la taille de la fenêtre de l'AK.
4. **La NNAK est le CA-OK réuni poussé à une particule par cluster**, avec un plancher de 4 au lieu de 25. Les quatre méthodes forment donc une famille à deux réglages : la finesse de la partition, et le fait de partager ou non un seul modèle de krigeage.
5. **Résultat de côté** : aux mailles lâches, les hyperparamètres globaux battent les locaux. À 7,5 km, la convergence passe de 50 à 72 % et le Mahalanobis de 88 à 80 % : le réajustement local achète de la cohérence au prix de la précision.
6. **Les fenêtres par mode se partagent leurs points** (mesuré le 2026-09-16) : une vingtaine de modes à l'acquisition, chaque point retenu servant à trois à sept fenêtres. Baisser $\alpha_{\mathrm{cov}}$ ou $\tilde{N}_{\min}$ n'y change rien ; élargir la bande du mean-shift réduit le recouvrement sans le supprimer.

## Proposé, non fait

- **Hyperparamètres partagés mais ajustés par pas sur l'union, résolutions toujours par cluster** : le coût du CA-OK avec la précision de l'AK. Une quinzaine de lignes dans `query_cak`.
- **La NNAK à 25 voisins par particule**, c'est-à-dire la même famille au plancher du chapitre : si sa cohérence remonte, son effondrement à 7,5 km venait bien du plancher.
- **Fusionner les clusters dont les fenêtres se touchent** : proposé le 2026-09-16, **rejeté par Bastien**, la covariance d'un cluster fusionné faisant enfler la fenêtre jusqu'à revenir à l'AK.

## Rédaction du 2026-09-22

Le code est mené dans l'autre conversation, la rédaction ici. Écrit et compilé, 172 pages sans référence non résolue ni débordement :
- **§4.4.2 Clustering**, en une seule section. Bastien voulait la structure du filtre optimal, puis de Kalman, puis des filtres particulaires, **mais l'analogie reste hors du texte** : elle ordonne la section, elle ne s'y dit pas. Contenu — HCA avec la récurrence de Lance-Williams, l'incrément d'inertie de Ward et les deux coupes du dendrogramme, puis $k$-means avec les deux minimisations partielles et l'algorithme de Lloyd, puis le mean shift avec la densité à noyau, la factorisation de son gradient et l'itération à point fixe, puis DBSCAN avec ses particules cœur et sa relation de connexité. La section se termine sur la partition retenue, chaque particule étant rattachée au mode le plus proche.
- **§4.4.3 Clustered Adaptive Kriging** — moments et fenêtre par cluster, puis les deux variantes, l'ICAK qui garde les fenêtres séparées et le MCAK qui les réunit, la différence de loi entre les deux, et les limites AK et NNAK selon la bande passante.
- **Les deux encadrés d'algorithme**, qui étaient des copies de celui de l'AK, réécrits.
- **Cinq références ajoutées** à vérifier : MacQueen 1967, Lloyd 1982, Ward 1963, Lance-Williams 1967, Ester et al. 1996. **HCA et DBSCAN** ajoutés à la liste des abréviations.

Restent à faire : la référence des travaux d'Achille sur les filtres clusterisés, les légendes de §4.4.4, et une éventuelle figure de dendrogramme.

## Où est le reste

`notes_localite_multimodalite.md` porte les tableaux complets, les coûts, les quantiles et le plan d'expériences. Les scripts d'analyse et les pilotes sont dans le scratchpad de la session `408c2db4-…`.
