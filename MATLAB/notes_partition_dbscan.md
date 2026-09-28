# DBSCAN contre mean shift pour l'ICAK et le MCAK — 2026-09-26

Note écrite à la demande de Bastien, pour l'annexe « sensibilité aux paramètres ». Tout ce qui
suit vient de lots à mille essais, appariés : mêmes relevés, mêmes états initiaux, mêmes bruits.
**Seule la partition change.**

## La question

Le mean shift impose une échelle par sa bande passante, et il fragmente : une vingtaine de modes
à l'acquisition, 5,3 en moyenne à $7{,}5$ km, dont les fenêtres se recouvrent — chaque point
retenu sert à près de quatre d'entre elles. DBSCAN lit les composantes connexes en densité : deux
amas reliés par une chaîne de particules restent un seul groupe, si éloignés soient leurs sommets.
L'intuition était qu'il supprimerait la redondance sans imposer de taille caractéristique.

## Ce qui a été fait

`krigeage.m`, fonction `partition_nuage` : le paramètre **`par.krig.partition`** vaut
`'meanshift'` (défaut, donc les lots antérieurs sont reproductibles) ou `'dbscan'`.

Trois choix d'implémentation, dont aucun n'est neutre :

- **$\varepsilon$ est proportionnel à la dispersion du nuage**, $0{,}3\,\sigma$ où $\sigma$ est la
  racine de la plus grande valeur propre de la covariance, et non une longueur en mètres. À huit
  kilomètres d'écart-type comme à cent cinquante, c'est le même nuage à l'échelle près, et un
  seuil absolu déclarerait tout connexe en poursuite. `dbscan_eps` impose une valeur fixe si on y
  tient, `dbscan_alpha` donne sinon le rapport.
- **`minpts` $=5$**, la règle usuelle $2d+1$ en dimension 2.
- **Le bruit n'est pas jeté** : un point étiqueté $-1$ rejoint le groupe le plus proche. Sinon on
  ne krigerait pas là où des particules se trouvent.

## Un paramètre non documenté qui conditionne tout

`mean_shift` **éclaircit le nuage à 500 points** avant de chercher les modes — un point sur dix,
pris à intervalle régulier — et `partition_dbscan` reprend la même règle. C'est dans le code
depuis le commit du 26 août, ça n'a jamais été discuté, et ça vaut pour **tous les lots du
chapitre**. Trois conséquences :

- $\varepsilon$ est calibré sur la densité de l'échantillon : $0{,}3\,\sigma$ vaut 2,7 distances
  entre voisins à cinq cents points, mais 8,6 à cinq mille. Le réglage n'est pas transposable.
- `minpts` $=5$ sur cinq cents points revient à cinquante sur cinq mille : un groupe doit porter
  environ 1 % du nuage pour exister.
- Un mode de moins de dix particules a une chance sur dix d'être échantillonné. Les modes faibles
  sont sous-comptés, systématiquement.

Correctif possible, non appliqué : $\varepsilon = c\,\sigma\sqrt{2\pi/n}$ avec $c\approx2{,}7$,
c'est-à-dire indexé sur la distance entre voisins. Le réglage devient alors insensible à
l'éclaircissage, qui redevient un pur paramètre de coût.

## Ce que ça change

**Le nombre de groupes s'effondre et se stabilise** : 2,0 à toutes les mailles, contre 2,5 à
$2{,}5$ km et 5,3 à $7{,}5$ km.

**La partition devient gratuite** : 3,36 Mflop contre 114 à 274 M, cinquante à quatre-vingts fois
moins, et la même valeur partout puisqu'elle ne dépend que du nombre de particules.

**Mais le krigeage explose**, en Gflop par essai :

| $\Delta_g$ | ICAK mean shift | ICAK DBSCAN | MCAK mean shift | MCAK DBSCAN |
|---|---|---|---|---|
| 2,5 km | 2,26 | **14,3** | 2,47 | **6,15** |
| 7,5 km | 0,50 | **0,73** | 0,15 | **0,25** |

La cause est l'acquisition. Tant que le nuage est large, DBSCAN n'y voit qu'un groupe et la
fenêtre monte au plafond de 2000 points, là où vingt modes donnaient vingt fenêtres d'une
centaine. Le coût étant quadratique en taille de fenêtre et cubique pour la factorisation, une
grande coûte bien plus que la somme de vingt petites, recouvrements compris. **La redondance
qu'on voulait supprimer était l'économie.**

## La précision : rien en médiane, beaucoup dans la queue

Convergence, cohérence et erreur finale sont indiscernables d'une partition à l'autre — 99,7 %
contre 99,7 % à $2{,}5$ km, 57,8 contre 58,3 % à $7{,}5$ km, Mahalanobis à 86–87 %, erreur finale
à $\pm20$ m. C'est attendu : le plafond est fixé par la carte et non par la partition, l'AK
atteint déjà la précision du krigeage complet, et cent vingt pas sur cent cinquante se passent
sur un nuage unimodal où la partition n'a pas d'objet.

L'écart est dans la **queue de l'acquisition**, à $7{,}5$ km, part des essais convergés encore
au-delà de cinq kilomètres :

| méthode | pas 20 | pas 30 | $q_{90}$ au pas 30 |
|---|---|---|---|
| AK | 4,2 % | 4,3 % | 2741 m |
| ICAK mean shift | **18,5 %** | 12,5 % | **6099 m** |
| ICAK DBSCAN | **10,5 %** | 7,9 % | **3776 m** |
| MCAK mean shift | 6,2 % | 3,3 % | 2717 m |
| MCAK DBSCAN | 4,6 % | 2,5 % | 2686 m |

**DBSCAN corrige près de la moitié du défaut de l'ICAK** : $-8$ points sur la queue, $-38\%$ sur
le neuvième décile. À $5{,}5$ km le même effet existe en petit (1,2 → 0,6 %) ; à $2{,}5$ km il n'y
a rien à corriger, aucun essai ne dépasse cinq kilomètres.

Ce qui reste — l'ICAK est encore deux fois et demie pire que l'AK — est ce que deux lois
distinctes suffisent à produire. Le MCAK, lui, n'a rien à gagner : une seule loi par construction,
quelle que soit la partition.

## Ce que ça dit du chapitre

1. **Le nombre de modes est un levier de coût**, pas une nuisance. Fragmenter rend les fenêtres
   petites, donc sert les troisième et quatrième critères du fenêtrage.
2. **Le défaut de l'ICAK est bien dans la multiplicité des lois**, et son ampleur suit le nombre de
   groupes : 5,3 modes donnent 18,5 % de queue, 2,0 en donnent 10,5, une seule loi en donne 4,6.
3. **Le MCAK est insensible au choix de l'algorithme de partition**, ce qui en fait l'argument le
   plus solide à écrire : il ne repose sur aucun réglage de clustering.

## Proposé, non fait

- **Balayer $\varepsilon\in\{0{,}15,\ 0{,}3,\ 0{,}5,\ 1\}\,\sigma$** : savoir si le résultat tient
  sur un plateau ou sur un point.
- **Balayer `minpts` $\in\{1,3,5\}$.** À 1, DBSCAN devient un pur union-find sur les composantes
  connexes : un paramètre au lieu de deux, mais plus aucun rempart contre l'effet de chaînage,
  qu'une seule particule posée entre deux modes suffit à déclencher.
- **Mesurer la partition sur le nuage complet**, une maille, pour savoir de combien
  l'éclaircissage à cinq cents points gonfle le compte des modes.
- **Relever la fraction de points noyaux et la distance d'un point de bruit à son centre
  d'accueil** : c'est elle qui dit si le rattachement élargit les fenêtres ou s'il est anodin.

## Reproduction

```
matlab -batch "cd('D:/bhubert/Desktop/Notes/MATLAB'); over.campagne='libre'; over.compare.models={'cak','mcak'}; over.krig.n_surveys=100; over.krig.partition='dbscan'; over.balayage=struct('champ','krig.pas','valeurs',[2500 3500 4500 5500 6500 7500]); over.mc.n_runs=1000; campagne_chapitre_4" > dbscan_1000.log 2>&1
```

Lots : `20260926_183545` (DBSCAN, ICAK et MCAK), `20260912_030842` (mean shift, ICAK), 
`20260925_095327` (mean shift, MCAK). Scripts d'analyse dans le scratchpad de la session
`ede9527c`, fichiers `cmp_partition.m` et `cmp_queue.m`.

**Attention au dépouillement** : `figures_campagnes` retient le lot le plus récent, donc
`20260926_183545` fournit désormais l'ICAK et le MCAK des figures du chapitre, en version DBSCAN.
Pour revenir au mean shift, il suffit de sortir ce lot de `resultats/` — les clés absentes sont
reprises dans les lots antérieurs, sans rien relancer.
