# Localité, multimodalité et ce que le CA-OK gagne vraiment

Note d'analyse, écrite le 2026-09-15. Elle part d'une question de Bastien — *la réunion marche aussi bien que le clustering ; si le champ était vraiment non stationnaire, le CA-OK, qui ajuste ses hyperparamètres par mode, n'aurait-il pas dû gagner ?* — et y répond avec les deux lots à 1000 essais qui portent les trois méthodes. Elle sert trois choses : garder les chiffres, **justifier plusieurs choix du chapitre**, et poser les expériences qui rendraient multimodalité et instationnarité visibles.

Les chiffres sont mesurés et citables. Les **hypothèses** sont signalées comme telles : aucune n'est testée.

## Les lots et les métriques

- `resultats/20260912_030842_reference_libre_1000runs.mat` — les quatre modèles contre la maille, scénario `reference`, $\sigma_0 = 8$ km, six mailles de 2,5 à 7,5 km.
- `resultats/20260914_235935_reference_libre_1000runs.mat` — le régime multimodal, $\sigma_0 = 12$ km, cinq mailles de 4,5 à 7,5 km, avec le statique et la carte connue.

RPF, 5000 particules, 150 pas, $\alpha_{\mathrm{cov}} = 2$, $\tilde{N}_{\min} = 25$, $\tilde{N}_{\max} = 2000$, bande passante du mean-shift 1000 m. **acq** : erreur médiane au pas 30. **fin** : médiane sur les essais de la médiane des dix derniers pas. **Mahal.** : taux de convergence au critère de Mahalanobis (99 %). **multi** : part des couples (pas, essai) où le mean-shift du CA-OK rend plus d'un mode. **pts** : nombre moyen de points retenus. Script de dépouillement : `scratchpad/analyse_cak_reunion.m` de la session du 2026-09-15 (temporaire ; il tient en quarante lignes et se réécrit).

## Les chiffres

### Scénario de référence, $\sigma_0 = 8$ km

| maille | AK acq / fin / Mahal. | CA-OK acq / fin / Mahal. (multi) | réunion acq / fin / Mahal. (pts) |
|---|---|---|---|
| 2,5 km | 364 / 239 / 100 % | 365 / 242 / 100 % (7 %) | 386 / 307 / 100 % (22) |
| 3,5 km | 473 / 381 / 100 % | 475 / 379 / 100 % (10 %) | 590 / 488 / 99 % (18) |
| 4,5 km | 598 / 663 / 99 % | 608 / 663 / 99 % (16 %) | 748 / 757 / 97 % (15) |
| 5,5 km | 772 / 1035 / 90 % | 812 / 1052 / 91 % (26 %) | 828 / 1086 / 86 % (14) |
| 6,5 km | 861 / 1357 / 87 % | 942 / 1388 / 87 % (39 %) | 926 / 1330 / 81 % (14) |
| 7,5 km | 1459 / 1684 / 87 % | 1699 / 1720 / 87 % (54 %) | 1450 / 1813 / 74 % (13) |

Au critère de seuil (2000 m), les trois restent à trois points les uns des autres à toutes les mailles.

### Régime multimodal, $\sigma_0 = 12$ km

| maille | statique acq / fin / Mahal. | AK | CA-OK (multi) | réunion (pts) |
|---|---|---|---|---|
| 4,5 km | 682 / 773 / 99 % | 627 / 667 / 98 % | 640 / 664 / 98 % (17 %) | 828 / 746 / 97 % (22) |
| 5,0 km | 692 / 910 / 95 % | 643 / 802 / 94 % | 665 / 804 / 94 % (23 %) | 784 / 803 / 92 % (21) |
| 5,5 km | 816 / 1142 / 89 % | 881 / 1042 / 90 % | 906 / 1050 / 90 % (28 %) | 944 / 1087 / 86 % (20) |
| 6,5 km | 851 / 1474 / 83 % | 1003 / 1375 / 86 % | 1162 / 1397 / 86 % (42 %) | 1055 / 1416 / 79 % (19) |
| 7,5 km | 1297 / 1704 / 78 % | 1673 / 1726 / 86 % | 2129 / 1790 / 84 % (57 %) | 1679 / 1808 / 75 % (18) |

## Ce que disent les chiffres

1. **Le CA-OK n'est meilleur que l'AK nulle part.** Aux mailles lâches il est **pire à l'acquisition** — +9 et +16 % à 6,5 et 7,5 km en référence, +16 et +27 % en multimodal — et rejoint l'AK en fin de vol. C'est précisément là qu'il est multimodal le plus souvent (39 à 57 % des pas).
2. **La réunion n'est pas aussi bonne.** Aux mailles fines, 14 à 28 % d'erreur en plus (307 contre 239 m à 2,5 km, 488 contre 381 à 3,5 km). Aux mailles lâches, elle retrouve la précision de l'AK mais perd jusqu'à **13 points de Mahalanobis** (74 contre 87 % à 7,5 km).
3. **L'AK fait déjà ce qu'on reproche à la réunion, sans en souffrir.** Pendant l'acquisition sa fenêtre unique, de rayon $6\sqrt{\lambda_{\max}}$ — 48 km au départ avec $\sigma_0 = 8$ km —, couvre tous les modes avec **un seul ajustement et une seule moyenne**. Il ne perd rien contre le CA-OK : les hyperparamètres varient donc peu à l'échelle qui sépare les modes, de l'ordre de $\sigma_0$.
4. **Le champ est non stationnaire, mais à grande échelle.** En multimodal, l'AK bat le statique de **7 à 14 % en erreur finale** de 4,5 à 6,5 km, sa NEES finale tombe de 0,67 à 0,42 à 4,5 km, et il gagne 8 points de Mahalanobis à 7,5 km. Les hyperparamètres locaux servent ; leur localité *par mode* n'ajoute rien.
5. **Mais le statique acquiert mieux que l'AK dès 5,5 km en multimodal** (816 contre 881 m, 851 contre 1003, 1297 contre 1673), avant de finir moins bien. Pendant l'acquisition, le fenêtrage coûte quelque chose que la fin de vol rembourse.

## Pourquoi, d'après le code

- **Sur un pas unimodal, le CA-OK est l'AK** : `query_cak` appelle `query_ak` quand le mean-shift rend un mode. C'est 90 % des pas aux mailles fines ; le CA-OK ne peut s'écarter de l'AK que sur les pas multimodaux.
- **La réunion** (`query_boules`) prend les quatre plus proches échantillons de chaque particule, sans plancher, et fait **un seul** ajustement et **une seule** résolution sur leur union — 13 à 22 points en moyenne, sous $\tilde{N}_{\min}$.
- **Là où le CA-OK est multimodal, il ne peut pas être local.** Une fenêtre au plancher de 25 points couvre un disque de rayon $\approx 2{,}8\,\Delta$ : 7 km à 2,5 km de maille, **21 km à 7,5 km**, soit quatre longueurs de corrélation ($\ell \approx 5$ km). La séparation $\ell \ll r \ll L_\theta$ de la note de stationnarité est intenable aux mailles lâches — qui sont justement celles où le CA-OK se partitionne.

**Mesuré le 2026-09-16**, sur le lot de référence : rayon équivalent de la fenêtre de l'AK, $\Delta\sqrt{\tilde{N}^r/\pi}$, calculé essai par essai puis moyenné. Il vaut **48,0 km au premier pas à toutes les mailles**, soit exactement $6\sigma_0$, ce qui valide l'équivalence. Il redescend ensuite jusqu'au plancher, atteint vers le pas 15 à 2,5 km mais seulement vers le pas 50 à 7,5 km : **7,1 / 9,9 / 12,7 / 15,5 / 18,3 / 21,2 km** de 2,5 à 7,5 km, soit $2{,}8\,\Delta$, ou de 1,4 à 4,2 longueurs de corrélation. La fenêtre est au plancher sur 91 à 95 % des couples (pas, essai) : **en poursuite, c'est la maille, et non l'incertitude du filtre, qui fixe la localité.**

## Hypothèses, non testées

- **H1 — la réunion perd la cohérence par son nombre de points** (critère 2), et non parce qu'elle enjambe les modes (critère 3). Point 3 à l'appui : l'AK enjambe les modes lui aussi, mais avec au moins 25 points, et ne perd rien.
- **H2 — le CA-OK perd à l'acquisition parce que le mean-shift fragmente des nuages larges mais unimodaux.** 54 % de pas « multimodaux » à 7,5 km dès $\sigma_0 = 8$ km, c'est beaucoup pour un vrai multimodal. Des morceaux d'un même mode sont alors krigés sous des ajustements différents, sur des fenêtres plus petites.
- **H3 — réajuster les hyperparamètres à chaque pas injecte un bruit d'estimation qui coûte pendant l'acquisition**, et la localité ne paie qu'une fois le nuage resserré. Point 5 à l'appui, et l'ordre statique < AK < CA-OK des erreurs d'acquisition aux mailles lâches suit l'ordre des effectifs d'ajustement.
- **H4 — les hyperparamètres varient à une échelle $L_\theta$ bien plus grande que l'écart entre modes**, mais plus petite que la carte.

## Ce que ça justifie dans le chapitre

- **$\tilde{N}_{\min}$** : si H1 tient, la réunion est la démonstration empirique du critère 2 — elle perd jusqu'à 28 % d'erreur et 13 points de cohérence pour être passée sous le plancher. C'était déjà son rôle déclaré dans `notes_campagne_E.md` (« le second échec motive $n_{\min}$ »).
- **Le réajustement local** (critère 3) : l'AK bat le statique en fin de vol en multimodal. La localité se justifie à l'échelle de dizaines de kilomètres, pas à celle d'un mode.
- **Le CA-OK se justifie par le coût et par la taille de fenêtre pendant l'acquisition** (critère 4), pas par la précision ni par la localité des hyperparamètres — cohérent avec le choix de Bastien d'introduire le CAK par le plafond de fenêtre. **Présenter la réunion comme la preuve qu'il faut clusteriser serait fragile.**
- **La double inégalité $\ell \ll r \ll L_\theta$** reçoit une conséquence mesurable : aux mailles lâches, la fenêtre minimale est déjà trop grande pour être locale.
- **Le critère de Mahalanobis** est à nouveau celui qui départage les méthodes ; le seuil ne voit presque rien.

## Balayage de $\tilde{N}_{\min}$ à 3,5 km (dépouillé le 2026-09-16)

Lot `resultats/20260908_063919_reference_libre_1000runs.mat` : **CA-OK**, maille 3,5 km, $\sigma_0 = 8$ km, $\alpha_{\mathrm{cov}} = 2$, 1000 essais, valeurs 15, 25, 40 et 60. Entre parenthèses, l'écart apparié à $\tilde{N}_{\min} = 25$ — mêmes graines, moyenne ± erreur-type, en mètres.

| $\tilde{N}_{\min}$ | points retenus | coût (Gflop) | erreur au pas 5 | au pas 30 | fin | Mahal. |
|---|---|---|---|---|---|---|
| 15 | 37,7 | 0,76 | 2734 m (−175 ± 46) | 557 m (+46 ± 20) | 397 m (+11 ± 11) | 99,7 % |
| 25 | 47,3 | 1,14 | 2905 m | 475 m | 379 m | 99,6 % |
| 40 | 61,8 | 2,11 | 3178 m (+130 ± 42) | 451 m (−107 ± 45) | 358 m (−46 ± 22) | 99,8 % |
| 60 | 81,1 | 4,31 | 3205 m (+311 ± 48) | 478 m (−16 ± 26) | 387 m (+69 ± 37) | 99,7 % |

- **Les taux de convergence et la NEES finale (0,19 à 0,23) ne bougent pas** : à 3,5 km, le plancher ne décide ni de la convergence ni de la cohérence.
- **Sous 25**, l'acquisition se dégrade (+46 m au pas 30, z = 2,3) sans effet en fin de vol ; le coût baisse d'un tiers.
- **À 40**, un gain modeste à partir du pas 30 (−46 m en fin de vol, z de −2,2 à −2,6) pour un coût multiplié par 1,85 ; **à 60**, plus aucun gain, pour un coût multiplié par 3,8. Avec 24 comparaisons, des z de 2,5 sont des indices, pas des preuves.
- **Au pas 5, un plancher haut coûte cher** (+311 m à 60, z = 6,5) alors que la fenêtre de l'AK y compte des centaines de points : pendant l'acquisition le CA-OK partitionne, et c'est le plancher *par mode* qui s'applique. Probablement propre au CA-OK — à vérifier sur l'AK.
- **Rayon de la fenêtre en poursuite**, $\approx\Delta\sqrt{\tilde{N}_{\min}/\pi}$ : 7,7 / 9,9 / 12,5 / 15,3 km, soit 1,5 à 3 longueurs de corrélation. Le gain s'arrête vers 2,5 ℓ, ce qui cadre avec l'effet d'écran — une lecture, pas une démonstration.
- **Manque** : le même balayage à 5,5 km (dans la file, arrêtée) et sur l'AK plutôt que le CA-OK.
- **En poursuite, ce balayage est celui de l'AK** : à 3,5 km le CA-OK n'est multimodal que sur 2,7 % des pas de poursuite, et sur un pas unimodal il *est* l'AK. Seule la partie acquisition (le pas 5) est propre au CA-OK.

## Taille et recouvrement des fenêtres du CA-OK (mesuré le 2026-09-16)

Question de Bastien : quelle est la taille des clusters, et se partagent-ils les mêmes points ? S'ils le font, la comparaison AK / CA-OK / réunion ne dit plus rien sur la localité.

**Mesure**, sans relancer : sur un pas multimodal, `cout(1)` est la somme $S$ des tailles des fenêtres (un point partagé compte autant de fois qu'il sert), `cout(3)` le nombre de fenêtres et `n_used` leur union $U$. `cout` est moyenné sur les essais, mais sur un pas unimodal `cout(1) = n_used` et `cout(3) = 1` : on retranche les pas unimodaux essai par essai, et il reste les sommes exactes sur les pas multimodaux. $S/U$ est le nombre moyen de fenêtres auxquelles appartient un point retenu : 1 si les fenêtres sont disjointes, le nombre de fenêtres si elles sont identiques. Rayon équivalent d'une fenêtre : $\Delta\sqrt{n/\pi}$. Chaque mode reçoit une fenêtre (colonnes modes et fenêtres identiques en RPF). Acquisition = pas 1 à 30, poursuite = la suite.

### Lot de référence, $\sigma_0 = 8$ km

| maille | acq. : pas multi | modes | points / fenêtre (rayon) | union | $S/U$ | pours. : pas multi | modes | points / fenêtre | union | $S/U$ |
|---|---|---|---|---|---|---|---|---|---|---|
| 2,5 km | 34 % | 23,0 | 81 (12,7 km) | 568 | 3,3 | 0,4 % | 2,02 | 25,0 | 28,9 | 1,74 |
| 3,5 km | 37 % | 24,5 | 49 (13,8 km) | 323 | 3,7 | 2,7 % | 2,04 | 25,0 | 28,2 | 1,81 |
| 4,5 km | 44 % | 23,5 | 37 (15,4 km) | 196 | 4,4 | 9,0 % | 2,12 | 25,1 | 28,4 | 1,88 |
| 5,5 km | 56 % | 21,5 | 31 (17,4 km) | 130 | 5,2 | 18,7 % | 2,32 | 25,1 | 28,5 | 2,04 |
| 6,5 km | 72 % | 20,5 | 28 (19,5 km) | 95 | 6,2 | 30,1 % | 2,64 | 25,2 | 29,1 | 2,28 |
| 7,5 km | 85 % | 21,1 | 27 (22,0 km) | 77 | 7,4 | 45,8 % | 3,35 | 25,3 | 30,3 | 2,80 |

### Lot multimodal, $\sigma_0 = 12$ km

| maille | acq. : pas multi | modes | points / fenêtre (rayon) | union | $S/U$ | pours. : pas multi | modes | points / fenêtre | union | $S/U$ |
|---|---|---|---|---|---|---|---|---|---|---|
| 4,5 km | 50 % | 26,7 | 53 (18,5 km) | 386 | 3,7 | 8,9 % | 2,15 | 25,4 | 29,1 | 1,88 |
| 5 km | 57 % | 26,2 | 46 (19,1 km) | 309 | 3,9 | 14,7 % | 2,31 | 25,3 | 29,2 | 2,00 |
| 5,5 km | 65 % | 24,8 | 41 (19,9 km) | 247 | 4,1 | 19,2 % | 2,51 | 25,9 | 31,3 | 2,07 |
| 6,5 km | 81 % | 24,9 | 34 (21,5 km) | 173 | 4,9 | 32,3 % | 3,06 | 26,2 | 33,3 | 2,40 |
| 7,5 km | 90 % | 26,9 | 31 (23,4 km) | 138 | 6,0 | 48,7 % | 4,03 | 25,9 | 33,9 | 3,08 |

**Mêmes (pas, essai), tous les pas multimodaux du CA-OK confondus** — les graines sont communes, les trajectoires divergent ensuite, c'est donc une comparaison de situations et non de pas identiques : union du CA-OK / fenêtre de l'AK / union de la réunion = 543 / 680 / 167 points à 2,5 km, 258 / 316 / 92 à 3,5 km, 45 / 44 / 17 à 7,5 km (référence) ; 238 / 292 / 81 à 4,5 km, 67 / 73 / 26 à 7,5 km (multimodal).

**Balayage de $\tilde{N}_{\min}$** (3,5 km) : à l'acquisition, 24 modes quelle que soit la valeur, et une union qui ne bouge pas (320, 323, 321, 333 points de 15 à 60) pendant que les fenêtres grossissent (43,4 / 48,9 / 59,4 / 75,5 points, $S/U$ de 3,3 à 5,5). **Relever le plancher n'élargit pas la zone couverte, il épaissit le recouvrement.**

**Balayage d'$\alpha_{\mathrm{cov}}$** (3,5 km, lot `20260907_115229`). Hors de ce balayage, tous les lots à 1000 essais font tourner le CA-OK à $\alpha_{\mathrm{cov}} = 2$ : c'est le même `par.krig.window_factor` que l'AK, appliqué à la covariance de chaque mode dans `window_select`. Erreur au pas 30 et fin de vol avec les définitions de ce document.

| $\alpha_{\mathrm{cov}}$ | acq. : modes | points / fenêtre (rayon) | union | $S/U$ | pours. : union | $S/U$ | pas 30 | fin | Mahal. |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 24,3 | 27,7 (10,4 km) | 148 | 4,6 | 28,2 | 1,81 | 480 m | 382 m | 99,7 % |
| 1,5 | 24,1 | 35,8 (11,8 km) | 221 | 3,9 | 28,3 | 1,83 | 477 m | 379 m | 99,8 % |
| 2 | 24,5 | 48,9 (13,8 km) | 323 | 3,7 | 28,2 | 1,81 | 475 m | 379 m | 99,6 % |
| 3 | 24,3 | 87,4 (18,5 km) | 538 | 4,0 | 28,5 | 1,83 | 476 m | 382 m | 99,7 % |

**Balayage de la bande passante du mean-shift** (3,5 km, lot `20260908_154926`) :

| bande | acq. : pas multi | modes | points / fenêtre (rayon) | union | $S/U$ | pours. : pas multi | pas 30 | fin | Mahal. |
|---|---|---|---|---|---|---|---|---|---|
| 750 m | 42 % | 25,4 | 45,8 (13,4 km) | 295 | 3,9 | 9,5 % | 480 m | 375 m | 99,8 % |
| 1000 m | 37 % | 24,5 | 48,9 (13,8 km) | 323 | 3,7 | 2,7 % | 475 m | 379 m | 99,6 % |
| 2000 m | 32 % | 17,1 | 63,0 (15,7 km) | 361 | 3,0 | 0,1 % | 477 m | 376 m | 99,8 % |
| 4000 m | 29 % | 8,2 | 113,5 (21,0 km) | 381 | 2,5 | 0,0 % | 475 m | 383 m | 99,6 % |

- **$\alpha_{\mathrm{cov}}$ n'agit qu'à l'acquisition**, les fenêtres de poursuite étant au plancher. À $\alpha_{\mathrm{cov}} = 1$ les fenêtres par mode tombent au plancher et l'union est divisée par deux, mais chaque point sert alors à 4,6 fenêtres : **baisser $\alpha_{\mathrm{cov}}$ réduit la zone couverte, pas la réutilisation**. Celle-ci vient du rayon du plancher (9,9 km) comparé à l'écart entre 24 modes.
- **Une bande plus large fragmente moins** (8 modes à 4 km) et réduit le recouvrement ($S/U$ = 2,5), **mais les fenêtres grossissent** (114 points) et l'union se rapproche de la fenêtre de l'AK : des clusters plus gros, donc des fenêtres plus grandes.
- **La précision ne bouge dans aucun des deux balayages** (475 à 480 m au pas 30, 375 à 383 m en fin de vol, 99,6 à 99,8 %). À 3,5 km la fragmentation ne coûte rien ; H2 reste ouverte à 7,5 km, là où le CA-OK perd à l'acquisition.

### Lecture

- **Ce ne sont pas deux ou trois modes mais une vingtaine** à chaque pas multimodal d'acquisition, dans les deux lots et à toutes les mailles (72 à 87 % de ces pas en ont au moins quatre). Pour un nuage de $\sigma_0 = 8$ à 12 km, un mean-shift de 1 km de bande découpe la densité des particules en morceaux. H2 décrit bien la situation : la fragmentation est mesurée ; qu'elle coûte en précision reste à tester.
- **En poursuite, les fenêtres sont au plancher (25 points) et presque confondues** : deux modes ont environ 22 points sur 25 en commun (union 28 à 29 points pour deux fenêtres de 25).
- **À l'acquisition, chaque point retenu sert à 3 à 7 fenêtres**, et d'autant plus que la maille est lâche. Le CA-OK réajuste une vingtaine de fois des hyperparamètres sur des sous-ensembles qui se chevauchent d'une même zone.
- **Cette zone est celle de l'AK** : l'union du CA-OK fait 0,8 à 1,0 fois la fenêtre de l'AK dans les mêmes situations. La réunion, elle, touche 2,5 à 4 fois moins de points que l'AK, et passe sous le plancher à 6,5 et 7,5 km dans le lot de référence (H1).
- **Conséquence pour la comparaison** : le CA-OK actuel ne teste pas la localité des hyperparamètres par mode. Il teste « une vingtaine d'ajustements sur la zone de l'AK » contre « un ajustement sur la même zone ». Que le CA-OK ne batte pas l'AK ne dit donc rien de la non-stationnarité à l'échelle des modes.

### Imposer des fenêtres disjointes

Exigence de Bastien : les clusters ne doivent pas réutiliser les mêmes points. Les chiffres écartent trois manières de l'obtenir, et une objection de Bastien écarte une quatrième :

- **Partager les points communs entre modes** (au centre le plus proche) : écarté. En poursuite, les 28 points de l'union iraient à deux modes, 14 chacun. À l'acquisition à 7,5 km, 77 points pour 21 modes, moins de 4 chacun. Tout passe sous le plancher que le critère 2 impose.
- **Élargir la bande du mean-shift** : écarté. À 4 km il reste 8 modes qui se recouvrent ($S/U$ = 2,5), avec des fenêtres de 114 points. La bande réduit la fragmentation sans garantir la disjonction, et elle grossit les fenêtres.
- **Baisser $\alpha_{\mathrm{cov}}$ pour le CA-OK** : écarté. À 1, les fenêtres sont au plancher et chaque point sert à 4,6 fenêtres. Rien ne change en poursuite, où le plancher décide déjà.
- **Fusionner puis refenêtrer** (proposé puis retiré le 2026-09-16) : deux modes dont les fenêtres partagent un point deviennent un cluster, dont les moments donnent une nouvelle fenêtre. **Objection de Bastien : les fenêtres grossissent et on revient à l'AK.** C'est exact, et quel que soit $\alpha_{\mathrm{cov}}$. Deux modes de même poids, d'étalement $\sigma$ et distants de $d$ ont une covariance fusionnée $\sigma^2 + d^2/4$ le long de leur axe, donc un rayon d'au moins $1{,}5\,\alpha_{\mathrm{cov}}\,d$. Or leurs deux fenêtres au plancher ne portent qu'à $d/2 + r_p$ du milieu. La fenêtre fusionnée déborde dès que $d > r_p/(1{,}5\,\alpha_{\mathrm{cov}} - 0{,}5)$, soit $0{,}4\,r_p$ à $\alpha_{\mathrm{cov}} = 2$ et $r_p$ à 1, c'est-à-dire presque toujours, puisque les fenêtres se touchent jusqu'à $d = 2r_p$. À $d \approx 2r_p$, elle atteint $6\,r_p$ (ou $3\,r_p$) là où les deux fenêtres couvraient $2\,r_p$ : elle avale les modes voisins, qui fusionnent à leur tour. La cascade vient de la règle, pas du scénario, et la lecture « tout fusionne, donc pas de localités séparées » aurait été fausse. C'est le défaut même de l'AK : la covariance d'un ensemble multimodal mesure l'écart entre les modes, pas la zone qu'occupent les particules.
- **Composantes connexes, sans refenêtrer** (recommandé) : on relie deux modes dont les fenêtres partagent un point, et chaque composante connexe forme un groupe. **Son ensemble d'apprentissage est l'union de ses fenêtres, jamais recalculée**, donc il ne peut pas dépasser ce que les modes ont retenu. Un ajustement et un krigeage par groupe ; chaque particule est krigée par le groupe de son centre le plus proche. Les groupes sont disjoints par construction, puisqu'un point commun les aurait reliés, et leur total est l'union actuelle. S'il n'y a qu'un groupe, `query_ak`, comme le cas unimodal aujourd'hui : le CA-OK ne s'écarte alors de l'AK que si le nuage se sépare en localités disjointes. Le chaînage — deux modes lointains reliés par des modes intermédiaires — n'est pas un défaut : il dit que les particules occupent l'espace entre eux à l'échelle d'une fenêtre, donc que la loi a posteriori n'y a pas de trou. `n_modes` compterait les groupes, ce qui lève aussi le point de vigilance « multi n'est pas multimodal ».

**Prévision, en poursuite** : deux fenêtres au plancher ne forment deux groupes que si leur union dépasse deux fenêtres pleines, donc la part des pas multimodaux qui en gardent deux est au plus (union − 25)/25 : **13 à 20 %** dans le lot de référence, **15 à 31 %** dans le multimodal. Rapporté à tous les pas de poursuite, le CA-OK disjoint serait l'AK sur au moins **99,6 %** des pas à 2,5 et 3,5 km, 97 % à 5,5 km, 91 % à 7,5 km (référence), et 85 % à 7,5 km (multimodal). **À l'acquisition**, impossible à prévoir sans enregistrer les fenêtres : $S/U$ de 3 à 7 est compatible avec une seule chaîne comme avec quelques groupes disjoints.

**Si l'acquisition ne forme qu'un groupe**, le CA-OK disjoint sera l'AK sur ces deux scénarios, et cette fois la lecture tient, puisque les groupes ne grossissent pas : leurs nuages n'offrent jamais deux localités séparées à l'échelle d'une fenêtre (plancher et $\alpha_{\mathrm{cov}}$). Pour montrer ce que le CA-OK apporte, il faudra alors C ou D.

**À faire, pas encore décidé** :
1. **Mesurer les groupes avant de coder la variante** : une copie de `query_cak`, sous le scratchpad, qui enregistre le nombre et la taille des composantes à chaque pas sans changer ce que le filtre krige. Cinquante essais à 2,5, 5,5 et 7,5 km sur les deux lots.
2. **S'il y a des pas d'acquisition à plusieurs groupes**, l'option dans `query_cak` (ancien comportement par défaut, pour que les lots existants restent reproductibles), puis le CA-OK seul sur les deux lots à 1000 essais — mêmes graines, mêmes tirages de relevé, mailles et $\sigma_0$ relus dans leur `par` — pour rester apparié à l'AK et à la réunion déjà calculés.

## Fenêtres resserrées : $\alpha_{\mathrm{cov}} = 1{,}2$ et $\tilde{N}_{\min} = 10$ (2026-09-17, 100 essais)

Demande de Bastien : AK et CA-OK à $\alpha_{\mathrm{cov}} = 1{,}2$ et $\tilde{N}_{\min} = 10$, et le CA-OK sur plusieurs bandes. Scénario de référence, six mailles, RPF, **100 essais sur 10 tirages** : ce sont exactement les 100 premiers essais de `20260912_030842` ($\alpha_{\mathrm{cov}} = 2$, $\tilde{N}_{\min} = 25$), mêmes relevés, mêmes graines, même code. Lot `resultats/20260917_111800_reference_libre_100runs.mat` (bande 1000 m). Le CA-OK à 4000 m tourne ; celui à 16 000 m a été retiré par Bastien.

**Statistiques appariées robustes**. Les écarts sont la médiane, essai par essai, de l'erreur nouvelle moins l'ancienne. Les moyennes sont inutilisables : quelques essais qui divergent font des écarts types de plusieurs centaines de mètres. Entre parenthèses, la part des essais qui empirent et le $p$ du test des rangs signés. Fin de vol = médiane des dix derniers pas.

| maille | AK, fin de vol | CA-OK, fin de vol | Mahal. AK | Mahal. CA-OK | coût AK | coût CA-OK |
|---|---|---|---|---|---|---|
| 2,5 km | −1 m (49 %, $p$ = 0,07) | +10 m (56 %, 0,02) | 99 → 99 % | 99 → 99 % | ÷ 6,8 | ÷ 6,0 |
| 3,5 km | **+76 m** (76 %, < 0,001) | **+93 m** (78 %, < 0,001) | 99 → 99 % | 99 → 99 % | ÷ 7,4 | ÷ 4,1 |
| 4,5 km | **+123 m** (75 %, < 0,001) | **+131 m** (73 %, < 0,001) | 99 → 99 % | 99 → 99 % | ÷ 6,8 | ÷ 3,1 |
| 5,5 km | **−122 m** (37 %, 0,03) | **−137 m** (40 %, 0,05) | 89 → 90 % | 89 → 91 % | ÷ 7,4 | ÷ 2,7 |
| 6,5 km | −153 m (38 %, 0,23) | **−257 m** (32 %, 0,003) | 84 → 76 % | 79 → 77 % | ÷ 6,7 | ÷ 2,6 |
| 7,5 km | **+134 m** (59 %, 0,03) | **+271 m** (60 %, 0,001) | **88 → 66 %** | **91 → 67 %** | ÷ 5,0 | ÷ 2,6 |

- **Aux mailles fines, on perd en précision** : +76 à +131 m en fin de vol à 3,5 et 4,5 km, trois essais sur quatre moins bons. À 2,5 et 4,5 km l'erreur au pas 30 augmente aussi (+25 à +141 m). La cohérence ne bouge pas, et la NEES remonte vers 1 (0,19 → 0,27 à 3,5 km) : le filtre déclare moins, mais il est moins précis.
- **À 5,5 et 6,5 km, on gagne en fin de vol** : −122 à −257 m. Mais à 5,5 km, 13 essais sur 100 repassent au-dessus du seuil de 2000 m (89 → 76 %) tout en restant cohérents : ils finissent entre 2000 et 3500 m.
- **À 7,5 km, la cohérence s'effondre** : le Mahalanobis perd 22 à 24 points (88 → 66 % pour l'AK, 91 → 67 % pour le CA-OK). La NEES de l'AK passe de 1,71 à 2,05, et l'erreur finale augmente de 134 à 271 m. Hypothèse non testée : des hyperparamètres ajustés sur 10 points, dans une fenêtre de 13 km de rayon, rendent la variance de krigeage trop confiante.
- **Le coût baisse beaucoup** : ÷ 5 à 7 pour l'AK, ÷ 2,6 à 6 pour le CA-OK, et le temps de calcul est divisé par 5.
- **Qui est responsable ?** Les deux paramètres ont bougé ensemble. Mais à 3,5 km, le balayage d'$\alpha_{\mathrm{cov}}$ à 1000 essais ne montrait aucun effet (379 à 382 m), et celui de $\tilde{N}_{\min}$ donnait déjà +11 m à 15. **La perte aux mailles fines vient donc probablement de $\tilde{N}_{\min} = 10$.** C'est une inférence ; pour trancher, il faut $\alpha_{\mathrm{cov}} = 1{,}2$ avec $\tilde{N}_{\min} = 25$.

**AK contre CA-OK, nouveaux réglages** (mêmes essais) : aucune différence en fin de vol, à aucune maille ($p \ge 0{,}37$). À 7,5 km, **le CA-OK acquiert toujours plus mal** : +221 m au pas 30, 67 % des essais moins bons, $p$ = 0,001, contre +211 m avec les anciens réglages. Resserrer les fenêtres ne supprime pas ce défaut.

**Fenêtres par mode du CA-OK**, bande inchangée (1000 m) :

| maille | acq. : modes | points / fenêtre | union | $S/U$ | pours. : modes | union | $S/U$ |
|---|---|---|---|---|---|---|---|
| 2,5 km | 23,9 | 32,4 | 331 | 2,34 | 2,02 | 12,5 | 1,61 |
| 3,5 km | 25,3 | 19,3 | 178 | 2,75 | 2,04 | 12,3 | 1,66 |
| 4,5 km | 24,2 | 14,2 | 100 | 3,44 | 2,12 | 12,0 | 1,78 |
| 5,5 km | 20,9 | 12,2 | 64 | 4,00 | 2,41 | 12,2 | 1,98 |
| 6,5 km | 20,3 | 11,2 | 48 | 4,77 | 2,74 | 12,7 | 2,18 |
| 7,5 km | 21,2 | 10,6 | 38 | 5,94 | 3,30 | 13,0 | 2,54 |

Toujours une vingtaine de modes à l'acquisition, puisque la bande n'a pas changé. Chaque point sert à 2,3 à 5,9 fenêtres, contre 3,3 à 7,4 avant. En poursuite, deux fenêtres de 10 points en ont encore environ 7,5 en commun, contre 22 sur 25 avant. **Le recouvrement diminue, mais la réutilisation des points demeure.**

**CA-OK à bande 4 km** (lot `20260917_113709_reference_libre_100runs.mat`, mêmes réglages, mêmes 100 essais), contre la bande de 1 km :

| maille | pas multimodaux | acq. : modes | points / fenêtre | union | $S/U$ | fin, écart apparié | pas 30, écart apparié | CA-OK 4 km contre AK, pas 30 |
|---|---|---|---|---|---|---|---|---|
| 2,5 km | 7 → 6 % | 8,0 | 82 | 360 | 1,83 | +1 m ($p$ = 0,88) | 0 m | −3 m ($p$ = 0,75) |
| 3,5 km | 11 → 6 % | 7,9 | 43 | 186 | 1,82 | 0 m | 0 m | −1 m |
| 4,5 km | 15 → 8 % | 7,1 | 26 | 101 | 1,83 | −4 m | −3 m | +2 m |
| 5,5 km | 27 → 9 % | 7,1 | 19 | 71 | 1,85 | −2 m | −18 m ($p$ = 0,12) | +1 m |
| 6,5 km | 44 → 12 % | 6,6 | 15 | 51 | 1,92 | +4 m | −28 m ($p$ = 0,38) | −22 m ($p$ = 0,40) |
| 7,5 km | 58 → 18 % | 6,6 | 13 | 40 | 2,09 | +14 m ($p$ = 0,82) | **−167 m** ($p$ = 0,003) | **+100 m** ($p$ = 0,02) |

- **Quatre fois la bande, trois fois moins de modes** à l'acquisition (7 au lieu de 22), et **un point ne sert plus qu'à deux fenêtres** environ, contre 2,3 à 5,9. L'union couverte ne change presque pas.
- **En poursuite, le partage devient rare et presque disjoint** : 0 à 7 % des pas, avec $S/U$ de 1,1 à 1,4. Ce sont de vrais modes séparés, et non plus des morceaux d'un même nuage.
- **Aucun effet sur la fin de vol**, à aucune maille.
- **À 7,5 km, l'acquisition du CA-OK s'améliore de 167 m, mais il reste moins bon que l'AK de 100 m.** La fragmentation expliquait donc une partie de son défaut d'acquisition (H2), pas tout. À 100 essais, sa cohérence y est un peu plus basse (Mahalanobis 67 → 61 %, convergence 42 → 34 %), un écart de l'ordre du bruit d'échantillonnage.
- **Le CA-OK à 4 km ne bat l'AK nulle part** ; il est 1,7 à 3,6 fois plus rapide que le CA-OK à 1 km.

## Statique contre AK au scénario de référence, et pourquoi l'article gagnait plus (2026-09-17)

Remarque de Bastien : dans le manuscrit, ni l'AK ni le CA-OK ne gagnent en précision sur le statique, alors que l'article EUSIPCO les montrait nettement meilleurs. Campagne A (statique) et lot des quatre modèles (AK) : mêmes 1000 essais sur les mêmes 100 relevés, RPF, $\alpha_{\mathrm{cov}} = 2$, $\tilde{N}_{\min} = 25$.

| maille | fin de vol, statique → AK | écart apparié | RMSE finale (figure) | pas 30, statique → AK | Mahal. |
|---|---|---|---|---|---|
| 2,5 km | 237 → 239 m | +2 m ($p$ = 0,02) | +4 % | 340 → 364 m | 100 → 100 % |
| 3,5 km | 438 → 381 m | −28 m ($p$ < 0,001) | −7 % | 508 → 473 m | 100 → 100 % |
| 4,5 km | 779 → 663 m | −66 m ($p$ < 0,001) | −13 % | 623 → 598 m | 99 → 99 % |
| 5,5 km | 1119 → 1035 m | −14 m ($p$ = 0,08) | −7 % | 705 → 772 m | 91 → 90 % |
| 6,5 km | 1473 → 1357 m | −57 m ($p$ < 0,001) | −3 % | 760 → 861 m | 85 → 87 % |
| 7,5 km | 1683 → 1684 m | −42 m ($p$ = 0,001) | −1 % | 1051 → 1459 m | 81 → 87 % |

Fin de vol = médiane des dix derniers pas ; RMSE finale = celle des courbes du manuscrit, sur les essais convergés. **Le gain existe mais il est faible** : au plus 13 % sur la RMSE, invisible sur des courbes en échelle logarithmique. **À l'acquisition, l'AK est pire que le statique dès 5,5 km.** Le CA-OK suit l'AK en fin de vol et acquiert encore plus mal (+454 m au pas 30 à 7,5 km). Seul gain net : +6 points de Mahalanobis à 7,5 km. La campagne B à 3,5 km redonne exactement les mêmes chiffres.

**L'avion ne volait pas plus lentement** (≈ 500 m/s dans les deux), mais tout le reste diffère :

| | EUSIPCO (tableau I) | référence du chapitre |
|---|---|---|
| pas de temps | 0,5 s, une mesure tous les 250 m ($\alpha_{\mathrm{corr}} \approx 10$) | 5 s, une tous les 2,5 km ($\alpha_{\mathrm{corr}} \approx 1$) |
| bruit de relevé / capteur | 80 / 80 nT | 10 / 5 nT |
| $P_0$ position | 3 km | 8 km |
| vol | ligne droite, 150 km, 600 pas | tondeuse, 375 km, 150 pas |
| filtre, noyau, plafond | APF, gaussien, aucun | RPF, Matérn 3/2, 2000 |
| code | 2A et ses deux défauts (noyau √2, division de l'APF commentée) | corrigé |

Trois mécanismes candidats, aucun testé :
1. **Le pas de temps.** Une dizaine de mesures consécutives lisent la même erreur de carte et le filtre les compte comme indépendantes (`th:correlated_map_error`) : toute différence de qualité de carte est amplifiée environ dix fois.
2. **Le relevé bruité.** À 80 nT, l'erreur de carte pèse beaucoup plus lourd face au champ, et le choix des hyperparamètres avec elle.
3. **La ligne droite à travers des régions contrastées.** C'est l'explication de l'article lui-même : hyperparamètres globaux valables des itérations 75 à 200 seulement. C'est H4.

La reproduction du 2026-08-25, faite depuis le tableau I de l'article et non depuis le code 2A, a déjà montré qu'aux valeurs de l'article le nuage se partage et que le CA-OK devient le meilleur des trois krigeages : **le gain tient au scénario, pas au code**. Le scénario `eusipco` existe toujours dans `campagne_chapitre_4.m`.

**Proposé, non lancé** : statique, AK et CA-OK sur `eusipco` à 100 essais pour retrouver le gain avec le code actuel, puis un facteur à la fois vers la référence (pas de temps, bruits, trajectoire) pour savoir lequel le crée. C'est ce qui dirait dans le manuscrit **quand** la localité paie.

## Pourquoi le CA-OK perd en RMSE à l'acquisition (2026-09-22)

Question de Bastien, en lisant la figure fusionnée : le CA-OK y passe au-dessus de l'AK et de la NNAK avant le premier virage, aux mailles lâches. Lot `20260912_030842`, 1000 essais, RMSE des essais convergés, comme la figure.

| maille | méthode | $k=10$ | $k=20$ | $k=30$ | $k=50$ | $k=150$ |
|---|---|---|---|---|---|---|
| 6,5 km | AK | 7199 | 2386 | 1509 | 977 | 1004 |
| | CA-OK | 10 556 | 2971 | 3207 | 992 | 1037 |
| | NNAK | 8474 | 3230 | 2127 | 1177 | 997 |
| 7,5 km | AK | 10 855 | 4052 | 2760 | 1261 | 1090 |
| | CA-OK | 14 311 | 7747 | 6538 | 1353 | 1100 |
| | NNAK | 11 551 | 4992 | 3394 | 2020 | 1284 |

**Ce n'est pas la distribution qui se décale, c'est sa queue qui s'épaissit.** À 7,5 km, au pas 20, sur les essais convergés :

| méthode | médiane | $q_{90}$ | $q_{95}$ | RMSE | part au-delà de 5 km |
|---|---|---|---|---|---|
| AK | 1199 m | 2748 m | 4130 m | 4052 m | 4 % |
| NNAK | 1446 m | 4528 m | 11 971 m | 4992 m | 8 % |
| CA-OK | 1681 m | 11 979 m | 20 350 m | 7747 m | 19 % |

La médiane du CA-OK n'est que 40 % au-dessus de celle de l'AK ; son $q_{90}$ est quatre fois plus haut. Les écarts appariés au pas 30 le confirment : +104 m de médiane contre l'AK à 7,5 km, +29 m à 6,5 km, moins de 10 m aux mailles fines. **Tout se referme au premier virage** : dès le pas 50 les trois méthodes sont à quelques pour cent, et les erreurs finales sont identiques.

**Deux mécanismes, aucun testé.** L'ordre AK, NNAK, CA-OK les oriente : la NNAK krige sur moins de points encore — 17 à 7,5 km contre 27 par fenêtre pour le CA-OK — et reste devant lui, donc ce n'est pas le nombre de points.
1. **La vraisemblance du CA-OK est une fonction par morceaux de la position**, chaque morceau ayant ses hyperparamètres et ses points, avec des sauts aux frontières entre clusters. Un cluster qui déclare une variance plus grande est moins pénalisé pour un même écart, donc il survit au rééchantillonnage alors qu'il est faux. L'AK et la NNAK notent tout le nuage sous une seule loi.
2. **Les particules proches d'une frontière sont krigées par une fenêtre centrée ailleurs**, donc au bord de leur domaine d'interpolation. Avec neuf à vingt modes, beaucoup de particules sont dans ce cas.

### La variante à fenêtre unique, et ce qu'elle coûterait

Proposition de Bastien : garder le découpage pour **choisir** les points, réunir les fenêtres des modes, ajuster une fois et kriger une fois. C'est le cas « un seul groupe » des composantes connexes ci-dessus, et la cousine de la NNAK, avec des disques par mode et un plancher au lieu des quatre plus proches voisins par particule. Elle supprime les deux mécanismes d'un coup.

Coût par essai, convention de `tab:cout_cak`, l'union étant lue dans `n_used` :

| maille | AK | CA-OK | NNAK | CA-OK sur l'union |
|---|---|---|---|---|
| 2,5 km | 46,5 G | 2,5 G | 2,7 G | 24,9 G |
| 3,5 km | 14,6 G | 1,1 G | 1,2 G | 7,9 G |
| 4,5 km | 5,9 G | 0,9 G | 0,7 G | 3,7 G |
| 5,5 km | 3,2 G | 0,9 G | 0,5 G | 2,4 G |
| 6,5 km | 2,0 G | 0,9 G | 0,4 G | 1,9 G |
| 7,5 km | 1,6 G | 0,9 G | 0,3 G | 1,6 G |

**Elle rend l'avantage de coût du CA-OK** — deux à dix fois son prix — parce que chaque particule est résolue contre l'union et non contre sa petite fenêtre. Aux mailles lâches, l'union fait la taille de la fenêtre de l'AK (45 points contre 44 à 7,5 km) et la variante devient l'AK, au même prix, là précisément où le CA-OK perd.

**Variante plus légère** : garder les fenêtres et les résolutions par cluster, mais un seul jeu d'hyperparamètres pour tous les clusters. Les poids redeviennent comparables sans toucher au coût. Réserve : la variance annoncée dépend aussi de la géométrie des points retenus, donc il restera des écarts entre clusters, plus petits.

### Les deux lots lancés le 2026-09-22 à 11h44

Réglages du manuscrit ($\alpha_{\mathrm{cov}} = 2$, $\tilde{N}_{\min} = 25$, bande 1000 m), six mailles, RPF, 100 essais sur 10 tirages, donc appariés aux cent premiers essais de `20260912_030842`.
1. **AK et CA-OK à hyperparamètres globaux** (`refit='global'`) : si la queue disparaît, elle venait des ajustements par cluster. L'AK tourne aussi, pour séparer ce que les hyperparamètres globaux font à tout le monde de ce qu'ils font au CA-OK.
2. **CA-OK sur l'union de ses fenêtres** (`par.krig.cak_union`, ajouté à `query_cak` le même jour) : la variante de Bastien.

Pilote et journal dans le scratchpad de la session `408c2db4-…`, `go_cak_union.sh` et `log_cak_union.txt`.

### Résultat du premier lot : ce sont bien les ajustements par cluster

Lot `20260922_120645_reference_libre_100runs.mat`, 22 minutes. Écarts appariés du CA-OK à l'AK, médiane des écarts essai par essai, au pas 20 et au pas 30 :

| maille | locaux, pas 20 | globaux, pas 20 | locaux, pas 30 | globaux, pas 30 |
|---|---|---|---|---|
| 2,5 km | +1 m | −5 m | +1 m | −2 m |
| 3,5 km | **+14 m** ($p$ = 0,03) | +1 m | +3 m | +3 m |
| 4,5 km | **+34 m** ($p$ = 0,001) | −5 m | −4 m | −3 m |
| 5,5 km | −9 m | −10 m | −13 m | −12 m |
| 6,5 km | **+156 m** ($p$ < 0,001) | −2 m | +11 m | +5 m |
| 7,5 km | **+215 m** ($p$ = 0,004) | +7 m | **+211 m** ($p$ = 0,006) | +9 m |

**Avec un seul jeu d'hyperparamètres, le CA-OK ne se distingue plus de l'AK à aucune maille.** Le premier mécanisme est donc le bon : c'est l'ajustement par cluster, sur vingt-cinq à trente points, qui rend les vraisemblances incomparables d'un cluster à l'autre. La queue s'allège sans disparaître tout à fait — à 7,5 km, au pas 20, le $q_{90}$ du CA-OK passe de 15,2 à 5,6 km et la part d'essais au-delà de 5 km de 30 à 12 %, contre 3,9 km et 7 % pour l'AK global —, donc les frontières entre clusters pèsent peut-être encore un peu. La variante à fenêtre unique le dira.

### Les quatre méthodes, mêmes 100 essais (lot `20260922_124202` pour la fenêtre unique)

Erreur médiane au pas 20 et au pas 30, part d'essais au-delà de 5 km au pas 20, erreur de fin de vol, les deux taux de convergence, et le coût par essai. Écart apparié à l'AK au pas 20 entre parenthèses.

| maille | méthode | pas 20 | $q_{90}$ | > 5 km | pas 30 | fin | seuil | Mahal. | Gflop |
|---|---|---|---|---|---|---|---|---|---|
| 2,5 km | AK | 285 | 495 | 0 % | 469 | 186 | 99 % | 99 % | 46,5 |
| | NNAK | 228 (+3) | 548 | 0 % | 452 | 204 | 98 % | 98 % | 2,7 |
| | CA-OK distincts | 284 (+1) | 491 | 0 % | 469 | 189 | 99 % | 99 % | 2,5 |
| | CA-OK réunis | 288 (−0) | 472 | 0 % | 462 | 182 | 99 % | 99 % | 25,8 |
| 4,5 km | AK | 588 | 1012 | 0 % | 498 | 758 | 99 % | 99 % | 5,9 |
| | NNAK | 459 (**−106**) | 985 | 0 % | 856 | 894 | 96 % | 95 % | 0,7 |
| | CA-OK distincts | 634 (**+34**) | 1069 | 0 % | 516 | 772 | 99 % | 99 % | 0,9 |
| | CA-OK réunis | 572 (−27) | 1065 | 0 % | 477 | 774 | 99 % | 99 % | 4,1 |
| 6,5 km | AK | 1061 | 1998 | 2 % | 961 | 1639 | 64 % | 84 % | 2,0 |
| | NNAK | 1088 (+1) | 10 985 | 11 % | 977 | 1244 | 70 % | 80 % | 0,4 |
| | CA-OK distincts | 1584 (**+156**) | 4150 | 8 % | 1057 | 1712 | 59 % | 79 % | 0,9 |
| | CA-OK réunis | 1057 (+12) | 2493 | 5 % | 942 | 1624 | 66 % | 87 % | 1,9 |
| 7,5 km | AK | 1960 | 5507 | 10 % | 1645 | 2028 | 50 % | 88 % | 1,6 |
| | NNAK | 2221 (**+156**) | 15 279 | 19 % | 1651 | 2249 | 42 % | **58 %** | 0,3 |
| | CA-OK distincts | 2048 (**+215**) | 15 243 | 30 % | 2018 | 2064 | 47 % | 91 % | 0,9 |
| | CA-OK réunis | 1955 (−49) | 11 967 | 15 % | 1747 | 2031 | 48 % | 91 % | 1,7 |

- **Aux mailles fines, les quatre sont équivalentes en précision** et ne se départagent que par le coût : le CA-OK par cluster et la NNAK à 2,5 Gflop, la fenêtre unique à 25,8 et l'AK à 46,5.
- **La fenêtre unique fait bien ce qu'on attendait d'elle** : plus d'écart à l'AK à aucune maille, et la queue de l'acquisition ramenée de 30 à 15 % d'essais au-delà de 5 km à 7,5 km. Elle coûte la moitié de l'AK aux mailles fines et autant que lui aux mailles lâches.
- **Le CA-OK par cluster garde le meilleur rapport** : dix-neuf fois moins cher que l'AK à 2,5 km, aussi précis jusqu'à 4,5 km, et c'est lui qui déclare la variance la plus juste à 7,5 km (91 %). Son seul défaut est la queue d'acquisition aux mailles lâches.
- **La NNAK est la moins chère partout**, mais elle perd la cohérence à 7,5 km (58 %) et se comporte de façon erratique : meilleure que l'AK en fin de vol à 5,5 et 6,5 km, moins bonne à 4,5 et 7,5.

**La variante qui reste à essayer** réunit les deux résultats du jour : garder les fenêtres et les résolutions par cluster, donc le coût du CA-OK, mais ajuster les hyperparamètres une seule fois par pas, sur l'union, et les partager. Le lot à hyperparamètres globaux montre que cela suffit à supprimer la queue ; l'ajustement par pas, lui, garde la localité que le réglage global perd.

**Résultat de côté, et il est gros : les hyperparamètres globaux valent mieux que les locaux aux mailles lâches.** Local contre global, sur les mêmes essais :

| maille | AK, pas 30 | AK, fin de vol | convergence AK | Mahalanobis AK |
|---|---|---|---|---|
| 2,5 km | −35 m ($p$ = 0,01) | +14 m | 99 → 98 % | 99 → 98 % |
| 4,5 km | −67 m | **+175 m** ($p$ < 0,001) | 99 → 99 % | 99 → 99 % |
| 5,5 km | −29 m | **−44 m** ($p$ = 0,02) | 89 → 90 % | 89 → 83 % |
| 6,5 km | **−242 m** ($p$ < 0,001) | −29 m | **64 → 69 %** | 84 → 85 % |
| 7,5 km | −17 m | **−146 m** ($p$ = 0,02) | **50 → 72 %** | **88 → 80 %** |

Le réajustement local paie en fin de vol à 4,5 km (+175 m sans lui) et coûte ailleurs. À 7,5 km il fait perdre vingt-deux points de convergence, tout en gagnant huit points de Mahalanobis : **il achète de la cohérence au prix de la précision**, ce qui va dans le sens de H3. Le CA-OK suit le même profil, en plus marqué.

## Rendre la multimodalité et l'instationnarité visibles

Par ordre de coût, du moins cher au plus démonstratif. Les commandes suivent l'idiome de `notes_campagne_E.md` et **n'ont pas été lancées** ; elles passent après les lots prioritaires de la file.

### A. Mesurer d'abord, sans filtre

1. **Carte des hyperparamètres locaux.** Ajuster $\ell$ et $\sigma_f$ par maximum de vraisemblance sur des fenêtres de 25, 50 et 250 points centrées tous les quelques kilomètres le long du corridor, puis sur toute la carte. Répéter sur plusieurs tirages de relevé pour séparer le bruit d'estimation de la variation réelle, et tracer le variogramme de $\log\ell$ : sa portée **est** $L_\theta$. Tranche H4, et dit si le corridor du scénario est simplement trop homogène. Script à écrire, une cinquantaine de lignes, quelques minutes de calcul.
2. **Distances entre modes pendant l'acquisition.** Enregistrer les centres du mean-shift sur quelques essais pour comparer l'écart entre modes à $L_\theta$ et au rayon des fenêtres.

### B. Isoler les mécanismes sur le scénario existant

1. **Hyperparamètres globaux contre locaux, AK et CA-OK** (tranche H3 et la part de la localité) — les résultats `fenetre` existent déjà dans le lot multimodal, sur les mêmes graines :
   ```bash
   matlab -batch "cd('D:/bhubert/Desktop/Notes/MATLAB'); over.campagne='libre'; over.compare.models={'ak','cak'}; over.compare.filters={'RPF'}; over.sigma_0=[12000 12000 2 2]; over.krig.refit='global'; over.balayage=struct('champ','krig.pas','valeurs',[4500 5500 7500]); over.mc.n_runs=1000; over.krig.n_surveys=100; campagne_chapitre_4"
   ```
2. **La réunion au-dessus du plancher** (tranche H1) : si le Mahalanobis remonte vers 87 % quand l'union dépasse 25 points, c'était le nombre.
   ```bash
   matlab -batch "cd('D:/bhubert/Desktop/Notes/MATLAB'); over.campagne='libre'; over.compare.models={'ak_boules'}; over.compare.filters={'RPF'}; over.krig.pas=7500; over.balayage=struct('champ','krig.k_voisins','valeurs',[4 8 16 32]); over.mc.n_runs=1000; over.krig.n_surveys=100; campagne_chapitre_4"
   ```
3. **La fragmentation du mean-shift** (tranche H2) : si l'erreur d'acquisition du CA-OK rejoint l'AK quand la bande passante grandit, c'était la fragmentation. **Déjà fait à 3,5 km** (lot `20260908_154926`, 750 à 4000 m) : précision plate, voir plus haut. Reste la maille de 7,5 km, où le CA-OK perd à l'acquisition.
   ```bash
   matlab -batch "cd('D:/bhubert/Desktop/Notes/MATLAB'); over.campagne='libre'; over.compare.models={'cak'}; over.compare.filters={'RPF'}; over.krig.pas=7500; over.balayage=struct('champ','krig.bandwidth','valeurs',[1000 2000 4000 8000]); over.mc.n_runs=1000; over.krig.n_surveys=100; campagne_chapitre_4"
   ```

### C. Faire durer la multimodalité

Le premier virage lève l'ambiguïté dès le pas 33 : le multimodal ne dure pas assez pour que le CA-OK puisse s'exprimer. Deux leviers, qui existent déjà :

1. **Le scénario `acquisition`** — vol en ligne, $\sigma_0 = 15$ km, 60 pas — n'a pas de virage pour lever l'ambiguïté. À vérifier qu'il tourne encore avec le code actuel.
   ```bash
   matlab -batch "cd('D:/bhubert/Desktop/Notes/MATLAB'); over.scenario='acquisition'; over.campagne='libre'; over.compare.models={'ok','ak','cak','ak_boules'}; over.compare.filters={'RPF'}; over.mc.n_runs=1000; over.krig.n_surveys=100; campagne_chapitre_4"
   ```
2. **Un $\sigma_0$ plus grand**, 20 km, sur la tondeuse — attention à la marge au bord, commentée dans le scénario `reference`.

Et **mesurer sur la phase multimodale seule** : erreur pendant les pas multimodaux, pas de levée d'ambiguïté, part des essais qui se verrouillent sur un mauvais mode. Les métriques de fin de vol la diluent.

### D. L'expérience démonstrative : une instationnarité contrôlée

La seule qui rendrait la localité **visible à coup sûr** : un champ synthétique dont on connaît les hyperparamètres en tout point. Par exemple deux tirages Matérn 3/2 stationnaires très différents — $\ell$ de 2 et 10 km, $\sigma_f$ de 50 et 250 nT — fondus le long d'une frontière par un masque lisse, ou un noyau non stationnaire de Paciorek–Schervish avec $\ell(x)$ variant sur une échelle $L_\theta$ réglable. Placer l'incertitude initiale à cheval sur la frontière, pour que les modes tombent dans des régimes différents, puis balayer $L_\theta$ contre l'écart entre modes et le rayon des fenêtres.

Prédiction : statique < AK < CA-OK quand $L_\theta$ est de l'ordre de l'écart entre modes ; les trois confondus quand $L_\theta$ le dépasse de beaucoup. C'est la validation directe de $\ell \ll r \ll L_\theta$. Il faut une carte au format de `carte_magnetometrie_anomalie_mexique.mat` (champs `mesures` et `pas`) et pointer `par.map.file` dessus. Variante sans synthèse : choisir un corridor réel qui traverse des régions contrastées, repérées grâce à la carte de A1.

### E. Figures pour le manuscrit

- La carte de $\ell(x)$ et $\sigma_f(x)$ avec la trajectoire par-dessus : l'instationnarité d'un coup d'œil.
- Un instantané du nuage à un pas multimodal, avec le disque unique de l'AK et les disques par mode du CA-OK sur la carte : pourquoi la fenêtre de l'AK enfle, ce que fait le CA-OK. Demande d'enregistrer les particules d'un essai.
- Les performances contre $L_\theta$, si D tourne.

## Points de vigilance

- **Appariement** : les comparaisons de B supposent les mêmes graines et les mêmes tirages de relevé que les lots existants — ne toucher ni à `run.seed` ni à `n_surveys`.
- **« multi » n'est pas « multimodal »** : c'est ce que rend le mean-shift à 1000 m de bande passante — une vingtaine de « modes » par pas d'acquisition, qui se partagent leurs points. Tant que H2 n'est pas tranchée, ne pas l'appeler part de pas multimodaux dans le manuscrit.
- **Noms des modèles** : `OK adaptatif`, `CA-OK`, `OK par réunion`, `OK éclairci`, `OK carte totale`, `carte connue`. Toujours apparier sur des motifs exclusifs.
- **Le temps de calcul n'est pas une contrainte** : 1000 essais d'emblée.
